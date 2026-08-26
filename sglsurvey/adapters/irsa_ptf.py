"""IRSA PTF level-1 (epochal science image) archive adapter.

Discovery runs against the IBE metadata search for ``ptf/images/level1``
(one row per CCD exposure); products come from the IBE data tree with
``pfilename``/``afilenameN`` served verbatim. Facts verified before
writing this module are recorded in
surveys/ptf-crossings/notes/ptf_recon_2026-08-26.md; mask bits are
Laher et al. 2014 Table 15.
"""

from __future__ import annotations

import csv
import hashlib
import io
from pathlib import Path
from typing import Iterator, Sequence

import requests
from astropy.wcs import WCS

from sglsurvey.adapters.base import (ConeRegion, CutoutSpec, LocalProduct,
                                     MjdRange, ProductSet)
from sglsurvey.records import Observation
from sglsurvey.snapshots import SnapshotStore

IBE_SEARCH = "https://irsa.ipac.caltech.edu/ibe/search/ptf/images/level1"
IBE_DATA = "https://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1"
COLLECTION = "ptf/images/level1"

#: fid 1 = g, 2 = R; Halpha fids are excluded by the frozen fid <= 2 gate.
FID_BAND = {1: "g", 2: "R"}
FILTER_WAVELENGTH_UM = {"g": 0.472, "R": 0.658}
#: Laher et al. 2014 Table 15: every bit except 2^1 (object detected).
MASK_FATAL_TEMPLATE = 65533

#: The IBE default column set omits fid, WCS, exptime, pid, version,
#: checksums and ancillary slots — request everything explicitly.
DISCOVERY_COLUMNS = (
    "expid,obsmjd,obsdate,exptime,aexptime,nid,fid,filter,ptffield,"
    "imgtype,pid,version,ccdid,pfilename,pchecksum,photcalflag,infobits,"
    "seeing,airmass,moonillf,crpix1,crpix2,ra,dec,cd1_1,cd1_2,cd2_1,"
    "cd2_2,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4,"
    "anciltype1,afilename1,achecksum1,anciltype2,afilename2,achecksum2,"
    "anciltype3,afilename3,achecksum3,anciltype4,afilename4,achecksum4,"
    "ipac_gid,ipac_pub_date")


def _f(row: dict, key: str) -> float | None:
    v = row.get(key, "")
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


class PtfNominalFootprint:
    """Containment test from the inventory linear WCS (metadata ra/dec
    columns are CRVAL at crpix; PV distortion omitted) — coarse stage
    only. CCD is 2048 x 4096 pixels at 1.01 arcsec/pix."""

    NX, NY = 2048, 4096

    def __init__(self, obs: Observation, pad_arcsec: float = 0.0):
        w = obs.wcs
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [w["crpix1"], w["crpix2"]]
        wcs.wcs.crval = [w["crval1"], w["crval2"]]
        wcs.wcs.cd = [[w["cd1_1"], w["cd1_2"]], [w["cd2_1"], w["cd2_2"]]]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        self._wcs = wcs
        self._nx, self._ny = self.NX, self.NY
        # TAN-vs-TPV error at the CCD corner is arcsec-scale; the coarse
        # pad must exceed it when edge behavior matters.
        self._pad_pix = pad_arcsec / (abs(w["cd1_1"]) * 3600.0)

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        x, y = self._wcs.wcs_world2pix([[ra_deg, dec_deg]], 0)[0]
        p = self._pad_pix
        return (-0.5 - p <= x <= self._nx - 0.5 + p
                and -0.5 - p <= y <= self._ny - 0.5 + p)


class PtfExactFootprint:
    """Exact usable-pixel test from a dmask product or cutout: full
    TAN-SIP WCS from the FITS header plus the Laher Table 15 bitplanes.

    Works on a cutout: the IBE cutout service preserves the full WCS
    with shifted CRPIX (verified at recon), so pixels outside the
    cutout are reported out-of-bounds (not unusable) — callers size
    cutouts to enclose the locus and treat out-of-bounds as "not
    evaluated".
    """

    FATAL_MASK = MASK_FATAL_TEMPLATE

    def __init__(self, msk_path: Path):
        import warnings

        from astropy.io import fits
        from astropy.utils.exceptions import AstropyWarning

        import numpy as np

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            with fits.open(msk_path) as hdul:
                hdu = next(h for h in hdul if h.data is not None)
                # dmask is signed int16 on disk; reinterpret as the
                # 16-bit flag word so the fatal template can be ANDed.
                self._mask = np.asarray(hdu.data).astype(np.int32) & 0xFFFF
                self._wcs = WCS(hdu.header)
                self.header = dict(hdu.header)
        self._ny, self._nx = self._mask.shape

    def status(self, radec_deg):
        """Per-point (in_bounds, usable) boolean arrays for (N,2) deg."""
        import numpy as np

        pix = self._wcs.wcs_world2pix(radec_deg, 0)
        with np.errstate(invalid="ignore"):
            x = np.rint(pix[:, 0]).astype(int)
            y = np.rint(pix[:, 1]).astype(int)
            inb = (x >= 0) & (x < self._nx) & (y >= 0) & (y < self._ny)
        usable = np.zeros(len(radec_deg), dtype=bool)
        if inb.any():
            usable[inb] = (self._mask[y[inb], x[inb]]
                           & self.FATAL_MASK) == 0
        return inb, usable

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        import numpy as np

        _, usable = self.status(np.array([[ra_deg, dec_deg]]))
        return bool(usable[0])


class PtfLevel1Adapter:
    archive_id = "irsa-ptf"
    collection = COLLECTION

    def __init__(self, session: requests.Session | None = None,
                 timeout_s: float = 900.0):
        self._session = session or requests.Session()
        self._timeout = timeout_s

    # -- discovery ---------------------------------------------------------
    def discovery_params(self, region: ConeRegion, time_range: MjdRange,
                         ) -> dict[str, str]:
        where = [f"obsmjd>={time_range.start_mjd_utc:.5f}",
                 f"obsmjd<={time_range.stop_mjd_utc:.5f}",
                 "imgtype='object'", "fid<=2"]
        return {
            "POS": f"{region.ra_deg:.6f},{region.dec_deg:.6f}",
            "SIZE": f"{2.0 * region.radius_deg:.6f}",
            "INTERSECT": "OVERLAPS",
            "WHERE": " AND ".join(where),
            "columns": DISCOVERY_COLUMNS,
            "ct": "csv",
        }

    def discover(self, region: ConeRegion, time_range: MjdRange,
                 store: SnapshotStore) -> Iterator[Observation]:
        from datetime import datetime, timezone

        params = self.discovery_params(region, time_range)
        request_utc = datetime.now(timezone.utc).isoformat()
        resp = self._session.get(IBE_SEARCH, params=params,
                                 timeout=self._timeout)
        resp.raise_for_status()
        rows = list(csv.DictReader(io.StringIO(resp.text)))
        query = "&".join(f"{k}={v}" for k, v in params.items())
        snapshot = store.store(service_url=IBE_SEARCH, query=query,
                               request_utc=request_utc,
                               response_bytes=resp.content,
                               row_count=len(rows),
                               http_status=resp.status_code)
        for row in rows:
            yield self._observation(row, snapshot.snapshot_id, request_utc)

    @staticmethod
    def _ancillaries(row: dict) -> dict[str, dict]:
        """Ancillary products keyed by anciltype — slot order varies
        between exposures (recon gotcha), so never index by slot."""
        out = {}
        for i in (1, 2, 3, 4):
            kind = (row.get(f"anciltype{i}") or "").strip()
            name = (row.get(f"afilename{i}") or "").strip()
            if kind and name:
                out[kind] = {"url": f"{IBE_DATA}/{name}",
                             "checksum_md5":
                                 (row.get(f"achecksum{i}") or "").strip()}
        return out

    def _observation(self, row: dict, snapshot_id: str,
                     discovered_utc: str) -> Observation:
        mjd_start = float(row["obsmjd"])   # shutter-open UTC (verified)
        exptime = _f(row, "aexptime") or _f(row, "exptime") or 60.0
        half = exptime / 2.0 / 86400.0
        fid = int(row["fid"])
        band = FID_BAND.get(fid, row.get("filter", "").strip())
        corners = tuple(
            (float(row[f"ra{i}"]), float(row[f"dec{i}"]))
            for i in (1, 2, 3, 4))
        wcs = {"naxis1": PtfNominalFootprint.NX,
               "naxis2": PtfNominalFootprint.NY,
               "crpix1": _f(row, "crpix1"), "crpix2": _f(row, "crpix2"),
               "crval1": _f(row, "ra"), "crval2": _f(row, "dec"),
               "ctype1": "RA---TPV", "ctype2": "DEC--TPV",
               "cd1_1": _f(row, "cd1_1"), "cd1_2": _f(row, "cd1_2"),
               "cd2_1": _f(row, "cd2_1"), "cd2_2": _f(row, "cd2_2")}
        products = {"sci": {"url": f"{IBE_DATA}/{row['pfilename']}",
                            "checksum_md5":
                                (row.get("pchecksum") or "").strip()}}
        anc = self._ancillaries(row)
        if "dmask" in anc:
            products["msk"] = anc["dmask"]
        if "sexcat" in anc:
            products["sexcat"] = anc["sexcat"]
        return Observation.build(
            archive_id=self.archive_id, collection=self.collection,
            release=f"v{row.get('version', '')}",
            native_key={"pid": int(row["pid"]),
                        "expid": int(row["expid"]),
                        "ccdid": int(row["ccdid"]),
                        "ptffield": int(row["ptffield"]),
                        "fid": fid},
            band=band,
            wavelength_um=FILTER_WAVELENGTH_UM.get(band),
            t_start_mjd_utc=mjd_start, t_mid_mjd_utc=mjd_start + half,
            t_stop_mjd_utc=mjd_start + 2 * half, exptime_s=exptime,
            corners_icrs_deg=corners, wcs=wcs,
            quality_flags={
                "infobits": int(_f(row, "infobits") or 0),
                "photcalflag": int(_f(row, "photcalflag") or 0),
                "seeing": _f(row, "seeing"), "airmass": _f(row, "airmass"),
                "moonillf": _f(row, "moonillf")},
            products=products,
            extra={"nid": int(_f(row, "nid") or 0),
                   "obsdate": row.get("obsdate"),
                   "pfilename": row.get("pfilename"),
                   "ipac_gid": row.get("ipac_gid"),
                   "ipac_pub_date": row.get("ipac_pub_date")},
            snapshot_id=snapshot_id, discovered_utc=discovered_utc,
        )

    # -- footprints --------------------------------------------------------
    def nominal_footprint(self, obs: Observation,
                          pad_arcsec: float = 0.0) -> PtfNominalFootprint:
        return PtfNominalFootprint(obs, pad_arcsec)

    def exact_footprint(self, obs: Observation,
                        products: ProductSet) -> PtfExactFootprint:
        msk = next(p for p in products.products if p.kind == "msk")
        return PtfExactFootprint(msk.path)

    # -- fetch -------------------------------------------------------------
    def fetch(self, obs: Observation, kinds: Sequence[str], dest: Path,
              *, cutout: CutoutSpec | None = None) -> ProductSet:
        """Download products or cutouts. Full-file downloads are
        verified against the archive-published MD5 (``pchecksum`` /
        ``achecksumN``); cutouts are re-encoded server-side, so their
        recorded checksum is the sha256 of the received bytes. A 404
        raises ``FileNotFoundError`` so callers can record
        availability."""
        dest.mkdir(parents=True, exist_ok=True)
        products = []
        for kind in kinds:
            meta = obs.products[kind]
            url = meta["url"]
            params = {}
            if cutout is not None:
                params = {"center": f"{cutout.ra_deg:.6f},{cutout.dec_deg:.6f}",
                          "size": f"{cutout.size_pix}pix"}
            name = url.rsplit("/", 1)[-1]
            if cutout is not None:
                name = f"cut{cutout.size_pix}-{name}"
            path = dest / name
            if not path.exists():
                resp = self._session.get(url, params=params,
                                         timeout=self._timeout)
                if resp.status_code == 404:
                    raise FileNotFoundError(url)
                resp.raise_for_status()
                path.write_bytes(resp.content)
            data = path.read_bytes()
            if cutout is None and meta.get("checksum_md5"):
                # Some achecksumN values are truncated by the archive
                # (observed: 11-char sexcat checksums) — verify the
                # published value as a prefix of the computed MD5.
                got = hashlib.md5(data).hexdigest()
                if not got.startswith(meta["checksum_md5"]):
                    raise IOError(f"MD5 mismatch for {url}: "
                                  f"{got} != {meta['checksum_md5']}")
                checksum = f"md5:{got}"
            else:
                checksum = f"sha256:{hashlib.sha256(data).hexdigest()}"
            products.append(LocalProduct(kind=kind, path=path,
                                         checksum=checksum))
        return ProductSet(observation_id=obs.observation_id,
                          products=tuple(products), cutout=cutout)
