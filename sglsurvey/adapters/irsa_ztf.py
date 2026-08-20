"""IRSA ZTF single-epoch (science-image) archive adapter.

Discovery runs against the IBE metadata search for ``ztf/products/sci``
(one row per CCD-quadrant exposure); products come from the IBE data
tree. Facts verified before writing this module are recorded in
surveys/ztf/notes/irsa_recon.md (ZSDS Explanatory Supplement section
numbers cited there).
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

IBE_SEARCH = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci"
IBE_DATA = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"
COLLECTION = "ztf/products/sci"

FILTER_WAVELENGTH_UM = {"zg": 0.472, "zr": 0.634, "zi": 0.789}
#: ZSDS §13.4 note 2: metadata INFOBITS bit 25 = archival bad-quality flag.
BAD_QUALITY_BIT = 1 << 25
#: ZSDS §10.3 "uncontaminated" template: every bit except 1 and 11.
MASK_FATAL_TEMPLATE = 6141

_PRODUCT_SUFFIX = {
    "sci": "sciimg.fits",
    "msk": "mskimg.fits",
    "diff": "scimrefdiffimg.fits.fz",
    "diffpsf": "diffimgpsf.fits",
    "psfcat": "psfcat.fits",
    "sexcat": "sexcat.fits",
}

JD_MJD = 2400000.5


def _f(row: dict, key: str) -> float | None:
    v = row.get(key, "")
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


class ZtfNominalFootprint:
    """Containment test from the inventory linear WCS (no distortion,
    no masks) — coarse stage only. Quadrant is 3072 x 3080 pixels."""

    NX, NY = 3072, 3080

    def __init__(self, obs: Observation, pad_arcsec: float = 0.0):
        w = obs.wcs
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [w["crpix1"], w["crpix2"]]
        wcs.wcs.crval = [w["crval1"], w["crval2"]]
        wcs.wcs.cd = [[w["cd1_1"], w["cd1_2"]], [w["cd2_1"], w["cd2_2"]]]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        self._wcs = wcs
        self._nx, self._ny = self.NX, self.NY
        # TPV distortion omitted here: worst-case TAN-vs-TPV error at the
        # quadrant corner is a few arcsec, so the coarse pad must exceed it.
        self._pad_pix = pad_arcsec / (abs(w["cd1_1"]) * 3600.0)

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        x, y = self._wcs.wcs_world2pix([[ra_deg, dec_deg]], 0)[0]
        p = self._pad_pix
        return (-0.5 - p <= x <= self._nx - 0.5 + p
                and -0.5 - p <= y <= self._ny - 0.5 + p)

    def contains_any(self, radec_deg) -> bool:
        import numpy as np

        pix = self._wcs.wcs_world2pix(radec_deg, 0)
        with np.errstate(invalid="ignore"):
            x, y = pix[:, 0], pix[:, 1]
            p = self._pad_pix
            ok = ((x >= -0.5 - p) & (x <= self._nx - 0.5 + p)
                  & (y >= -0.5 - p) & (y <= self._ny - 0.5 + p))
        return bool(np.any(ok))


class ZtfExactFootprint:
    """Exact usable-pixel test from a mskimg product or cutout: full TPV
    WCS from the FITS header plus the ZSDS mask bitplanes (§10.3).

    Works on a cutout: the IBE cutout service preserves the full WCS
    with shifted CRPIX, so pixels outside the cutout are reported as
    out-of-bounds (not unusable) — callers must size cutouts to enclose
    the locus, and treat out-of-bounds as "not evaluated".
    """

    FATAL_MASK = MASK_FATAL_TEMPLATE

    def __init__(self, msk_path: Path):
        import warnings

        from astropy.io import fits
        from astropy.utils.exceptions import AstropyWarning

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            with fits.open(msk_path) as hdul:
                hdu = next(h for h in hdul if h.data is not None)
                self._mask = hdu.data
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


class ZtfSciAdapter:
    archive_id = "irsa-ztf"
    collection = COLLECTION
    # A quadrant is ~52' on a side; the IBE box search already uses the
    # footprint (OVERLAPS), so no half-diagonal inflation is needed.

    def __init__(self, session: requests.Session | None = None,
                 timeout_s: float = 900.0, public_only: bool = True):
        self._session = session or requests.Session()
        self._timeout = timeout_s
        self._public_only = public_only

    # -- discovery ---------------------------------------------------------
    def discovery_params(self, region: ConeRegion, time_range: MjdRange,
                         ) -> dict[str, str]:
        jd0 = time_range.start_mjd_utc + JD_MJD
        jd1 = time_range.stop_mjd_utc + JD_MJD
        where = [f"obsjd>={jd0:.5f}", f"obsjd<={jd1:.5f}",
                 "imgtype='object'"]
        if self._public_only:
            where.append("ipac_gid=1")
        return {
            "POS": f"{region.ra_deg:.6f},{region.dec_deg:.6f}",
            "SIZE": f"{2.0 * region.radius_deg:.6f}",
            "INTERSECT": "OVERLAPS",
            "WHERE": " AND ".join(where),
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
    def product_base(row_like: dict) -> str:
        ffd = str(row_like["filefracday"])
        return (f"{IBE_DATA}/{ffd[:4]}/{ffd[4:8]}/{ffd[8:]}/"
                f"ztf_{ffd}_{int(row_like['field']):06d}_"
                f"{row_like['filtercode']}_c{int(row_like['ccdid']):02d}_o_"
                f"q{int(row_like['qid'])}_")

    def _observation(self, row: dict, snapshot_id: str,
                     discovered_utc: str) -> Observation:
        jd_start = float(row["obsjd"])
        exptime = _f(row, "exptime") or 30.0
        mjd_start = jd_start - JD_MJD
        half = exptime / 2.0 / 86400.0
        corners = tuple(
            (float(row[f"ra{i}"]), float(row[f"dec{i}"]))
            for i in (1, 2, 3, 4))
        wcs = {"naxis1": ZtfNominalFootprint.NX,
               "naxis2": ZtfNominalFootprint.NY,
               "crpix1": _f(row, "crpix1"), "crpix2": _f(row, "crpix2"),
               "crval1": _f(row, "crval1"), "crval2": _f(row, "crval2"),
               "ctype1": "RA---TPV", "ctype2": "DEC--TPV",
               "cd1_1": _f(row, "cd11"), "cd1_2": _f(row, "cd12"),
               "cd2_1": _f(row, "cd21"), "cd2_2": _f(row, "cd22")}
        base = self.product_base(row)
        products = {kind: {"url": base + suffix}
                    for kind, suffix in _PRODUCT_SUFFIX.items()}
        infobits = int(_f(row, "infobits") or 0)
        return Observation.build(
            archive_id=self.archive_id, collection=self.collection,
            release=f"ipac_gid{row.get('ipac_gid', '')}",
            native_key={"pid": int(row["pid"]),
                        "filefracday": int(row["filefracday"]),
                        "field": int(row["field"]),
                        "ccdid": int(row["ccdid"]),
                        "qid": int(row["qid"]),
                        "filtercode": row["filtercode"]},
            band=row["filtercode"],
            wavelength_um=FILTER_WAVELENGTH_UM.get(row["filtercode"]),
            t_start_mjd_utc=mjd_start, t_mid_mjd_utc=mjd_start + half,
            t_stop_mjd_utc=mjd_start + 2 * half, exptime_s=exptime,
            corners_icrs_deg=corners, wcs=wcs,
            quality_flags={
                "infobits": infobits,
                "bad_quality": bool(infobits & BAD_QUALITY_BIT),
                "seeing": _f(row, "seeing"), "airmass": _f(row, "airmass"),
                "moonillf": _f(row, "moonillf"),
                "maglimit": _f(row, "maglimit")},
            products=products,
            extra={"rcid": int(row["rcid"]), "fid": int(row["fid"]),
                   "nid": int(row["nid"]), "expid": int(row["expid"]),
                   "obsdate": row.get("obsdate"),
                   "ipac_pub_date": row.get("ipac_pub_date")},
            snapshot_id=snapshot_id, discovered_utc=discovered_utc,
        )

    # -- footprints --------------------------------------------------------
    def nominal_footprint(self, obs: Observation,
                          pad_arcsec: float = 0.0) -> ZtfNominalFootprint:
        return ZtfNominalFootprint(obs, pad_arcsec)

    def exact_footprint(self, obs: Observation,
                        products: ProductSet) -> ZtfExactFootprint:
        msk = next(p for p in products.products if p.kind == "msk")
        return ZtfExactFootprint(msk.path)

    # -- fetch -------------------------------------------------------------
    def fetch(self, obs: Observation, kinds: Sequence[str], dest: Path,
              *, cutout: CutoutSpec | None = None) -> ProductSet:
        """Download products or cutouts. ZTF publishes no checksums, so
        the recorded checksum is the sha256 of the received bytes. A
        404 (product not served — see recon note) raises
        ``FileNotFoundError`` so callers can record availability."""
        dest.mkdir(parents=True, exist_ok=True)
        products = []
        for kind in kinds:
            url = obs.products[kind]["url"]
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
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            products.append(LocalProduct(kind=kind, path=path,
                                         checksum=f"sha256:{digest}"))
        return ProductSet(observation_id=obs.observation_id,
                          products=tuple(products), cutout=cutout)
