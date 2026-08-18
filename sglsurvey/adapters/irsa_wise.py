"""IRSA WISE merged-L1b archive adapter.

Discovery runs against ``wise.neowiser_merge_p1bm_frm`` (verified union
of all four mission phases, 62.7M frame-band rows); image products come
from the unified IBE ``merge/merge_p1bm_frm`` tree (verified to host
cryo and NEOWISE frames alike). See surveys/wise/notes/irsa_recon.md.
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

TAP_SYNC = "https://irsa.ipac.caltech.edu/TAP/sync"
IBE_BASE = ("https://irsa.ipac.caltech.edu/ibe/data/wise/merge/"
            "merge_p1bm_frm")
TABLE = "wise.neowiser_merge_p1bm_frm"

BAND_WAVELENGTH_UM = {1: 3.4, 2: 4.6, 3: 12.0, 4: 22.0}

DISCOVERY_COLUMNS = [
    "scan_id", "scangrp", "frame_num", "band", "wrelease", "image_set",
    "mjd_obs", "date_obs", "exptime",
    "ra1", "dec1", "ra2", "dec2", "ra3", "dec3", "ra4", "dec4",
    "naxis1", "naxis2", "crpix1", "crpix2", "crval1", "crval2",
    "ctype1", "ctype2", "cd1_1", "cd1_2", "cd2_1", "cd2_2",
    "pxscal1", "pxscal2", "pa",
    "qual_frame", "qual_scan", "qa_status", "moon_sep", "saa_sep",
    "magzp", "magzpunc",
]

_PRODUCT_SUFFIX = {"int": "-int-1b.fits", "msk": "-msk-1b.fits.gz",
                   "unc": "-unc-1b.fits.gz"}


def _f(row: dict, key: str) -> float | None:
    v = row.get(key, "")
    return float(v) if v not in ("", None, "null") else None


class WiseNominalFootprint:
    """Containment test from the inventory WCS (nominal: no masks)."""

    def __init__(self, obs: Observation, pad_arcsec: float = 0.0):
        w = obs.wcs
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [w["crpix1"], w["crpix2"]]
        wcs.wcs.crval = [w["crval1"], w["crval2"]]
        wcs.wcs.cd = [[w["cd1_1"], w["cd1_2"]], [w["cd2_1"], w["cd2_2"]]]
        wcs.wcs.ctype = [str(w["ctype1"])[:8], str(w["ctype2"])[:8]]
        self._wcs = wcs
        self._nx, self._ny = w["naxis1"], w["naxis2"]
        self._pad_pix = pad_arcsec / abs(w["pxscal1"])

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        x, y = self._wcs.wcs_world2pix([[ra_deg, dec_deg]], 0)[0]
        p = self._pad_pix
        return (-0.5 - p <= x <= self._nx - 0.5 + p
                and -0.5 - p <= y <= self._ny - 0.5 + p)

    def contains_any(self, radec_deg) -> bool:
        """Vectorized test: True if any (ra, dec) row lands in-frame."""
        import numpy as np

        pix = self._wcs.wcs_world2pix(radec_deg, 0)
        with np.errstate(invalid="ignore"):
            x, y = pix[:, 0], pix[:, 1]
            p = self._pad_pix
            ok = ((x >= -0.5 - p) & (x <= self._nx - 0.5 + p)
                  & (y >= -0.5 - p) & (y <= self._ny - 0.5 + p))
        return bool(np.any(ok))


class WiseExactFootprint:
    """Exact usable-pixel test from a -msk product: full FITS WCS
    (including SIP distortion when present) plus the mask bitplanes.

    Fatal bits follow the All-Sky Explanatory Supplement IV.4.a Table 2
    NaN set — {0,1,2,3,4,9,10..18}: broken/noisy/dead hardware and
    sample-read saturation — plus the transient-condition bits
    {21 dynamic bad pixel, 27 cosmic ray/outlier, 28 spike-outlier},
    since a transient at the predicted position corrupts that frame's
    photometry even though the pixel is healthy.
    """

    FATAL_BITS = (0, 1, 2, 3, 4, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                  21, 27, 28)
    FATAL_MASK = sum(1 << b for b in FATAL_BITS)

    def __init__(self, msk_path: Path):
        import warnings

        from astropy.io import fits
        from astropy.utils.exceptions import AstropyWarning

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            with fits.open(msk_path) as hdul:
                self._mask = hdul[0].data
                self._wcs = WCS(hdul[0].header)
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


class WiseMergeL1bAdapter:
    archive_id = "irsa-wise"
    collection = TABLE
    # 47' x 47' frames: half-diagonal, used to convert "locus inside
    # frame" into a frame-center cone constraint.
    FRAME_HALF_DIAGONAL_DEG = 0.56

    def __init__(self, session: requests.Session | None = None,
                 timeout_s: float = 900.0):
        self._session = session or requests.Session()
        self._timeout = timeout_s

    def discovery_query(self, region: ConeRegion, time_range: MjdRange,
                        ) -> str:
        radius = region.radius_deg + self.FRAME_HALF_DIAGONAL_DEG
        cols = ", ".join(DISCOVERY_COLUMNS)
        return (
            f"SELECT {cols} FROM {TABLE} WHERE "
            f"CONTAINS(POINT('ICRS',crval1,crval2),"
            f"CIRCLE('ICRS',{region.ra_deg:.6f},{region.dec_deg:.6f},"
            f"{radius:.6f}))=1 "
            f"AND mjd_obs >= {time_range.start_mjd_utc:.5f} "
            f"AND mjd_obs <= {time_range.stop_mjd_utc:.5f}"
        )

    def discover(self, region: ConeRegion, time_range: MjdRange,
                 store: SnapshotStore) -> Iterator[Observation]:
        from datetime import datetime, timezone

        query = self.discovery_query(region, time_range)
        request_utc = datetime.now(timezone.utc).isoformat()
        resp = self._session.get(TAP_SYNC,
                                 params={"QUERY": query, "FORMAT": "CSV"},
                                 timeout=self._timeout)
        resp.raise_for_status()
        rows = list(csv.DictReader(io.StringIO(resp.text)))
        snapshot = store.store(service_url=TAP_SYNC, query=query,
                               request_utc=request_utc,
                               response_bytes=resp.content,
                               row_count=len(rows),
                               http_status=resp.status_code)
        for row in rows:
            yield self._observation(row, snapshot.snapshot_id, request_utc)

    def _observation(self, row: dict, snapshot_id: str,
                     discovered_utc: str) -> Observation:
        band_num = int(row["band"])
        mjd = float(row["mjd_obs"])
        exptime = _f(row, "exptime") or 0.0
        half = exptime / 2.0 / 86400.0
        corners = tuple(
            (float(row[f"ra{i}"]), float(row[f"dec{i}"]))
            for i in (1, 2, 3, 4))
        wcs = {k: (_f(row, k) if k not in ("ctype1", "ctype2") else row[k])
               for k in ("naxis1", "naxis2", "crpix1", "crpix2", "crval1",
                         "crval2", "ctype1", "ctype2", "cd1_1", "cd1_2",
                         "cd2_1", "cd2_2", "pxscal1", "pxscal2", "pa")}
        stem = f"{row['scan_id']}{int(row['frame_num']):03d}"
        base = (f"{IBE_BASE}/{row['scangrp']}/{row['scan_id']}/"
                f"{int(row['frame_num']):03d}/{stem}")
        products = {
            kind: {"url": f"{base}-w{band_num}{suffix}",
                   "md5_url": f"{base}-w{band_num}{suffix}.md5"}
            for kind, suffix in _PRODUCT_SUFFIX.items()}
        return Observation.build(
            archive_id=self.archive_id, collection=self.collection,
            release=row.get("wrelease", ""),
            native_key={"scan_id": row["scan_id"],
                        "frame_num": int(row["frame_num"]),
                        "band": band_num},
            band=f"W{band_num}",
            wavelength_um=BAND_WAVELENGTH_UM[band_num],
            t_start_mjd_utc=mjd - half, t_mid_mjd_utc=mjd,
            t_stop_mjd_utc=mjd + half, exptime_s=exptime,
            corners_icrs_deg=corners, wcs=wcs,
            quality_flags={k: _f(row, k) for k in
                           ("qual_frame", "qual_scan", "moon_sep",
                            "saa_sep", "magzp", "magzpunc")}
            | {"qa_status": row.get("qa_status")},
            products=products,
            extra={"image_set": row.get("image_set"),
                   "scangrp": row["scangrp"],
                   "date_obs": row.get("date_obs")},
            snapshot_id=snapshot_id, discovered_utc=discovered_utc,
        )

    def nominal_footprint(self, obs: Observation,
                          pad_arcsec: float = 0.0) -> WiseNominalFootprint:
        return WiseNominalFootprint(obs, pad_arcsec)

    def exact_footprint(self, obs: Observation,
                        products: ProductSet) -> WiseExactFootprint:
        msk = next(p for p in products.products if p.kind == "msk")
        return WiseExactFootprint(msk.path)

    def fetch(self, obs: Observation, kinds: Sequence[str], dest: Path,
              *, cutout: CutoutSpec | None = None) -> ProductSet:
        dest.mkdir(parents=True, exist_ok=True)
        products = []
        for kind in kinds:
            info = obs.products[kind]
            url = info["url"]
            params = {}
            if cutout is not None:
                params = {"center": f"{cutout.ra_deg},{cutout.dec_deg}",
                          "size": f"{cutout.size_pix}pix"}
            resp = self._session.get(url, params=params,
                                     timeout=self._timeout)
            resp.raise_for_status()
            name = url.rsplit("/", 1)[-1]
            if cutout is not None:
                name = f"cut{cutout.size_pix}-{name}"
            path = dest / name
            path.write_bytes(resp.content)
            if cutout is None:
                md5 = self._session.get(info["md5_url"],
                                        timeout=self._timeout)
                md5.raise_for_status()
                expected = md5.text.split()[0]
                actual = hashlib.md5(resp.content).hexdigest()
                if actual != expected:
                    raise IOError(f"md5 mismatch for {url}: "
                                  f"{actual} != {expected}")
                checksum = f"md5:{actual}"
            else:
                checksum = ("sha256:"
                            + hashlib.sha256(resp.content).hexdigest())
            products.append(LocalProduct(kind=kind, path=path,
                                         checksum=checksum))
        return ProductSet(observation_id=obs.observation_id,
                          products=tuple(products), cutout=cutout)
