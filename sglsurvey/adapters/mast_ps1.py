"""MAST Pan-STARRS1 DR2 single-epoch ("warp") archive adapter.

Discovery uses the `ps1filenames.py` listing service (one row per warp
covering the skycell that contains a query point), sampled over the
discovery cone; products come from the ps1images tree (full files) or
the `fitscut.cgi` cutout service. Facts verified before writing this
module are in surveys/panstarrs/notes/mast_recon.md.

Fourth adapter; first one that does not touch IRSA.
"""

from __future__ import annotations

import hashlib
import io
import json
import warnings
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np
import requests
from astropy.wcs import WCS

from sglsurvey.adapters.base import (ConeRegion, CutoutSpec, LocalProduct,
                                     MjdRange, ProductSet)
from sglsurvey.records import Observation
from sglsurvey.snapshots import SnapshotStore

LIST_URL = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
FITSCUT_URL = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
FILE_BASE = "https://ps1images.stsci.edu"
COLLECTION = "ps1/dr2/warp"

FILTER_WAVELENGTH_UM = {"g": 0.481, "r": 0.617, "i": 0.752,
                        "z": 0.866, "y": 0.962}
#: Nominal 3pi exposure times (s); the precise pass reads header EXPTIME.
NOMINAL_EXPTIME_S = {"g": 43.0, "r": 40.0, "i": 45.0, "z": 30.0, "y": 30.0}
PIX_ARCSEC = 0.25
NX, NY = 6240, 6243

# Mask bits (header MSKNAMnn / MSKVALnn, recon note).
MASK_BITS = {"DETECTOR": 1, "FLAT": 2, "DARK": 4, "BLANK": 8, "CTE": 16,
             "SAT": 32, "LOW": 64, "SUSPECT": 128, "CR": 256, "SPIKE": 512,
             "GHOST": 1024, "STREAK": 2048, "STARCORE": 4096,
             "CONV.BAD": 8192, "CONV.POOR": 16384, "MARK": 32768}
#: IPP's own bad template MASK.VALUE (8575) plus SPIKE, GHOST, STREAK and
#: STARCORE. SUSPECT and CONV.POOR are not fatal (hypotheses v1.0 section 9).
MASK_FATAL_TEMPLATE = 8575 | 512 | 1024 | 2048 | 4096  # = 16255

_PRODUCT_SUFFIX = {"img": ".fits", "msk": ".mask.fits", "wt": ".wt.fits"}


def _fits_open_any(path: Path):
    from astropy.io import fits
    from astropy.utils.exceptions import AstropyWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        hdul = fits.open(path)
        hdu = next(h for h in hdul if h.data is not None)
        return hdul, hdu


class Ps1SkycellFootprint:
    """Containment test from the cached skycell WCS (all warps of a
    skycell share it). Coarse stage only: the exposure covers only part
    of the skycell and GPC1 gaps are invisible here."""

    def __init__(self, wcs_dict: dict, pad_arcsec: float = 0.0):
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [wcs_dict["crpix1"], wcs_dict["crpix2"]]
        wcs.wcs.crval = [wcs_dict["crval1"], wcs_dict["crval2"]]
        wcs.wcs.cdelt = [wcs_dict["cdelt1"], wcs_dict["cdelt2"]]
        wcs.wcs.pc = [[wcs_dict["pc1_1"], wcs_dict["pc1_2"]],
                      [wcs_dict["pc2_1"], wcs_dict["pc2_2"]]]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        self._wcs = wcs
        self._nx, self._ny = wcs_dict["naxis1"], wcs_dict["naxis2"]
        self._pad_pix = pad_arcsec / PIX_ARCSEC

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        return self.contains_any(np.array([[ra_deg, dec_deg]]))

    def contains_any(self, radec_deg) -> bool:
        pix = self._wcs.wcs_world2pix(radec_deg, 0)
        with np.errstate(invalid="ignore"):
            x, y = pix[:, 0], pix[:, 1]
            p = self._pad_pix
            ok = ((x >= -0.5 - p) & (x <= self._nx - 0.5 + p)
                  & (y >= -0.5 - p) & (y <= self._ny - 0.5 + p))
        return bool(np.any(ok))


class Ps1ExactFootprint:
    """Exact usable-pixel test from a warp mask (full skycell file or a
    fitscut cutout). A pixel is usable when it is inside the array, the
    mask has no fatal bit, and — for full masks — the mask is finite.
    fitscut renders mask 0 as NaN, so NaN is treated as 0 here; pixels
    the exposure did not touch carry BLANK/CONV.BAD bits and are fatal.
    """

    FATAL_MASK = MASK_FATAL_TEMPLATE

    def __init__(self, msk_path: Path):
        hdul, hdu = _fits_open_any(msk_path)
        with hdul:
            data = np.asarray(hdu.data)
            self._wcs = WCS(hdu.header)
            self.header = dict(hdu.header)
        if data.dtype.kind == "f":
            data = np.where(np.isfinite(data), data, 0.0)
        self._mask = data.astype(np.int64)
        self._ny, self._nx = self._mask.shape

    def status(self, radec_deg):
        """Per-point (in_bounds, usable) boolean arrays for (N,2) deg."""
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
        _, usable = self.status(np.array([[ra_deg, dec_deg]]))
        return bool(usable[0])


class Ps1WarpAdapter:
    archive_id = "mast-ps1"
    collection = COLLECTION

    def __init__(self, session: requests.Session | None = None,
                 timeout_s: float = 600.0,
                 skycell_cache: Path | None = None):
        self._session = session or requests.Session()
        self._timeout = timeout_s
        self._cache_path = skycell_cache
        self._skycells: dict[str, dict] = {}
        if skycell_cache is not None and skycell_cache.exists():
            self._skycells = json.loads(skycell_cache.read_text())

    # -- discovery ---------------------------------------------------------
    @staticmethod
    def sample_points(region: ConeRegion, spacing_deg: float = 0.12
                      ) -> list[tuple[float, float]]:
        """Hexagonal sampling of the cone so that every skycell (~0.4 deg,
        overlapping) touching it is listed at least once."""
        pts = [(region.ra_deg, region.dec_deg)]
        r = region.radius_deg
        cosd = max(np.cos(np.deg2rad(region.dec_deg)), 1e-3)
        n = int(np.ceil(r / spacing_deg)) + 1
        for j in range(-n, n + 1):
            for i in range(-n, n + 1):
                dx = spacing_deg * (i + 0.5 * (j % 2))
                dy = spacing_deg * j * np.sqrt(3) / 2
                if (dx == 0 and dy == 0) or np.hypot(dx, dy) > r + 0.5 * spacing_deg:
                    continue
                pts.append(((region.ra_deg + dx / cosd) % 360.0,
                            float(np.clip(region.dec_deg + dy, -90, 90))))
        return pts

    def list_warps(self, ra_deg: float, dec_deg: float, store: SnapshotStore,
                   filters: str = "grizy") -> list[dict]:
        from datetime import datetime, timezone

        params = {"ra": f"{ra_deg:.6f}", "dec": f"{dec_deg:.6f}",
                  "filters": filters, "type": "warp"}
        request_utc = datetime.now(timezone.utc).isoformat()
        resp = self._session.get(LIST_URL, params=params,
                                 timeout=self._timeout)
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        rows = []
        if lines:
            header = lines[0].split()
            for line in lines[1:]:
                vals = line.split()
                if len(vals) == len(header):
                    rows.append(dict(zip(header, vals)))
        query = "&".join(f"{k}={v}" for k, v in params.items())
        snap = store.store(service_url=LIST_URL, query=query,
                           request_utc=request_utc,
                           response_bytes=resp.content, row_count=len(rows),
                           http_status=resp.status_code)
        for r in rows:
            r["_snapshot_id"] = snap.snapshot_id
            r["_request_utc"] = request_utc
        return rows

    def discover(self, region: ConeRegion, time_range: MjdRange,
                 store: SnapshotStore) -> Iterator[Observation]:
        seen = set()
        for ra, dec in self.sample_points(region):
            for row in self.list_warps(ra, dec, store):
                if row["filename"] in seen:
                    continue
                seen.add(row["filename"])
                mjd = float(row["mjd"])
                if not (time_range.start_mjd_utc <= mjd
                        <= time_range.stop_mjd_utc):
                    continue
                yield self._observation(row)

    def _observation(self, row: dict) -> Observation:
        filt = row["filter"]
        projcell, subcell = int(row["projcell"]), int(row["subcell"])
        stem = row["filename"]
        assert stem.endswith(".fits")
        stem = stem[:-5]
        exptime = NOMINAL_EXPTIME_S[filt]
        mjd_start = float(row["mjd"])
        half = exptime / 2.0 / 86400.0
        skycell = f"{projcell:04d}.{subcell:03d}"
        # The WCS is a property of the skycell, fetched lazily into the
        # adapter cache (skycell_wcs); it is deliberately NOT embedded in
        # the record so that observation identity is independent of it.
        wcs = {"skycell": skycell, "naxis1": NX, "naxis2": NY,
               "ctype1": "RA---TAN", "ctype2": "DEC--TAN",
               "cdelt_arcsec": PIX_ARCSEC}
        products = {kind: {"url": FILE_BASE + stem + suf}
                    for kind, suf in _PRODUCT_SUFFIX.items()}
        return Observation.build(
            archive_id=self.archive_id, collection=self.collection,
            release="dr2",
            native_key={"projcell": projcell, "subcell": subcell,
                        "filter": filt, "mjd_tag": stem.rsplit(".", 1)[-1]},
            band=filt, wavelength_um=FILTER_WAVELENGTH_UM[filt],
            t_start_mjd_utc=mjd_start, t_mid_mjd_utc=mjd_start + half,
            t_stop_mjd_utc=mjd_start + 2 * half, exptime_s=exptime,
            corners_icrs_deg=(), wcs=wcs,
            quality_flags={"badflag": int(row.get("badflag", 0))},
            products=products,
            extra={"skycell": skycell, "shortname": row.get("shortname"),
                   "list_ra": float(row["ra"]), "list_dec": float(row["dec"]),
                   "exptime_nominal": True},
            snapshot_id=row.get("_snapshot_id"),
            discovered_utc=row.get("_request_utc"),
        )

    # -- skycell WCS cache -------------------------------------------------
    def skycell_wcs(self, obs: Observation, dest: Path) -> dict:
        """WCS dict for the observation's skycell, fetched once per
        skycell from a full mask file and cached on disk."""
        key = obs.extra["skycell"]
        if key in self._skycells:
            return self._skycells[key]
        pset = self.fetch(obs, ["msk"], dest)
        hdul, hdu = _fits_open_any(pset.products[0].path)
        with hdul:
            h = hdu.header
            d = {"naxis1": int(h["NAXIS1"]), "naxis2": int(h["NAXIS2"]),
                 "crpix1": float(h["CRPIX1"]), "crpix2": float(h["CRPIX2"]),
                 "crval1": float(h["CRVAL1"]), "crval2": float(h["CRVAL2"]),
                 "cdelt1": float(h["CDELT1"]), "cdelt2": float(h["CDELT2"]),
                 "pc1_1": float(h.get("PC001001", h.get("PC1_1", -1.0))),
                 "pc1_2": float(h.get("PC001002", h.get("PC1_2", 0.0))),
                 "pc2_1": float(h.get("PC002001", h.get("PC2_1", 0.0))),
                 "pc2_2": float(h.get("PC002002", h.get("PC2_2", 1.0))),
                 "ctype1": "RA---TAN", "ctype2": "DEC--TAN",
                 "source_mask": pset.products[0].path.name}
        self._skycells[key] = d
        if self._cache_path is not None:
            self._cache_path.parent.mkdir(parents=True, exist_ok=True)
            self._cache_path.write_text(json.dumps(self._skycells, indent=1))
        return d

    # -- footprints --------------------------------------------------------
    def nominal_footprint(self, obs: Observation,
                          pad_arcsec: float = 0.0) -> Ps1SkycellFootprint:
        return Ps1SkycellFootprint(self._skycells[obs.extra["skycell"]],
                                   pad_arcsec)

    def exact_footprint(self, obs: Observation,
                        products: ProductSet) -> Ps1ExactFootprint:
        msk = next(p for p in products.products if p.kind == "msk")
        return Ps1ExactFootprint(msk.path)

    # -- fetch -------------------------------------------------------------
    def fetch(self, obs: Observation, kinds: Sequence[str], dest: Path,
              *, cutout: CutoutSpec | None = None) -> ProductSet:
        """Download full files or fitscut cutouts. No published per-file
        checksum applies to cutouts; recorded checksum = sha256 of the
        received bytes. A 404 raises FileNotFoundError."""
        dest.mkdir(parents=True, exist_ok=True)
        products = []
        for kind in kinds:
            url = obs.products[kind]["url"]
            name = url.rsplit("/", 1)[-1]
            if cutout is not None:
                name = f"cut{cutout.size_pix}-{name}"
            path = dest / name
            if not path.exists():
                if cutout is None:
                    resp = self._session.get(url, timeout=self._timeout)
                else:
                    rel = url[len(FILE_BASE):]
                    resp = self._session.get(
                        FITSCUT_URL,
                        params={"ra": f"{cutout.ra_deg:.6f}",
                                "dec": f"{cutout.dec_deg:.6f}",
                                "size": str(int(cutout.size_pix)),
                                "format": "fits", "red": rel},
                        timeout=self._timeout)
                if resp.status_code == 404:
                    raise FileNotFoundError(url)
                resp.raise_for_status()
                if not resp.content.startswith(b"SIMPLE"):
                    raise IOError(f"non-FITS response for {url}: "
                                  f"{resp.content[:80]!r}")
                path.write_bytes(resp.content)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            products.append(LocalProduct(kind=kind, path=path,
                                         checksum=f"sha256:{digest}"))
        return ProductSet(observation_id=obs.observation_id,
                          products=tuple(products), cutout=cutout)


# -- catalogs --------------------------------------------------------------
CATALOG_URL = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2"
DETECTION_COLS = ["objID", "detectID", "filterID", "obsTime", "ra", "dec",
                  "psfFlux", "psfFluxErr", "psfQfPerfect", "infoFlag",
                  "infoFlag2", "infoFlag3", "zp", "expTime", "imageID",
                  "psfMajorFWHM", "psfMinorFWHM"]
MEAN_COLS = ["objID", "raMean", "decMean", "nDetections", "qualityFlag",
             "gMeanPSFMag", "rMeanPSFMag", "iMeanPSFMag", "zMeanPSFMag",
             "yMeanPSFMag", "gMeanPSFMagErr", "rMeanPSFMagErr",
             "iMeanPSFMagErr", "zMeanPSFMagErr", "yMeanPSFMagErr",
             "gMeanPSFMagNpt", "rMeanPSFMagNpt", "iMeanPSFMagNpt",
             "zMeanPSFMagNpt", "yMeanPSFMagNpt"]
FILTER_ID = {1: "g", 2: "r", 3: "i", 4: "z", 5: "y"}


def catalog_cone(table: str, ra_deg: float, dec_deg: float, radius_deg: float,
                 columns: Sequence[str], store: SnapshotStore,
                 session: requests.Session | None = None,
                 pagesize: int = 500000, timeout_s: float = 900.0):
    """Paged cone query of a DR2 table; returns (rows, snapshot_ids)."""
    import csv
    from datetime import datetime, timezone

    session = session or requests.Session()
    rows, snaps = [], []
    page = 1
    while True:
        params = {"ra": f"{ra_deg:.6f}", "dec": f"{dec_deg:.6f}",
                  "radius": f"{radius_deg:.5f}", "pagesize": str(pagesize),
                  "page": str(page), "columns": "[" + ",".join(columns) + "]"}
        url = f"{CATALOG_URL}/{table}.csv"
        request_utc = datetime.now(timezone.utc).isoformat()
        resp = session.get(url, params=params, timeout=timeout_s)
        resp.raise_for_status()
        got = list(csv.DictReader(io.StringIO(resp.text)))
        snap = store.store(service_url=url,
                           query="&".join(f"{k}={v}" for k, v in params.items()),
                           request_utc=request_utc, response_bytes=resp.content,
                           row_count=len(got), http_status=resp.status_code)
        snaps.append(snap.snapshot_id)
        rows.extend(got)
        if len(got) < pagesize:
            break
        page += 1
    return rows, snaps
