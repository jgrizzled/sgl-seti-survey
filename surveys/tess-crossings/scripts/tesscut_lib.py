"""Minimal TESScut access layer for the tess-crossings survey.

Fetches astrocut FFI cutout cubes with local caching and a sha256
manifest (runs/tess-crossings/products/manifest.jsonl), and loads
them into a small Cube object. Facts verified in
notes/tesscut_recon.md. The full sglsurvey adapter (Observation
records etc.) is deferred to the search stage; coverage refinement
and the freezes only need cadence bookkeeping.

Time handling: PIXELS.TIME is BTJD (TDB, barycentre-corrected via
TIMECORR). Window intersection happens in UTC at the observer, so
cadences convert as  jd_tdb_sc = TIME - TIMECORR + 2457000  ->
Time(..., scale="tdb").utc.mjd.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

REPO = Path(__file__).resolve().parents[3]
RUN = REPO / "runs" / "tess-crossings"
PRODUCTS = RUN / "products"
ASTROCUT = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut"


@dataclass
class Cube:
    label: str
    sector: int
    camera: int
    ccd: int
    mjd_utc: np.ndarray        # per cadence, spacecraft UTC
    btjd: np.ndarray
    flux: np.ndarray           # (n, ny, nx) e-/s
    flux_err: np.ndarray
    flux_bkg: np.ndarray
    quality: np.ndarray
    wcs: WCS
    path: Path


def fetch_cube(label: str, ra: float, dec: float, sector: int,
               size_px: int, session: requests.Session | None = None
               ) -> Path:
    """Download (or reuse) one astrocut zip; returns the FITS path."""
    PRODUCTS.mkdir(parents=True, exist_ok=True)
    dest = PRODUCTS / f"{label}-s{sector:04d}-{size_px}px"
    dest.mkdir(exist_ok=True)
    existing = sorted(dest.glob("*.fits"))
    if existing:
        return existing[0]
    session = session or requests.Session()
    r = session.get(ASTROCUT, params={
        "ra": f"{ra:.6f}", "dec": f"{dec:.6f}",
        "y": str(size_px), "x": str(size_px), "units": "px",
        "sector": str(sector)}, timeout=1200)
    r.raise_for_status()
    z = zipfile.ZipFile(io.BytesIO(r.content))
    names = z.namelist()
    if not names:
        raise FileNotFoundError(f"no cutout returned for {label} "
                                f"s{sector}")
    z.extractall(dest)
    path = dest / names[0]
    rec = {"label": label, "sector": sector, "ra": ra, "dec": dec,
           "size_px": size_px, "file": str(path.relative_to(REPO)),
           "sha256": "sha256:" + hashlib.sha256(
               path.read_bytes()).hexdigest(),
           "fetched_utc": Time.now().isot}
    with (PRODUCTS / "manifest.jsonl").open("a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return path


def load_cube(path: Path, label: str = "") -> Cube:
    with fits.open(path) as f:
        d = f[1].data
        hdr0 = f[0].header
        w = WCS(f[2].header)
        tc = np.nan_to_num(np.asarray(d["TIMECORR"], float))
        jd_tdb_sc = np.asarray(d["TIME"], float) - tc + 2457000.0
        ok = np.isfinite(jd_tdb_sc)
        mjd = np.full(len(jd_tdb_sc), np.nan)
        if ok.any():
            mjd[ok] = Time(jd_tdb_sc[ok], format="jd",
                           scale="tdb").utc.mjd
        return Cube(
            label=label or path.stem,
            sector=int(hdr0.get("SECTOR", 0)),
            camera=int(hdr0.get("CAMERA", 0)),
            ccd=int(hdr0.get("CCD", 0)),
            mjd_utc=mjd, btjd=np.asarray(d["TIME"], float),
            flux=np.asarray(d["FLUX"], np.float32),
            flux_err=np.asarray(d["FLUX_ERR"], np.float32),
            flux_bkg=np.asarray(d["FLUX_BKG"], np.float32),
            quality=np.asarray(d["QUALITY"], np.int64),
            wcs=w, path=path)


def tic_cone(ra: float, dec: float, radius_deg: float,
             session: requests.Session | None = None):
    """TIC cone via the MAST invoke API; returns row dicts."""
    session = session or requests.Session()
    req = {"service": "Mast.Catalogs.Filtered.Tic.Position.Rows",
           "format": "json",
           "params": {"columns": "ID,ra,dec,Tmag,GAIAmag,pmRA,pmDEC",
                      "ra": ra, "dec": dec, "radius": radius_deg,
                      "filters": []}}
    r = session.post("https://mast.stsci.edu/api/v0/invoke",
                     data={"request": json.dumps(req)}, timeout=300)
    r.raise_for_status()
    doc = r.json()
    RUN.mkdir(parents=True, exist_ok=True)
    with (RUN / "tic_queries.jsonl").open("a") as fh:
        fh.write(json.dumps({"ra": ra, "dec": dec,
                             "radius": radius_deg,
                             "utc": Time.now().isot,
                             "n": len(doc.get("data", []))}) + "\n")
    return doc.get("data", [])
