"""Pre-freeze flux-scale check (draft hypotheses §4): measure NSC DR2
stars on fetched instcal CCDs with plain aperture photometry and
determine (a) the MAGZERO convention empirically (per-second vs
per-exposure counts), (b) the per-frame offset and scatter of
header MAGZERO vs a star-derived zero point.

This is a validation, not the survey flux scale: the frozen pipeline
will star-calibrate per frame through the identical matched filter
(PS1 lesson 1). Target here: convention resolved, star ZP scatter
quantified, gross errors (> 0.2 mag systematic per frame after
calibration) excluded.

Usage: uv run python surveys/decam/scripts/flux_scale_check.py
Writes surveys/decam/results/flux_scale_check.json (+ snapshots under
runs/decam/fluxcheck_v1/).
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))

REPO = Path(__file__).resolve().parents[3]
CUT_DIR = REPO / "runs" / "decam" / "products" / "cut"
RUN_DIR = REPO / "runs" / "decam" / "fluxcheck_v1"
OUT = REPO / "surveys" / "decam" / "results" / "flux_scale_check.json"

APERTURE_PIX = 5.0
ANNULUS_PIX = (10.0, 15.0)
STAR_MAG = (16.0, 19.5)
MIN_EXPTIME_S = 30.0
MAX_FRAMES_PER_BAND = 4
NSC_BAND_COL = {"g": "gmag", "r": "rmag", "i": "imag", "z": "zmag",
                "Y": "ymag", "u": "umag"}


def aperture_phot(img, x, y, r_ap, r_in, r_out):
    """Background-annulus aperture sum at (x, y); returns counts or
    None near edges / non-finite backgrounds."""
    h, w = img.shape
    if not (r_out < x < w - r_out and r_out < y < h - r_out):
        return None
    x0, x1 = int(x - r_out - 1), int(x + r_out + 2)
    y0, y1 = int(y - r_out - 1), int(y + r_out + 2)
    sub = img[y0:y1, x0:x1]
    yy, xx = np.mgrid[y0:y0 + sub.shape[0], x0:x0 + sub.shape[1]]
    rr = np.hypot(xx - x, yy - y)
    ann = sub[(rr >= r_in) & (rr < r_out) & np.isfinite(sub)]
    if len(ann) < 30:
        return None
    bkg = np.median(ann)
    ap = (rr < r_ap) & np.isfinite(sub)
    if ap.sum() < np.pi * r_ap**2 * 0.9:
        return None
    return float(np.sum(sub[ap] - bkg))


def main() -> None:
    from astropy.io import fits
    from astropy.wcs import WCS
    from dl import queryClient as qc

    store = SnapshotStore(RUN_DIR)
    manifest = [json.loads(l) for l in (CUT_DIR / "manifest.jsonl").open()
                if '"schema": "v1"' in l]
    by_band = defaultdict(list)
    for m in manifest:
        hdr = m.get("header") or {}
        band = hdr.get("band")
        if (not m["files"] or m["missing"] or band not in NSC_BAND_COL
                or (hdr.get("exptime") or 0) < MIN_EXPTIME_S):
            continue
        by_band[band].append(m)
    picks = []
    for band, ms in sorted(by_band.items()):
        ms = sorted(ms, key=lambda m: m["observation_id"])
        step = max(len(ms) // MAX_FRAMES_PER_BAND, 1)
        picks.extend(ms[::step][:MAX_FRAMES_PER_BAND])
    print(f"{len(picks)} frames selected "
          f"({', '.join(f'{b}:{len(v)}' for b, v in sorted(by_band.items()))} usable)")

    frames = []
    per_band = defaultdict(list)
    for m in picks:
        band = m["header"]["band"]
        img_key = next(k for k in m["files"] if k.startswith("image:"))
        path = CUT_DIR / m["files"][img_key]
        with fits.open(path) as hdul:
            img = np.asarray(hdul[1].data, dtype=float)
            wcs = WCS(hdul[1].header)
            magzero = float(hdul[0].header.get("MAGZERO", np.nan))
            exptime = float(hdul[0].header.get("EXPTIME", np.nan))
            ny, nx = img.shape
        corners = wcs.all_pix2world(
            np.array([[1, 1], [nx, 1], [nx, ny], [1, ny]]), 1)
        cra, cdec = corners[:, 0].mean(), corners[:, 1].mean()
        rad = 0.17  # covers a 9' x 18' CCD from centre
        col = NSC_BAND_COL[band]
        sql = (f"SELECT ra,dec,{col},class_star,ndet FROM nsc_dr2.object "
               f"WHERE q3c_radial_query(ra,dec,{cra:.6f},{cdec:.6f},{rad}) "
               f"AND {col} BETWEEN {STAR_MAG[0]} AND {STAR_MAG[1]} "
               f"AND class_star > 0.8 AND ndet >= 5")
        request_utc = datetime.now(timezone.utc).isoformat()
        text = qc.query(sql=sql, fmt="csv", timeout=120)
        store.store(service_url="datalab:queryClient/query", query=sql,
                    request_utc=request_utc, response_bytes=text.encode(),
                    row_count=text.count("\n") - 1, http_status=200)
        lines = text.strip().splitlines()[1:]
        offs = []
        for line in lines:
            ra, dec, mag = (float(v) for v in line.split(",")[:3])
            if mag > 99:
                continue
            px, py = wcs.all_world2pix([[ra, dec]], 0)[0]
            counts = aperture_phot(img, px, py, APERTURE_PIX,
                                   *ANNULUS_PIX)
            if counts is None or counts <= 0:
                continue
            offs.append(mag + 2.5 * np.log10(counts))
        if len(offs) < 10:
            print(f"  {path.name}: only {len(offs)} stars, skipped")
            continue
        offs = np.array(offs)
        zp_star = float(np.median(offs))
        mad = float(1.4826 * np.median(np.abs(offs - zp_star)))
        frames.append({
            "file": path.name, "band": band, "exptime_s": exptime,
            "magzero_header": magzero, "n_stars": len(offs),
            "zp_star_counts": round(zp_star, 3),
            "zp_scatter_mad": round(mad, 3),
            "zp_star_minus_magzero": round(zp_star - magzero, 3),
            "zp_star_minus_magzero_pers": round(
                zp_star - magzero - 2.5 * np.log10(exptime), 3)})
        per_band[band].append(frames[-1])
        print(f"  {path.name} ({band}, {exptime:.0f}s): {len(offs)} stars, "
              f"ZP* {zp_star:.3f} (MAD {mad:.3f}); "
              f"ZP*-MAGZERO {zp_star - magzero:+.3f}, "
              f"-2.5logEXP {zp_star - magzero - 2.5 * np.log10(exptime):+.3f}",
              flush=True)

    band_summary = {}
    for band, fs in sorted(per_band.items()):
        d_raw = [f["zp_star_minus_magzero"] for f in fs]
        d_per = [f["zp_star_minus_magzero_pers"] for f in fs]
        band_summary[band] = {
            "n_frames": len(fs),
            "median_offset_counts_conv": round(float(np.median(d_raw)), 3),
            "spread_counts_conv": round(float(np.ptp(d_raw)), 3),
            "median_offset_per_second_conv": round(
                float(np.median(d_per)), 3),
            "spread_per_second_conv": round(float(np.ptp(d_per)), 3)}
        print(f"[{band}] counts-conv offset {band_summary[band]['median_offset_counts_conv']:+.3f} "
              f"(spread {band_summary[band]['spread_counts_conv']:.3f}); "
              f"per-second-conv {band_summary[band]['median_offset_per_second_conv']:+.3f} "
              f"(spread {band_summary[band]['spread_per_second_conv']:.3f})")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "method": ("aperture photometry (r=5 pix, annulus 10-15) of "
                   "NSC DR2 stars 16-19.5 mag, class_star>0.8, ndet>=5; "
                   "ZP* = median(mag_cat + 2.5 log10 counts)"),
        "frames": frames, "band_summary": band_summary}, indent=2))
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
