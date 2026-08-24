"""Per-frame star-calibrated zero points through the IDENTICAL matched
filter (PS1 lesson 1; hypotheses v1.0 §4): for every complete
stage-2 exposure, build the v2 flux map and measure NSC DR2 stars
(band mag 16.0-21.0, class_star > 0.7, ndet >= 5, from the screen_v1
object snapshots); zp_star = median(mag_cat + 2.5 log10 F_mf).

Writes runs/decam/zeropoints.jsonl (one row per observation id) — the
profile's build_map uses zp_star and never the header MAGZERO
(measured unreliable, results/flux_scale_check.json).

Usage: uv run python surveys/decam/scripts/calibrate_zeropoints.py
       [--corridors c1 c2 ...] [--limit N]
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time as _time
from pathlib import Path

import numpy as np

from sglsurvey.adapters.noirlab_decam import DQ_FATAL_DEFAULT
from sglsurvey.photometry import build_flux_map_decam
from sglsurvey.records import read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from decam_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
CUT = REPO / "runs" / "decam" / "products" / "cut"
DQ = REPO / "runs" / "decam" / "products" / "dqmask"
SCREEN = REPO / "runs" / "decam" / "screen_v1"
OUT = REPO / "runs" / "decam" / "zeropoints.jsonl"
BAND_COL = {"g": "gmag", "r": "rmag", "i": "imag", "z": "zmag",
            "Y": "ymag", "u": "umag", "VR": "vrmag"}
STAR_MAG = (16.0, 21.0)
MIN_STARS = 6


def object_snapshots():
    """corridor-cone object CSVs: [(ra, dec, rows_dict_list)]"""
    cones = []
    for s in read_records(SCREEN / "records" / "query_snapshot.jsonl"):
        q = s["query"]
        if "nsc_dr2.object" not in q or "COUNT" in q.upper():
            continue
        m = re.search(r"q3c_radial_query\(ra,dec,([-\d.]+),([-\d.]+),", q)
        if m:
            rows = list(csv.DictReader(io.StringIO(
                (SCREEN / s["response_path"]).read_text())))
            cones.append((float(m.group(1)), float(m.group(2)), rows))
    return cones


def stars_near(cones, ra, dec, band):
    col = BAND_COL.get(band)
    best, bd = None, 1e9
    for cra, cdec, rows in cones:
        d = np.hypot((cra - ra) * np.cos(np.deg2rad(dec)), cdec - dec)
        if d < bd:
            best, bd = rows, d
    if best is None or col is None or bd > 0.3:
        return np.empty((0, 3))
    out = []
    for r in best:
        try:
            mag = float(r[col])
            if not (STAR_MAG[0] <= mag <= STAR_MAG[1]):
                continue
            if float(r.get("class_star") or 0) <= 0.7:
                continue
            if int(float(r.get("ndet") or 0)) < 5:
                continue
            out.append((float(r["ra"]), float(r["dec"]), mag))
        except (KeyError, TypeError, ValueError):
            continue
    return np.array(out) if out else np.empty((0, 3))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corridors", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    only = set(args.corridors) if args.corridors else None

    best = {}
    for line in open(CUT / "manifest.jsonl"):
        m = json.loads(line)
        if m.get("schema") == "v1":
            best[m["observation_id"]] = m
    rows = [m for m in best.values() if m["files"] and not m["missing"]]
    if only:
        rows = [m for m in rows if set(m["corridor"]) & only]
    done = set()
    if OUT.exists():
        done = {json.loads(l)["observation_id"] for l in open(OUT)}
    todo = [m for m in rows if m["observation_id"] not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(rows)} exposures, {len(todo)} to calibrate")

    cones = object_snapshots()
    t0 = _time.monotonic()
    with OUT.open("a") as fh:
        for i, m in enumerate(todo):
            band = m["header"]["band"]
            ccd_files = [
                (CUT / m["files"][f"image:{c}"],
                 CUT / m["files"][f"wtmap:{c}"]
                 if f"wtmap:{c}" in m["files"] else None, c)
                for c in m["ccds"] if f"image:{c}" in m["files"]]
            rec = {"observation_id": m["observation_id"], "band": band,
                   "zp_star": None, "n_stars": 0, "zp_mad": None,
                   "magzero_header": m["header"].get("magzero")}
            fm = None
            try:
                fm = build_flux_map_decam(
                    ccd_files, DQ / m["dqmask"], DQ_FATAL_DEFAULT, band,
                    m["header"].get("mjd_obs") or 0.0,
                    (m["center_ra_deg"], m["center_dec_deg"]),
                    m["size_pix"])
            except Exception as exc:
                rec["error"] = str(exc)[:120]
            if fm is not None:
                st = stars_near(cones, m["center_ra_deg"],
                                m["center_dec_deg"], band)
                if len(st):
                    f, v, g = fm.sample(st[:, 0], st[:, 1])
                    ok = (np.isfinite(f) & (f > 0) & (g > 0.9)
                          & np.isfinite(v)
                          & (f / np.sqrt(np.maximum(v, 1e-30)) > 7))
                    if ok.sum() >= MIN_STARS:
                        zps = st[ok, 2] + 2.5 * np.log10(f[ok])
                        med = float(np.median(zps))
                        mad = float(1.4826 * np.median(np.abs(zps - med)))
                        keep = np.abs(zps - med) < 3 * max(mad, 0.02)
                        rec.update(
                            zp_star=round(float(np.median(zps[keep])), 4),
                            n_stars=int(keep.sum()),
                            zp_mad=round(mad, 4),
                            fwhm_arcsec=round(float(fm.fwhm_arcsec), 3),
                            var_source=fm.var_source)
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if (i + 1) % 20 == 0:
                print(f"  {i + 1}/{len(todo)} "
                      f"({(i + 1) / (_time.monotonic() - t0):.2f}/s)",
                      flush=True)
    n_ok = sum(1 for l in open(OUT)
               if json.loads(l).get("zp_star") is not None)
    print(f"done: {n_ok} calibrated of {len(rows)}")


if __name__ == "__main__":
    main()
