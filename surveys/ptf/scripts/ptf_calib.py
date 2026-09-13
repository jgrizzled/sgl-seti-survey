"""Shared PTF flux-scale helpers (hypotheses §4): PS1 DR2 mean-object
calibrators through the frozen transforms (g = gMeanPSFMag;
R = r - 0.153 (r - i) - 0.117, Jordi et al. 2006), the dwarf-locus
colour restriction and quality gates fixed at the PTF crossings dev
stage (surveys/ptf-crossings/scripts/dev_search.py), and the per-frame
zero point through the identical matched filter.

Calibrator source for the corridor survey: the PS1 corridor screen
snapshots (runs/panstarrs/screen_v1, one 0.2-deg `mean` cone per
corridor — the same 62 corridors, PS1 and PTF share the Dec > -30 sky).
"""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import numpy as np

from sglsurvey.records import read_records

REPO = Path(__file__).resolve().parents[3]
PS1_SCREEN = REPO / "runs" / "panstarrs" / "screen_v1"
STAR_MAG = {"R": (15.5, 19.5), "g": (16.0, 20.0)}
MIN_STARS = 5
MAX_SCATTER = 0.2
_cache: dict = {}


def mould_r(r, i):
    return r - 0.153 * (r - i) - 0.117


def calibrators_from_rows(rows, band):
    """[(ra, dec, predicted PTF mag)] under the frozen calibrator rules."""
    lo, hi = STAR_MAG[band]
    out = []
    for r in rows:
        try:
            nd = int(float(r["nDetections"]))
            rm, im = float(r["rMeanPSFMag"]), float(r["iMeanPSFMag"])
            gm = float(r["gMeanPSFMag"])
            if band == "R":
                re_, ie = (float(r["rMeanPSFMagErr"]),
                           float(r["iMeanPSFMagErr"]))
                if min(rm, im) < -100 or max(re_, ie) > 0.1:
                    continue
                if not (0.0 <= rm - im <= 0.8):
                    continue
                m = mould_r(rm, im)
            else:
                if (min(gm, rm) < -100
                        or float(r["gMeanPSFMagErr"]) > 0.1):
                    continue
                if not (0.2 <= gm - rm <= 1.2):
                    continue
                m = gm
        except (KeyError, TypeError, ValueError):
            continue
        if nd >= 5 and lo <= m <= hi:
            out.append((float(r["raMean"]), float(r["decMean"]), m))
    return out


def ps1_mean_rows(corridor):
    """PS1 DR2 mean objects of the corridor's screen cone."""
    if ("rows", corridor) not in _cache:
        stats = json.loads((PS1_SCREEN / "catalog_stats.json").read_text())
        snaps = _cache.setdefault("snaps", {
            s["snapshot_id"]: s for s in read_records(
                PS1_SCREEN / "records" / "query_snapshot.jsonl")})
        rows = []
        for sid in stats[corridor]["mean_snapshots"]:
            p = PS1_SCREEN / snaps[sid]["response_path"]
            rows.extend(csv.DictReader(io.StringIO(p.read_text())))
        _cache[("rows", corridor)] = rows
    return _cache[("rows", corridor)]


def corridor_calibrators(corridor, band):
    key = ("cal", corridor, band)
    if key not in _cache:
        _cache[key] = np.array(calibrators_from_rows(
            ps1_mean_rows(corridor), band)).reshape(-1, 3)
    return _cache[key]


def star_zeropoint(fm, stars, min_stars=MIN_STARS,
                   max_scatter=MAX_SCATTER):
    """Per-frame ZP = median(m_pred + 2.5 log10 F_mf) over calibrators
    with good_frac >= 0.7 and S/N >= 5; a 3-sigma-clipped robust scatter
    <= max_scatter and >= min_stars survivors gate the frame. Returns
    (zp or None, n_used, scatter, n_candidates_in_frame)."""
    if len(stars) == 0:
        return None, 0, None, 0
    f, v, g = fm.sample(stars[:, 0], stars[:, 1])
    inb = np.isfinite(f)
    ok = (inb & (f > 0) & (g >= 0.7) & np.isfinite(v)
          & (f / np.sqrt(np.maximum(v, 1e-30)) >= 5))
    if ok.sum() < min_stars:
        return None, int(ok.sum()), None, int(inb.sum())
    zps = stars[ok, 2] + 2.5 * np.log10(f[ok])
    med = float(np.median(zps))
    mad = float(1.4826 * np.median(np.abs(zps - med)))
    keep = np.abs(zps - med) < 3 * max(mad, 0.02)
    if keep.sum() < min_stars:
        return None, int(keep.sum()), mad, int(inb.sum())
    zp = float(np.median(zps[keep]))
    mad = float(1.4826 * np.median(np.abs(zps[keep] - zp)))
    if mad > max_scatter:
        return None, int(keep.sum()), mad, int(inb.sum())
    return zp, int(keep.sum()), mad, int(inb.sum())
