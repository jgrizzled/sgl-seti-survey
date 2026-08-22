"""PS1 overlay of the shared registry: grade + order corridors for the
Pan-STARRS1 scale-up without changing membership (targets/README.md
policy). Mirrors surveys/ztf/targets/build_ztf_overlay.py.

Per corridor (WISE corridor keys): antipode centre at MJD 56000 (3pi
midpoint), PS1 visibility (corridor Dec > -30), warp count per filter
from one snapshotted `ps1filenames` listing at the antipode, epoch
span, the parallax-phase split of those epochs (the pilot's key
lesson: 3pi revisits a field at the same season), and the calibrator
density (DR2 mean-table stars with 15 < r < 20.5 within 0.05 deg) that
decides whether the per-warp star zero point will be available.
Writes overlay_v1.json/.md.

Usage: uv run python surveys/panstarrs/targets/build_ps1_overlay.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, load_target_registry

from sglsurvey.adapters.mast_ps1 import Ps1WarpAdapter, catalog_cone
from sglsurvey.geometry import GeometryContext, locus_radec
from sglsurvey.snapshots import SnapshotStore

REPO = Path(__file__).resolve().parents[3]
from sglsurvey.corridors import MEMBERS  # noqa: E402

OUT_JSON = REPO / "surveys" / "panstarrs" / "targets" / "overlay_v1.json"
OUT_MD = REPO / "surveys" / "panstarrs" / "targets" / "overlay_v1.md"
SNAP_DIR = REPO / "runs" / "panstarrs" / "overlay_v1"
DEC_LIMIT = -30.0
PILOT = {"ross128", "epsind", "proxima"}
CAL_RADIUS_DEG = 0.05
CAL_MAG = (15.0, 20.5)


def phase_split(mjds: np.ndarray) -> tuple[int, int]:
    """Epochs within +-91.3 d of the circular-median day-of-year vs the
    rest (same definition as sample_tensor.py)."""
    ang = (mjds % 365.25) / 365.25 * 2 * np.pi
    ref = (np.arctan2(np.sin(ang).mean(), np.cos(ang).mean())
           % (2 * np.pi)) / (2 * np.pi) * 365.25
    d = (mjds % 365.25) - ref
    d = (d + 182.625) % 365.25 - 182.625
    p0 = int((np.abs(d) < 91.3).sum())
    return p0, int(len(mjds) - p0)


def main() -> None:
    registry = load_target_registry(REPO / "registries" / "pilot_wise_2026.yaml")
    ctx = GeometryContext.ps1_default()
    adapter = Ps1WarpAdapter()
    store = SnapshotStore(SNAP_DIR)
    t = Time(56000.0, format="mjd")
    rows = []
    for corridor, members in MEMBERS.items():
        members = [m for m in members if m in registry.ids]
        if not members:
            continue
        pts, _ = locus_radec(ctx, registry[members[0]], Role.RX, t,
                             tolerance_arcsec=30.0)
        ra, dec = float(np.mean(pts[:, 0])), float(np.mean(pts[:, 1]))
        row = {"corridor": corridor, "endpoints": members,
               "antipode_ra_deg": round(ra, 4), "antipode_dec_deg": round(dec, 4),
               "ps1_visible": dec > DEC_LIMIT, "pilot": corridor in PILOT}
        if dec > DEC_LIMIT:
            try:
                warps = adapter.list_warps(ra, dec, store)
            except Exception as exc:
                print(f"{corridor:14s} listing FAILED {exc}", flush=True)
                warps = []
            mjds = np.array([float(w["mjd"]) for w in warps])
            p0, p1 = phase_split(mjds) if len(mjds) else (0, 0)
            try:
                stars, _ = catalog_cone(
                    "mean", ra, dec, CAL_RADIUS_DEG,
                    ["objID", "rMeanPSFMag", "rMeanPSFMagNpt"], store)
                n_cal = sum(1 for s in stars
                            if CAL_MAG[0] < float(s["rMeanPSFMag"]) < CAL_MAG[1]
                            and int(float(s["rMeanPSFMagNpt"])) >= 3)
            except Exception as exc:
                print(f"{corridor:14s} catalog FAILED {exc}", flush=True)
                n_cal = None
            row.update({
                "warps": len(warps),
                "warps_by_filter": dict(Counter(w["filter"] for w in warps)),
                "skycell": (f"{int(warps[0]['projcell']):04d}."
                            f"{int(warps[0]['subcell']):03d}") if warps else None,
                "epoch_span_mjd": ([round(float(mjds.min()), 1),
                                    round(float(mjds.max()), 1)] if len(mjds) else None),
                "phase_split": [p0, p1],
                "minor_phase_frac": round(min(p0, p1) / max(len(mjds), 1), 3),
                "calibrators_r_per_0p05deg": n_cal,
                "calibrator_flag": (None if n_cal is None else
                                    "sparse" if n_cal < 40 else "ok"),
            })
            print(f"{corridor:14s} dec={dec:+6.1f} warps={len(warps):4d} "
                  f"phase={p0}:{p1} cal={n_cal}", flush=True)
        else:
            print(f"{corridor:14s} dec={dec:+6.1f} NOT VISIBLE", flush=True)
        rows.append(row)

    vis = [r for r in rows if r["ps1_visible"]]
    vis.sort(key=lambda r: -r["warps"])
    queue = [r["corridor"] for r in vis if not r["pilot"] and r["warps"] > 0]
    OUT_JSON.write_text(json.dumps({"dec_limit": DEC_LIMIT, "rows": rows,
                                    "queue": queue}, indent=1))
    lines = ["# PS1 overlay v1 (2026-08-20)", "",
             f"{len(vis)} of {len(rows)} corridors have Dec > {DEC_LIMIT}. "
             "Grading only; membership is the universal list's. `phase` = "
             "epochs within +-91 d of the median season vs the rest (the "
             "parallax-phase veto needs both); `cal` = DR2 stars 15 < r < 20.5 "
             "within 0.05 deg (sparse < 40: per-warp star ZP will often fall "
             "back to the filter median).",
             "", "| corridor | endpoints | Dec | warps (g/r/i/z/y) | skycell | phase | cal |",
             "|---|---|---|---|---|---|---|"]
    for r in vis:
        fb = r["warps_by_filter"]
        lines.append(f"| {r['corridor']}{' (pilot)' if r['pilot'] else ''} | "
                     f"{', '.join(r['endpoints'])} | {r['antipode_dec_deg']:+.1f} | "
                     f"{r['warps']} ({'/'.join(str(fb.get(b, 0)) for b in 'grizy')}) | "
                     f"{r['skycell']} | {r['phase_split'][0]}:{r['phase_split'][1]} | "
                     f"{r['calibrators_r_per_0p05deg']} ({r['calibrator_flag']}) |")
    lines += ["", "Not visible: " + ", ".join(r["corridor"] for r in rows if not r["ps1_visible"])]
    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"\nqueue ({len(queue)}): {queue}")


if __name__ == "__main__":
    main()
