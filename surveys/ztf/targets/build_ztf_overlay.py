"""ZTF overlay of the shared registry: grade + order corridors for the
ZTF scale-up without changing membership (targets/README.md policy).

Per corridor (WISE corridor keys): antipode centre at MJD 60000,
ZTF-visibility (corridor Dec > -28), public good-quality frame count
in a 0.2 deg box around the antipode (one snapshotted IBE query),
dominant field/quadrant and the antipode's pixel position in it
(primary-grid CCD-gap flag, the pilot's key lesson), secondary-field
frame count, filter counts, and epoch span. Writes overlay_v1.json/.md.

Usage: uv run python surveys/ztf/targets/build_ztf_overlay.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, load_target_registry

from sglsurvey.adapters.base import ConeRegion, MjdRange
from sglsurvey.adapters.irsa_ztf import ZtfNominalFootprint, ZtfSciAdapter
from sglsurvey.geometry import GeometryContext, locus_radec
from sglsurvey.snapshots import SnapshotStore

REPO = Path(__file__).resolve().parents[3]
from sglsurvey.corridors import CORRIDOR_OF, MEMBERS  # noqa: E402

OUT_JSON = REPO / "surveys" / "ztf" / "targets" / "overlay_v1.json"
OUT_MD = REPO / "surveys" / "ztf" / "targets" / "overlay_v1.md"
SNAP_DIR = REPO / "runs" / "ztf" / "overlay_v1"
DEC_LIMIT = -28.0
NX, NY = ZtfNominalFootprint.NX, ZtfNominalFootprint.NY
PILOT = {"ross128", "epsind", "proxima"}


def main() -> None:
    registry = load_target_registry(REPO / "registries" / "pilot_wise_2026.yaml")
    ctx = GeometryContext.ztf_default()
    adapter = ZtfSciAdapter()
    store = SnapshotStore(SNAP_DIR)
    t = Time(60000.0, format="mjd")
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
               "ztf_visible": dec > DEC_LIMIT, "pilot": corridor in PILOT}
        if dec > DEC_LIMIT:
            obs = list(adapter.discover(ConeRegion(ra, dec, 0.1),
                                        MjdRange(58150.0, 61300.0), store))
            good = [o for o in obs if not o.quality_flags.get("bad_quality")]
            quads = Counter((o.native_key["field"], o.native_key["ccdid"],
                             o.native_key["qid"]) for o in good)
            fields = Counter(o.native_key["field"] for o in good)
            prim = fields.most_common(1)[0][0] if fields else None
            # antipode pixel position in the primary field's quadrants
            pos = {}
            for (f, c, q), n in quads.items():
                if f != prim:
                    continue
                o = next(o for o in good if (o.native_key["field"],
                                             o.native_key["ccdid"],
                                             o.native_key["qid"]) == (f, c, q))
                x, y = ZtfNominalFootprint(o)._wcs.wcs_world2pix([[ra, dec]], 0)[0]
                pos[f"c{c:02d}q{q}"] = {"n": n, "x": round(float(x)), "y": round(float(y)),
                                        "inside": bool(0 <= x < NX and 0 <= y < NY)}
            inside_any = any(v["inside"] for v in pos.values())
            # distance (pix) from the antipode to the nearest live quadrant edge
            def edge_margin(v):
                return min(v["x"], NX - 1 - v["x"], v["y"], NY - 1 - v["y"])
            margin = max((edge_margin(v) for v in pos.values()), default=None)
            mjds = np.array([o.t_mid_mjd_utc for o in good])
            row.update({
                "frames_public_good": len(good),
                "frames_by_filter": dict(Counter(o.band for o in good)),
                "primary_field": prim,
                "primary_field_frames": fields.get(prim, 0),
                "secondary_fields": {str(k): v for k, v in fields.items() if k != prim},
                "antipode_in_primary_quadrant": inside_any,
                "antipode_margin_pix": margin,
                "grid_flag": ("gap" if not inside_any else
                              "edge" if (margin is not None and margin < 400) else "ok"),
                "primary_quadrants": pos,
                "epoch_span_mjd": ([round(float(mjds.min()), 1),
                                    round(float(mjds.max()), 1)] if len(mjds) else None),
            })
            print(f"{corridor:14s} dec={dec:+6.1f} frames={len(good):5d} "
                  f"primary={prim} {row['grid_flag']:4s} margin={margin}", flush=True)
        else:
            print(f"{corridor:14s} dec={dec:+6.1f} NOT VISIBLE", flush=True)
        rows.append(row)

    # ordering: ok before edge before gap; then by frame count
    order = {"ok": 0, "edge": 1, "gap": 2}
    vis = [r for r in rows if r["ztf_visible"]]
    vis.sort(key=lambda r: (order[r["grid_flag"]], -r["frames_public_good"]))
    queue = [r["corridor"] for r in vis if not r["pilot"]]
    OUT_JSON.write_text(json.dumps({"dec_limit": DEC_LIMIT, "rows": rows,
                                    "queue": queue}, indent=1))
    lines = ["# ZTF overlay v1 (2026-08-20)", "",
             f"{len(vis)} of {len(rows)} corridors have Dec > {DEC_LIMIT}. "
             "Grading only; membership is the universal list's. `grid_flag`: "
             "antipode inside a primary-grid quadrant (ok), within 400 pix of "
             "its edge (edge), or in a CCD gap (gap — only secondary fields cover it).",
             "", "| corridor | endpoints | Dec | frames (g/r/i) | primary field | flag | margin px | secondary frames |",
             "|---|---|---|---|---|---|---|---|"]
    for r in vis:
        fb = r["frames_by_filter"]
        lines.append(f"| {r['corridor']}{' (pilot)' if r['pilot'] else ''} | "
                     f"{', '.join(r['endpoints'])} | {r['antipode_dec_deg']:+.1f} | "
                     f"{r['frames_public_good']} ({fb.get('zg',0)}/{fb.get('zr',0)}/{fb.get('zi',0)}) | "
                     f"{r['primary_field']} | {r['grid_flag']} | {r['antipode_margin_pix']} | "
                     f"{sum(r['secondary_fields'].values())} |")
    lines += ["", "Not visible: " + ", ".join(r["corridor"] for r in rows if not r["ztf_visible"])]
    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"\nqueue ({len(queue)}): {queue}")


if __name__ == "__main__":
    main()
