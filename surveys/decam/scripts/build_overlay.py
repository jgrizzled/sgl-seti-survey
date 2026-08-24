"""Build the DECam southern overlay (surveys/decam/targets/overlay_v1):
per-corridor archive coverage grades for the 15 southern corridors, on
the pattern of the ZTF/PS1 overlays — grades and ordering only, never
membership (plan targets rules).

Per corridor, from one snapshotted adv_search count query over the
discovery-cone box: exposure count, era span, distinct calendar months
(parallax-phase proxy), per-band counts after the draft quality cut
(EXPTIME >= 30 s, grizY), and a crowding note from the NSC object
density. Writes overlay_v1.json + overlay_v1.md.

Usage: uv run python surveys/decam/scripts/build_overlay.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from sglseti import Role, load_target_registry

from sglsurvey.geometry import GeometryContext, discovery_cone
from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from decam_corridors import MEMBERS, PILOT_CORRIDORS  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
RUN_DIR = REPO / "runs" / "decam" / "overlay_v1"
OUT_DIR = REPO / "surveys" / "decam" / "targets"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
TIME_RANGE_MJD = (56000.0, 62000.0)
MIN_EXPTIME_S = 30.0
BANDS = ("g", "r", "i", "z", "Y")
#: NSC objects per sq deg above which a corridor is graded crowded
#: (Galactic plane / LMC); pilot cones run 30k-500k per sq deg.
CROWDED_PER_SQDEG = 3.0e5


def main() -> None:
    import math
    import urllib.request

    from dl import queryClient as qc

    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.decam_default()
    store = SnapshotStore(RUN_DIR)

    rows = []
    for corridor, members in sorted(MEMBERS.items()):
        cone = discovery_cone(ctx, registry[members[0]],
                              (Role.RX, Role.TX), *TIME_RANGE_MJD)
        r_deg = cone.radius_deg + 1.1
        cosd = max(np.cos(np.deg2rad(cone.dec_deg)), 1e-3)
        q = {"outfields": ["caldat", "ifilter", "exposure", "obs_type",
                           "EXPNUM", "prod_type"],
             "search": [["instrument", "decam"], ["proc_type", "instcal"],
                        ["prod_type", "image"],
                        ["ra_center", (cone.ra_deg - r_deg / cosd) % 360,
                         (cone.ra_deg + r_deg / cosd) % 360],
                        ["dec_center", cone.dec_deg - r_deg,
                         cone.dec_deg + r_deg]]}
        body = json.dumps(q).encode()
        url = "https://astroarchive.noirlab.edu/api/adv_search/find/?limit=100000"
        request_utc = datetime.now(timezone.utc).isoformat()
        req = urllib.request.Request(
            url, data=body, headers={"Content-Type": "application/json"})
        raw = urllib.request.urlopen(req, timeout=120).read()
        exps = json.loads(raw)[1:]
        store.store(service_url=url, query=body.decode(),
                    request_utc=request_utc, response_bytes=raw,
                    row_count=len(exps), http_status=200)

        good = [e for e in exps
                if (e.get("exposure") or 0) >= MIN_EXPTIME_S
                and (e.get("ifilter") or "?").split()[0] in BANDS
                and e.get("obs_type") in (None, "object")]
        dates = sorted(e["caldat"] for e in good if e.get("caldat"))
        months = sorted({d[5:7] for d in dates})
        band_n = Counter((e.get("ifilter") or "?").split()[0]
                         for e in good)

        sql = (f"SELECT COUNT(*) FROM nsc_dr2.object WHERE "
               f"q3c_radial_query(ra,dec,{cone.ra_deg:.5f},"
               f"{cone.dec_deg:.5f},0.2)")
        request_utc = datetime.now(timezone.utc).isoformat()
        text = qc.query(sql=sql, fmt="csv", timeout=120)
        store.store(service_url="datalab:queryClient/query", query=sql,
                    request_utc=request_utc, response_bytes=text.encode(),
                    row_count=1, http_status=200)
        n_obj = int(text.strip().splitlines()[1])
        density = n_obj / (math.pi * 0.2**2)

        # grade: ok / sparse / crowded (both possible; crowded wins for
        # ordering, sparse for feasibility)
        n_good = len(good)
        grade = ("crowded" if density > CROWDED_PER_SQDEG
                 else "sparse" if n_good < 40 else "ok")
        rows.append({
            "corridor": corridor, "members": members,
            "cone": {"ra_deg": round(cone.ra_deg, 4),
                     "dec_deg": round(cone.dec_deg, 4),
                     "radius_deg": round(cone.radius_deg, 4)},
            "n_exposures_raw": len(exps), "n_exposures_quality": n_good,
            "bands": {b: band_n.get(b, 0) for b in BANDS},
            "date_range": [dates[0], dates[-1]] if dates else None,
            "n_calendar_months": len(months),
            "nsc_objects_per_sqdeg": int(density),
            "grade": grade,
            "pilot": corridor in PILOT_CORRIDORS})
        print(f"{corridor:12s} {grade:8s} {n_good:6d} quality exposures "
              f"({len(exps)} raw), months {len(months):2d}, "
              f"density {density / 1e3:.0f}k/deg2, bands "
              + " ".join(f"{b}:{band_n.get(b, 0)}" for b in BANDS),
              flush=True)

    # queue: ok grades first by quality-exposure count desc, then
    # crowded, then sparse; pilot corridors excluded (already running)
    order = {"ok": 0, "crowded": 1, "sparse": 2}
    queue = [r["corridor"] for r in sorted(
        rows, key=lambda r: (order[r["grade"]],
                             -r["n_exposures_quality"]))
        if not r["pilot"]]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    doc = {"version": "overlay_v1",
           "built_utc": datetime.now(timezone.utc).isoformat(),
           "quality_cut": {"min_exptime_s": MIN_EXPTIME_S,
                           "bands": list(BANDS),
                           "obs_type": "object-or-null"},
           "rows": rows, "queue": queue}
    (OUT_DIR / "overlay_v1.json").write_text(json.dumps(doc, indent=2))

    lines = ["# DECam southern overlay v1", "",
             f"Built {doc['built_utc'][:10]}; quality cut EXPTIME >= "
             f"{MIN_EXPTIME_S:.0f} s, grizY, obs_type object. Grades "
             "order the queue and never change membership.", "",
             "| corridor | grade | quality exp | months | density/deg2 |"
             " g/r/i/z/Y | span |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in sorted(rows, key=lambda r: r["corridor"]):
        b = r["bands"]
        lines.append(
            f"| {r['corridor']}{' (pilot)' if r['pilot'] else ''} "
            f"| {r['grade']} | {r['n_exposures_quality']} "
            f"| {r['n_calendar_months']} "
            f"| {r['nsc_objects_per_sqdeg']:,} "
            f"| {b['g']}/{b['r']}/{b['i']}/{b['z']}/{b['Y']} "
            f"| {r['date_range'][0][:4]}-{r['date_range'][-1][:4]} |"
            if r["date_range"] else
            f"| {r['corridor']} | {r['grade']} | 0 | 0 | - | - | - |")
    lines += ["", "Queue (non-pilot): " + ", ".join(queue), ""]
    (OUT_DIR / "overlay_v1.md").write_text("\n".join(lines))
    print(f"\nwrote {OUT_DIR.relative_to(REPO)}/overlay_v1.{{json,md}}")


if __name__ == "__main__":
    main()
