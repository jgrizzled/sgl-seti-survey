"""LASCO sunward-channels coverage intersect (freeze v1.0 stage 2).

For every sunward event x rung of `crossings/soho_v1`: the filled-cone
window (D2), the visibility-gated fraction (D7: occulter annuli from
`results/occulter_radii_v1.json`), and the per-day LASCO frame counts
from the archive daily-directory listings (SDAC primary, NRL for days
past the SDAC end; every listing snapshotted). Day-granularity frame
counts are the coverage measure — per-frame times come from headers at
the search stage.

Outputs: results/coverage_v1.json (per-event rows + aggregates),
listing snapshots under runs/heliospheric-crossings/coverage/listings/.
"""

from __future__ import annotations

import json
import re
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
EVENTS = REPO / "crossings" / "soho_v1" / "events.ecsv"
OCC = REPO / "surveys" / "heliospheric-crossings" / "results" / "occulter_radii_v1.json"
OUT = REPO / "surveys" / "heliospheric-crossings" / "results" / "coverage_v1.json"
SNAP = REPO / "runs" / "heliospheric-crossings" / "coverage" / "listings"

RSUN_AU = 0.00465047
AU_KM = 1.495978707e8
SDAC = "https://umbra.nascom.nasa.gov/pub/lasco_level05/"
NRL = "https://lasco-www.nrl.navy.mil/lz/level_05/"
SDAC_END_MJD = Time("2025-02-28").mjd  # measured at recon

RUNGS = {  # channel -> [(label, radius_au, camera)]
    "S1": [("1.2Rsun", 1.2 * RSUN_AU, "c2"), ("2.5Rsun", 2.5 * RSUN_AU, "c2"),
           ("0.1AU", 0.1, "c3")],
    "S2": [("0.1AU", 0.1, "c3")],
}
COMBOS = {"S1": ("outbound", "target"), "S2": ("inbound", "anti_target")}


def visible_days(b_au: float, r_au: float, v_kms: float,
                 ann_in_au: float, ann_out_au: float) -> tuple[float, float]:
    """(window_full_days, visible_days) for the filled-cone window of
    radius r and the camera visibility annulus [ann_in, ann_out]."""
    def t_at(radius):  # days from t_ca to b_e = radius (one side)
        if radius <= b_au:
            return 0.0
        return float(np.sqrt(radius**2 - b_au**2) * AU_KM / (v_kms * 86400.0))

    half = t_at(r_au)
    lo, hi = t_at(min(ann_in_au, r_au)), t_at(min(ann_out_au, r_au))
    return 2 * half, 2 * max(0.0, hi - lo)


class Listings:
    def __init__(self):
        SNAP.mkdir(parents=True, exist_ok=True)
        self.counts: dict[tuple[str, str], dict] = {}

    def fetch(self, cam: str, yymmdd: str, mjd: float) -> dict:
        key = (cam, yymmdd)
        if key in self.counts:
            return self.counts[key]
        snap = SNAP / f"{cam}_{yymmdd}.html"
        rec = {"n_frames": 0, "source": None}
        if snap.exists():
            html = snap.read_text()
            rec = {"n_frames": len(re.findall(r'href="\d+\.fts"', html)),
                   "source": html.splitlines()[0][4:] if html.startswith("<!--") else "snapshot"}
        else:
            order = [("sdac", SDAC), ("nrl", NRL)]
            if mjd > SDAC_END_MJD:
                order = [("nrl", NRL), ("sdac", SDAC)]
            for src, base in order:
                url = f"{base}{yymmdd}/{cam}/"
                for attempt in range(3):
                    try:
                        req = urllib.request.Request(url)
                        html = urllib.request.urlopen(req, timeout=30).read().decode(
                            "utf-8", "replace")
                        n = len(re.findall(r'href="\d+\.fts"', html))
                        if n > 0:
                            snap.write_text(f"<!--{src}-->\n" + html)
                            rec = {"n_frames": n, "source": src}
                        break
                    except urllib.error.HTTPError as e:
                        if e.code == 404:
                            break
                        time.sleep(1 + attempt)
                    except Exception:
                        time.sleep(1 + attempt)
                if rec["n_frames"]:
                    break
            if not rec["n_frames"]:
                snap.write_text("<!--empty-->\n")
                rec["source"] = "empty"
        self.counts[key] = rec
        return rec


def main() -> None:
    t = Table.read(EVENTS)
    occ = json.loads(OCC.read_text())["adopted"]
    ann = {"c2": [x * RSUN_AU for x in occ["c2_rsun"]],
           "c3": [x * RSUN_AU for x in occ["c3_rsun"]]}

    rows = []
    dayset: set[tuple[str, str, float]] = set()
    for ch, (direction, side) in COMBOS.items():
        m = (np.asarray(t["link_direction"]) == direction) & (np.asarray(t["side"]) == side)
        sub = t[m]
        for label, r_au, cam in RUNGS[ch]:
            sel = sub[np.asarray(sub["b_min_au"]) <= r_au]
            for ev in sel:
                b, v = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
                full_d, vis_d = visible_days(b, r_au, v, *ann[cam])
                mjd_ca = Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd
                d0, d1 = int(np.floor(mjd_ca - full_d / 2)), int(np.floor(mjd_ca + full_d / 2))
                days = [Time(d, format="mjd").strftime("%y%m%d") for d in range(d0, d1 + 1)]
                rows.append({
                    "event_id": ev["event_id"], "target_id": ev["target_id"],
                    "channel": ch, "rung": label, "camera": cam,
                    "b_min_rsun": round(b / RSUN_AU, 3),
                    "t_ca_utc": ev["t_ca_utc"][:16],
                    "window_days": round(full_d, 2),
                    "visible_days": round(vis_d, 2),
                    "visible_fraction": round(vis_d / full_d, 3) if full_d else 0.0,
                    "days": days, "mjd_ca": mjd_ca,
                })
                for i, d in enumerate(days):
                    dayset.add((cam, d, mjd_ca))

    print(f"{len(rows)} event-rung rows; {len({(c, d) for c, d, _ in dayset})} unique (cam,day) listings")
    lst = Listings()
    uniq = sorted({(c, d) for c, d, _ in dayset})
    mjd_by = {(c, d): m for c, d, m in dayset}
    with ThreadPoolExecutor(max_workers=6) as ex:
        list(ex.map(lambda cd: lst.fetch(cd[0], cd[1], mjd_by[cd]), uniq))

    for row in rows:
        counts = [lst.counts[(row["camera"], d)]["n_frames"] for d in row["days"]]
        row["n_frames_window_days"] = int(sum(counts))
        row["n_days_with_frames"] = int(sum(1 for c in counts if c))
        row["covered"] = bool(row["n_frames_window_days"] > 0)
        row["searchable_geom"] = bool(row["visible_fraction"] > 0)
        del row["days"], row["mjd_ca"]

    # aggregates per channel x rung x target
    agg: dict[str, dict] = {}
    for row in rows:
        k = f'{row["channel"]}_{row["rung"]}'
        a = agg.setdefault(k, {}).setdefault(row["target_id"], {
            "events": 0, "covered": 0, "covered_and_visible": 0,
            "frames": [], "vis_frac": []})
        a["events"] += 1
        a["covered"] += row["covered"]
        a["covered_and_visible"] += row["covered"] and row["searchable_geom"]
        a["frames"].append(row["n_frames_window_days"])
        a["vis_frac"].append(row["visible_fraction"])
    for k, per_t in agg.items():
        for tid, a in per_t.items():
            a["median_frames_per_window"] = float(np.median(a.pop("frames")))
            a["visible_fraction"] = float(np.median(a.pop("vis_frac")))

    sources = {}
    for rec in lst.counts.values():
        sources[rec["source"]] = sources.get(rec["source"], 0) + 1
    out = {
        "stage": "coverage_v1", "events_input": str(EVENTS.relative_to(REPO)),
        "occulter_annuli_rsun": occ,
        "listing_sources": sources,
        "n_event_rung_rows": len(rows),
        "aggregates": agg,
        "rows": rows,
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"listing_sources": sources}, indent=1))
    for k in sorted(agg):
        tot = {s: sum(a[s] for a in agg[k].values()) for s in ("events", "covered", "covered_and_visible")}
        print(k, tot)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
