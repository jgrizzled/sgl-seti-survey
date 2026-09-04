"""Which registry targets did K2 itself observe, and at what impact parameter?

Addendum to footprint_intersect.py: the sky-only test found six targets
whose *on-star* position lies on K2 silicon in some campaign (channel A,
any b). This confirms against MAST (K2 timeseries products within 1'
of the target's campaign-epoch position) and records b(t) at the
campaign mid-time from the kepler_v1 geometry: b = r_helio x sin(sun
elongation) for an observer on the axis-free side, so a quadrature-
observed star sits at b ~ 0.6-1.0 AU (the wide rung, plan §3 'A 1.0 AU').
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from footprint_intersect import (REPO, RUN, EVENTS, SC_TABLE, angsep, mast,  # noqa: E402
                                 radec_to_vec, load_fovs, campaign_ranges, PLAN_CLIP_D)
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

OUT = HERE.parent / "results" / "onstar_k2_v1.json"
IN = HERE.parent / "results" / "footprint_intersect_v1.json"


def main():
    store = SnapshotStore(RUN)
    d = json.load(open(IN))
    sc = np.load(SC_TABLE)
    fovs = load_fovs()
    ranges = campaign_ranges(store)
    t = Table.read(EVENTS)
    pairs = sorted({(e["target_id"], k) for e in d["on_silicon_any_time"]
                    if e["channel"] == "A" for k in e["on_silicon_any_time"]})
    out = {}
    for tid, camp in pairs:
        ts = ranges[camp]["timeseries"]
        p0, p1 = (Time(x, scale="utc").mjd for x in fovs[camp]["planned"])
        lo, hi = max(ts["t_min_mjd"], p0 - PLAN_CLIP_D), min(ts["t_max_mjd"], p1 + PLAN_CLIP_D)
        mid = 0.5 * (lo + hi)
        # star position at the campaign epoch: nearest kepler_v1 inbound event's star_icrs
        rows = t[(t["target_id"] == tid) & (t["link_direction"] == "inbound")]
        ev_mjd = Time(list(rows["t_ca_utc"]), format="isot", scale="utc").mjd
        j = int(np.argmin(np.abs(ev_mjd - mid)))
        ra, dec = float(rows["star_icrs_ra_deg"][j]), float(rows["star_icrs_dec_deg"][j])
        # geometry at campaign start/mid/end
        geo = {}
        for name, m in (("start", lo), ("mid", mid), ("end", hi)):
            xyz = np.array([np.interp(m, sc["mjd_utc"], sc["xyz_au"][:, k]) for k in range(3)])
            sun = get_body_barycentric("sun", Time(m, format="mjd", scale="utc")).xyz.to_value("AU")
            o = xyz - sun
            r = np.linalg.norm(o)
            s = radec_to_vec(ra, dec)
            cos_e = -np.dot(o, s) / r
            b = r * np.sqrt(max(1 - cos_e ** 2, 0.0))
            geo[name] = {"mjd": round(float(m), 2), "sun_elongation_deg": round(float(np.degrees(np.arccos(np.clip(cos_e, -1, 1)))), 1),
                         "b_au": round(float(b), 3)}
        dra = 1 / 60 / max(np.cos(np.radians(dec)), 0.2)
        rows_m = mast(store, "SELECT obs_collection, sequence_number, target_name, dataproduct_type, COUNT(*) AS n, "
                             "MIN(t_min) AS tmin, MAX(t_max) AS tmax FROM dbo.obspointing WHERE obs_collection='K2' "
                             f"AND s_ra BETWEEN {ra - dra:.5f} AND {ra + dra:.5f} AND s_dec BETWEEN {dec - 1/60:.5f} AND {dec + 1/60:.5f} "
                             "GROUP BY obs_collection, sequence_number, target_name, dataproduct_type")
        out[f"{tid}:{camp}"] = {"ra": ra, "dec": dec, "campaign": camp, "geometry": geo,
                                "k2_products_within_1arcmin": [{"sequence": int(r[1]), "target_name": r[2], "type": r[3],
                                                                "n": int(r[4]), "t_min": float(r[5]), "t_max": float(r[6])} for r in rows_m]}
        print(tid, camp, geo["mid"], [(r[1], r[2], r[3], r[4]) for r in rows_m])
    OUT.write_text(json.dumps(out, indent=1))
    print("wrote", OUT.relative_to(REPO))


if __name__ == "__main__":
    main()
