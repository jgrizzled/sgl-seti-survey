"""Exceedance adjudication helpers (freeze §6 veto ladder).

`census(mjd, ra, dec)`: separations of the major planets and the
brightest asteroids (H < 7 list, Horizons OBSERVER ephemeris with
CENTER='@-144') from a patch at one epoch. `event_detail(role, unit,
event)`: the arc series of the source and control patches for
persistence / split-half / anticorrelation tests, from the reduce
series dump, plus the Gaia-template excess (rule 7: static content).
`unit_report(role, unit)`: every included event's z, z_template,
split halves, pulse and the per-orbit pattern.

    python adjudicate.py <role> <unit_key> [event_id]
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import solo_geometry as G
import series as S

BRIGHT_ASTEROIDS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23,
                    24, 25, 27, 28, 29, 30, 31, 32, 37, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49,
                    51, 52, 63, 64, 65, 68, 69, 88, 89, 115, 192, 216, 230, 324, 349, 354, 387,
                    471, 511, 532, 704]


def horizons_radec(body: str, mjd: float):
    t0 = Time(mjd - 0.02, format="mjd").iso[:16]
    t1 = Time(mjd + 0.03, format="mjd").iso[:16]
    q = {"format": "text", "COMMAND": f"'{body}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
         "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'@-144'", "START_TIME": f"'{t0}'",
         "STOP_TIME": f"'{t1}'", "STEP_SIZE": "'1h'", "QUANTITIES": "'1,9'",
         "CSV_FORMAT": "'YES'", "ANG_FORMAT": "'DEG'"}
    try:
        txt = urllib.request.urlopen("https://ssd.jpl.nasa.gov/api/horizons.api?"
                                     + urllib.parse.urlencode(q), timeout=120).read().decode()
        line = txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()[0]
        p = [x.strip() for x in line.split(",")]
        return float(p[3]), float(p[4]), float(p[5]) if p[5] not in ("n.a.", "") else np.nan
    except Exception:
        return None


def census(mjd: float, ra: float, dec: float, radius_deg: float = 0.5, asteroids: bool = True) -> dict:
    s = G.radec_to_vec(ra, dec)
    out = {"planets": {}, "asteroids": {}}
    for body, dirs in G.planet_directions(mjd, bodies=("mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune")).items():
        out["planets"][body] = round(float(np.degrees(np.arccos(np.clip(np.dot(dirs[0], s), -1, 1)))), 3)
    if asteroids:
        for n in BRIGHT_ASTEROIDS:
            r = horizons_radec(f"{n};", mjd)
            if r is None:
                continue
            sep = float(np.degrees(np.arccos(np.clip(np.dot(G.radec_to_vec(r[0], r[1]), s), -1, 1))))
            if sep < radius_deg:
                out["asteroids"][n] = {"sep_deg": round(sep, 3), "V": r[2]}
    return out


def event_detail(role: str, ukey: str, eid: str, version: str = "v1") -> dict:
    ser = json.loads((S.SER / f"series_{role}_{version}.json").read_text())
    idx = json.loads((S.SER / f"index_{role}_v1.json").read_text())
    key = f"{ukey}|{eid}"
    se = ser[key]
    mj = np.array(se["mjd"])
    d0 = np.array(se["D"]["0"])
    sd = S._robust_sd(d0)
    n = len(d0)
    h = n // 2
    z1 = np.mean(d0[:h]) / (sd / np.sqrt(h))
    z2 = np.mean(d0[h:]) / (sd / np.sqrt(n - h))
    imax = int(np.argmax(np.minimum(d0[:-1], d0[1:]))) if n > 1 else 0
    neigh = d0[max(0, imax - 2): imax + 3]
    ctrl = {k: np.array(v) for k, v in se["D"].items() if k != "0"}
    dt = np.array(se.get("D_template_src", []))
    return {"unit": ukey, "event": eid, "n_arc": n, "z": round(float(np.mean(d0) / (sd / np.sqrt(n))), 2),
            "z_template": round(float(np.mean(dt) / (S._robust_sd(dt) / np.sqrt(len(dt)) + 1e-30)), 2) if len(dt) >= 10 else None,
            "split_half_z": [round(float(z1), 2), round(float(z2), 2)],
            "slope_D_vs_time": round(float(np.polyfit(mj - mj.mean(), d0, 1)[0] if n > 2 else 0), 4),
            "pulse": {"mjd": float(mj[imax]), "utc": Time(mj[imax], format="mjd").iso[:16],
                      "D": round(float(d0[imax]), 3), "sigma": round(float((d0[imax] - np.median(d0)) / sd), 2),
                      "neighbours_D": [round(float(x), 3) for x in neigh],
                      "controls_same_frame_median_D": round(float(np.median([c[imax] for c in ctrl.values()])), 3) if ctrl else None},
            "mean_D": round(float(np.mean(d0)), 4), "sd_D": round(sd, 4),
            "patch_radec": idx["patches"][key]["radec"][0], "orbit": idx["events"][key]["orbit"],
            "tile": idx["patches"][key]["tile"], "bet0": idx["patches"][key]["bet0"], "ladder": idx["patches"][key]["ladder"],
            "mjd_ca": idx["events"][key]["mjd_ca"], "n_gaia_src": idx["patches"][key]["n_gaia"][0],
            "template_src_stars": idx["patches"][key]["stars"][0][:6]}


def unit_report(role: str, ukey: str, version: str = "v1") -> dict:
    res = json.loads((S.SURV / "results" / f"{role}_search_{version}.json").read_text())
    r = res["units"][ukey]
    evs = {}
    for eid, ev in r["events"].items():
        if ev["status"] != "included":
            evs[eid] = {"status": ev["status"]}
            continue
        evs[eid] = event_detail(role, ukey, eid, version)
    return {"unit": ukey, "stats": {s: r[s] for s in ("S_stack", "S_event", "S_pulse")},
            "per_patch": r["per_patch"], "events": evs}


if __name__ == "__main__":
    role, ukey = sys.argv[1:3]
    if len(sys.argv) > 3:
        d = event_detail(role, ukey, sys.argv[3])
        print(json.dumps(d, indent=1))
        ra, dec = d["patch_radec"]
        print(json.dumps(census(d["pulse"]["mjd"], ra, dec), indent=1))
    else:
        print(json.dumps(unit_report(role, ukey), indent=1))
