"""Exceedance adjudication helpers (freeze §6 veto ladder).

`census(mjd, ra, dec)`: separations of the major planets, Earth and
the Moon (astropy ephemeris from the PSP position) and of the
brightest asteroids (H < 7 list, JPL Horizons OBSERVER ephemeris with
CENTER='@-96'; SkyBoT has no PSP code) from a patch at one
epoch. `event_detail(role, unit, event)`: the arc series of the source
and control patches for persistence / split-half / anticorrelation
tests, from the reduce series dump.

    python adjudicate.py <role> <unit_key> <event_id>
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
import psp_geometry as G
import series as S

# numbered asteroids with H < ~7 (the ones that can reach V < 11 near conjunction)
BRIGHT_ASTEROIDS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23,
                    24, 25, 27, 28, 29, 30, 31, 32, 37, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49,
                    51, 52, 63, 64, 65, 68, 69, 88, 89, 115, 192, 216, 230, 324, 349, 354, 387,
                    471, 511, 532, 704]


def horizons_radec(body: str, mjd: float) -> tuple[float, float, float] | None:
    t0 = Time(mjd - 0.02, format="mjd").iso[:16]
    t1 = Time(mjd + 0.03, format="mjd").iso[:16]
    q = {"format": "text", "COMMAND": f"'{body}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
         "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'@-96'", "START_TIME": f"'{t0}'",
         "STOP_TIME": f"'{t1}'", "STEP_SIZE": "'1h'", "QUANTITIES": "'1,9'",
         "CSV_FORMAT": "'YES'", "ANG_FORMAT": "'DEG'"}
    try:
        txt = urllib.request.urlopen("https://ssd.jpl.nasa.gov/api/horizons.api?"
                                     + urllib.parse.urlencode(q), timeout=120).read().decode()
        line = txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()[0]
        p = [x.strip() for x in line.split(",")]
        return float(p[3]), float(p[4]), float(p[5])
    except Exception:
        return None


def census(mjd: float, ra: float, dec: float, radius_deg: float = 0.5,
           asteroids: bool = True) -> dict:
    s = G.radec_to_vec(ra, dec)
    out = {"planets": {}, "asteroids": {}}
    for body, dirs in G.planet_directions(mjd, bodies=("mercury", "venus", "earth", "mars",
                                                        "jupiter", "saturn", "uranus", "neptune")).items():
        sep = float(np.degrees(np.arccos(np.clip(np.dot(dirs[0], s), -1, 1))))
        out["planets"][body] = round(sep, 3)
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
    template = None
    sd = S._robust_sd(d0)
    n = len(d0)
    h = n // 2
    z1 = np.mean(d0[:h]) / (sd / np.sqrt(h))
    z2 = np.mean(d0[h:]) / (sd / np.sqrt(n - h))
    imax = int(np.argmax(d0))
    neigh = d0[max(0, imax - 2): imax + 3]
    ctrl = {k: np.array(v) for k, v in se["D"].items() if k != "0"}
    same_frame_ctrl = np.array([np.median([c[imax] for c in ctrl.values()])])
    return {"unit": ukey, "event": eid, "n_arc": n, "z": round(float(np.mean(d0) / (sd / np.sqrt(n))), 2),
            "split_half_z": [round(float(z1), 2), round(float(z2), 2)],
            "pulse": {"mjd": float(mj[imax]), "utc": Time(mj[imax], format="mjd").iso[:16],
                      "D": round(float(d0[imax]), 3), "sigma": round(float((d0[imax] - np.median(d0)) / sd), 2),
                      "neighbours_D": [round(float(x), 3) for x in neigh],
                      "controls_same_frame_median_D": round(float(same_frame_ctrl[0]), 3)},
            "mean_D": round(float(np.mean(d0)), 4), "sd_D": round(sd, 4),
            "patch_radec": idx["patches"][key]["radec"][0], "encounter": idx["events"][key]["encounter"],
            "mjd_ca": idx["events"][key]["mjd_ca"]}


if __name__ == "__main__":
    role, ukey, eid = sys.argv[1:4]
    d = event_detail(role, ukey, eid)
    print(json.dumps(d, indent=1))
    ra, dec = d["patch_radec"]
    print(json.dumps(census(d["pulse"]["mjd"], ra, dec), indent=1))
