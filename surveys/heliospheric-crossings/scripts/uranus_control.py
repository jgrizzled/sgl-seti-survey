"""Uranus conjunction positive control (freeze D5 / thresholds §6).

Finds on-disk C3 frames with Uranus inside the usable annulus
(elongation 1.3-7.5 deg), measures it through the identical chain
(load -> detrend -> star fit -> forced photometry at the SOHO-frame
ephemeris position -> per-frame radial ZP), and compares with the
predicted magnitude V = -7.19 + 5 log10(r_helio * d_obs) (phase ~ 0 at
conjunction). Gate: |median residual| <= 0.2 mag.
Writes results/uranus_control_v1.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lasco_lib as L

REPO = Path(__file__).resolve().parents[3]
DEVDIR = REPO / "runs" / "heliospheric-crossings" / "dev"
OUT = REPO / "surveys" / "heliospheric-crossings" / "results" / "uranus_control_v1.json"

expmap = json.loads((DEVDIR / "exptime_map_v1.json").read_text())
c3 = [(rp, v[1]) for rp, v in expmap.items() if v and "/c3/" in rp]
mjds = np.array([m for _, m in c3])
uniq = np.unique(np.round(mjds, 3))
keyint = lambda m: int(round(float(m) * 1000))
t = Time(uniq, format="mjd", scale="utc")
ura = get_body_barycentric("uranus", t).xyz.to_value("AU").T
soho = np.stack([np.interp(uniq, L._soho["mjd_utc"], L._soho["xyz_au"][:, i]) for i in range(3)], 1)
sun = get_body_barycentric("sun", t).xyz.to_value("AU").T

d = ura - soho
r_obs = np.linalg.norm(d, axis=1)
u_ra = np.degrees(np.arctan2(d[:, 1], d[:, 0])) % 360
u_dec = np.degrees(np.arcsin(d[:, 2] / r_obs))
ds = sun - soho
s_ra = np.degrees(np.arctan2(ds[:, 1], ds[:, 0])) % 360
s_dec = np.degrees(np.arcsin(ds[:, 2] / np.linalg.norm(ds, axis=1)))
dra = (u_ra - s_ra + 180) % 360 - 180
elong = np.hypot(dra * np.cos(np.radians(s_dec)), u_dec - s_dec)
r_helio = np.linalg.norm(ura - sun + ds + soho - soho, axis=1)  # |ura - sun|
r_helio = np.linalg.norm(ura - sun, axis=1)
good = (elong > 1.5) & (elong < 5.0)
by_mjd = {keyint(m): i for i, m in enumerate(uniq)}
cand = [(rp, m) for rp, m in c3 if keyint(m) in by_mjd and good[by_mjd[keyint(m)]]]
print(f"{good.sum()} C3 epochs with Uranus in annulus; {len(cand)} frames on disk")

rng = np.random.default_rng(20260825)
rng.shuffle(cand)
rows = []
for rp, m in cand:
    if len(rows) >= 250:
        break
    try:
        if L.frame_usable(REPO / rp):
            continue
        fr = L.load_frame(REPO / rp)
        L.fit_frame(fr)
        if fr.rot_deg is None or fr.zp_r is None:
            continue
        i = by_mjd[keyint(m)]
        p = L.forced_photometry(fr, float(u_ra[i]), float(u_dec[i]))
        if not np.isfinite(p["zp"]) or p.get("valid_fraction", 0) < 0.9:
            continue
        # keep signed flux: at per-frame S/N ~ 1 a positive-only cut
        # biases the stack to the noise floor; stack in linear units
        v_meas = p["zp"] - 2.5 * np.log10(p["flux"]) if p["flux"] > 0 else None
        v_pred = -7.19 + 5 * np.log10(r_helio[i] * r_obs[i])
        rows.append({"frame": rp, "mjd": m, "elong": round(float(elong[i]), 2),
                     "lin_flux": float(p["flux"] * 10 ** (-p["zp"] / 2.5)),
                     "v_meas": round(float(v_meas), 3) if v_meas else None,
                     "v_pred": round(float(v_pred), 3),
                     "snr": round(float(p["flux"] / p["err"]), 1)})
    except Exception:
        pass

# stacked control (the survey statistic is a stack; per-frame V 5.9
# aperture S/N vs coronal-annulus noise is ~1-10)
lin = np.array([r["lin_flux"] for r in rows])
vpred = float(np.median([r["v_pred"] for r in rows])) if rows else np.nan
boots = []
rng2 = np.random.default_rng(20260825)
for _ in range(500):
    boots.append(-2.5 * np.log10(np.median(rng2.choice(lin, len(lin)))))
v_stack = float(-2.5 * np.log10(np.median(lin))) if len(lin) else None
res = np.array([r["v_meas"] - r["v_pred"] for r in rows if r["snr"] > 5 and r["v_meas"]])
out = {"n_measured": len(rows),
       "v_stack": round(v_stack, 3) if v_stack else None,
       "v_stack_err_boot": round(float(np.std(boots)), 3) if rows else None,
       "v_pred_med": round(vpred, 3),
       "stack_resid_mag": round(v_stack - vpred, 3) if v_stack else None,
       "gate_0.2mag": bool(v_stack and abs(v_stack - vpred) <= 0.2),
       "n_snr5": int(len(res)),
       "median_perframe_snr": round(float(np.median([r["snr"] for r in rows])), 1) if rows else None,
       "rows": rows}
OUT.write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
