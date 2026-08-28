"""ZP era validation (thresholds §6): (a) per-year stability of the
level-0.5 star-ZP chain from the colour-ensemble records (the per-frame
star ZP absorbs CCD degradation by construction — this checks for
anomalies/discontinuities, esp. across the level-1 era end 2017-08);
(b) matched-exposure spot check: the same 2010-06-15 C3 exposure
measured through the level-0.5 chain and the level-1 chain (both
star-calibrated to Hipparcos V) — median per-star magnitude difference,
0.1 mag gate. Writes results/zp_validation_v1.json."""

from __future__ import annotations

import gzip
import io
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lasco_lib as L

REPO = Path(__file__).resolve().parents[3]
DEVDIR = REPO / "runs" / "heliospheric-crossings" / "dev"
RECON = REPO / "runs" / "heliospheric-crossings" / "recon"
OUT = REPO / "surveys" / "heliospheric-crossings" / "results" / "zp_validation_v1.json"

# (a) per-year ensemble stability (residuals about the global radial+colour model)
rec = np.load(DEVDIR / "color_ensemble_records.npz")
c3 = rec["cam"] == "c3"
zp, r, bv, mjd = rec["zp"][c3], rec["r"][c3], rec["bv"][c3], rec["mjd"][c3]
# global radial profile + colour term (same construction as the ensemble fit)
prof = {}
for lo in np.arange(0, r.max() + 3, 3):
    m = (r >= lo) & (r < lo + 3)
    if m.sum() >= 10:
        prof[lo + 1.5] = np.median(zp[m])
keys = sorted(prof)
resid = zp - np.interp(r, keys, [prof[k] for k in keys]) - 0.429 * (bv - 0.65)
year = 1858.88 + mjd / 365.25
per_year = {}
for y in range(1996, 2027):
    m = (year >= y) & (year < y + 1)
    if m.sum() >= 20:
        per_year[y] = {"n": int(m.sum()), "med": round(float(np.median(resid[m])), 3),
                       "mad": round(float(np.median(np.abs(resid[m] - np.median(resid[m]))) * 1.4826), 3)}
meds = [v["med"] for v in per_year.values()]
pre = [v["med"] for y, v in per_year.items() if y < 2017]
post = [v["med"] for y, v in per_year.items() if y >= 2018]

# (b) matched-exposure L0.5 vs L1 spot check (2010-06-15 C3)
def chain_mags(fr):
    if not (L.fit_frame(fr) or fr.rot_deg is not None) or fr.zp_r is None:
        return {}
    dyr = (fr.mjd_mid - 48348.5625) / 365.25
    ra = L._cat["ra"] + L._cat["pmra_masyr"] * dyr / 3.6e6 / np.cos(np.radians(L._cat["dec"]))
    dec = L._cat["dec"] + L._cat["pmde_masyr"] * dyr / 3.6e6
    sep = np.hypot((ra - fr.sun_ra) * np.cos(np.radians(dec)), dec - fr.sun_dec)
    m = (sep < 8.0) & (L._cat["vmag"] > 5.0) & (L._cat["vmag"] < 9.0)
    out = {}
    for i in np.nonzero(m)[0]:
        p = L.forced_photometry(fr, float(ra[i]), float(dec[i]))
        if p["flux"] > 0 and np.isfinite(p["zp"]) and p["flux"] / p["err"] > 8:
            out[i] = p["zp"] - 2.5 * np.log10(p["flux"])
    return out

fr05 = L.load_frame(RECON / "c3_100615_32227512.fts")
m05 = chain_mags(fr05)
blob = gzip.decompress((RECON / "c3_L1_35227512.fts.gz").read_bytes())
from astropy.io import fits
hdu = fits.open(io.BytesIO(blob))[0]
import tempfile
with tempfile.NamedTemporaryFile(suffix=".fts", delete=False) as tf:
    fits.HDUList([hdu]).writeto(tf.name, overwrite=True)
    fr1 = L.load_frame(Path(tf.name))
m1 = chain_mags(fr1)
common = sorted(set(m05) & set(m1))
d = np.array([m05[i] - m1[i] for i in common])

out = {"per_year_c3": per_year,
       "year_to_year_scatter_of_medians": round(float(np.std(meds)), 3),
       "pre2017_vs_post2018_offset": round(float(np.median(post) - np.median(pre)), 3) if pre and post else None,
       "l1_spot_check": {"n_stars": len(common),
                         "median_dmag_L05_minus_L1": round(float(np.median(d)), 3) if len(d) else None,
                         "mad": round(float(np.median(np.abs(d - np.median(d))) * 1.4826), 3) if len(d) else None,
                         "gate_0.1mag": bool(len(d) and abs(np.median(d)) <= 0.1)}}
OUT.write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "per_year_c3"}, indent=1))
print("years:", {y: v["med"] for y, v in per_year.items()})
