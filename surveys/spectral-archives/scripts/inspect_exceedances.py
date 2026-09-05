"""Manual-adjudication evidence for every exceedance of a family
(thresholds §4 ladder, human step). For each exceedance: the z profile
(+-12 px), the template and null-ensemble behaviour at the peak
(template depth, MAD, fraction of nulls with |z| > 3 there), the
observer-frame wavelength and the SPIRou OH sky model intensity there
(empirical NIR sky-line reference, observer frame), nearby strong
template features, and whether the same barycentric wavelength shows
S > T/2 in other in-window spectra of the unit or in any null.

usage: inspect_exceedances.py dev|confirmatory
writes results/<family>_exceedance_evidence_v1.json + prints a digest.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spec_lib import CFG, C_KMS, Ensemble, Grid, load_spectrum, normalise, stellar_line_mask, to_grid, gauss_fit, despike  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
RUN = REPO / "runs" / "spectral-archives" / "v1"
DATA = RUN / "data"; ENS = RUN / "ensembles"
fam = sys.argv[1]
res = json.load(open(HERE / "results" / f"{fam}_search_v1.json"))
units = json.load(open(HERE / "results" / "units_v1.json"))
import yaml
REG = yaml.safe_load(open(REPO / "registries" / "pilot_wise_2026.yaml"))["targets"]

# empirical OH sky reference from a SPIRou t file (observer frame, vacuum nm)
sky_ref = None
for f in sorted((DATA / "SPIRou").glob("*.fits"))[:20]:
    h = fits.open(f)
    if "OHLine" in [x.name for x in h]:
        W = h["WaveAB"].data.ravel(); OH = h["OHLine"].data.ravel()
        ok = np.isfinite(W) & np.isfinite(OH)
        o = np.argsort(W[ok]); sky_ref = (W[ok][o], OH[ok][o]); break


def oh_intensity(lam_obs_nm, halfwidth_nm=0.03):
    if sky_ref is None:
        return None
    w, oh = sky_ref
    sel = (w > lam_obs_nm - halfwidth_nm) & (w < lam_obs_nm + halfwidth_nm)
    if not sel.any():
        return None
    return {"max": float(np.nanmax(oh[sel])), "p99_all": float(np.nanpercentile(oh, 99)), "p90_all": float(np.nanpercentile(oh, 90))}


out = []
for u in res["units"]:
    if not u.get("exceedances"):
        continue
    tid, inst = u["target"], u["instrument"]
    grid = Grid(inst)
    E = np.load(ENS / f"{tid}_{inst}.npz")
    template, mad, pix_ok, extra = E["template"], E["mad"], E["pix_ok"], E["extra_mask"]
    NF, NE = E["F"], E["E"]
    # load in-window spectra
    specs = {}
    for s in u["spectra"]:
        if s.get("status") != "ok":
            continue
        sp = load_spectrum(inst, DATA / inst / s["file"])
        F, Er, M = to_grid(sp, grid); nf, ne, nm = normalise(F, Er, M)
        z = (nf - template) / np.fmax(mad, ne); z[~pix_ok | extra] = np.nan; z, _ = despike(z)
        specs[s["file"]] = {"z": z, "nf": nf, "berv": sp["meta"].get("berv"), "utc": s["utc"], "raw": F, "snr": sp["meta"].get("snr")}
    for e in u["exceedances"]:
        lam = e["evidence"]["lambda_bary_vac_nm"]; p = int(grid.idx(lam))
        fn = e["spectrum"]; sp = specs.get(fn)
        sl = slice(p - 12, p + 13)
        zn = (NF - template) / np.fmax(mad, NE)
        frac_null_hot = float(np.nanmean(np.abs(zn[:, p]) > 3)) if np.isfinite(zn[:, p]).any() else None
        berv = sp["berv"] if sp else None
        lam_obs = lam / (1 + (berv or 0) / C_KMS)
        rv = float(REG[tid]["state"]["astrometry"]["radial_velocity_km_s"])
        lam_star = lam / (1 + rv / C_KMS)
        others = {}
        for k, v in specs.items():
            if k == fn:
                continue
            others[k[:30]] = float(np.nanmax(v["z"][p - 2:p + 3])) if np.isfinite(v["z"][p - 2:p + 3]).any() else None
        ev = {"unit": u["unit_id"], "statistic": e["statistic"], "cell": e["cell"], "S": e["S"], "T": e["T"], "disposition_auto": e["disposition"],
              "lambda_bary_nm": lam, "lambda_obs_nm": lam_obs, "lambda_star_rest_nm": lam_star, "spectrum": fn, "utc": sp["utc"] if sp else None,
              "snr": sp["snr"] if sp else None, "berv": berv,
              "z_profile": [None if not np.isfinite(x) else round(float(x), 1) for x in (sp["z"][sl] if sp else [])],
              "template_profile": [None if not np.isfinite(x) else round(float(x), 3) for x in template[sl]],
              "mad_profile": [None if not np.isfinite(x) else round(float(x), 4) for x in mad[sl]],
              "norm_flux_profile": [None if not np.isfinite(x) else round(float(x), 3) for x in (sp["nf"][sl] if sp else [])],
              "frac_nulls_abs_z_gt3_at_peak": frac_null_hot,
              "template_depth_at_peak": float(template[p]) if np.isfinite(template[p]) else None,
              "template_min_within_5px": float(np.nanmin(template[p - 5:p + 6])),
              "oh_sky_model_at_obs_lambda": oh_intensity(lam_obs), "other_inwindow_max_z_pm2px": others,
              "ladder": e["evidence"]}
        out.append(ev)
        print(f"\n{u['unit_id']} {e['statistic']} {e['cell']} S={e['S']:.1f} T={e['T']:.1f} auto={e['disposition']}")
        print(f"  λ_bary {lam:.3f} nm, λ_obs {lam_obs:.3f}, λ_star-rest {lam_star:.3f}; spectrum {fn[:34]} {ev['utc']} snr {ev['snr']}")
        print(f"  z  : {ev['z_profile']}")
        print(f"  tpl: {ev['template_profile']}")
        print(f"  nf : {ev['norm_flux_profile']}")
        print(f"  frac nulls |z|>3 at peak {frac_null_hot}; template min ±5px {ev['template_min_within_5px']:.3f}; OH sky {ev['oh_sky_model_at_obs_lambda']}")
        print(f"  other in-window max z ±2px: {others}")
(HERE / "results" / f"{fam}_exceedance_evidence_v1.json").write_text(json.dumps(out, indent=1))
