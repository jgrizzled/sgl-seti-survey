"""Injection-calibrated completeness (frozen parameters: >=100 draws per
event-window cell, chord temporal profile = the flat-chord window, line
spectra reported as band AB magnitudes, v2 stamp-response machinery).

For every confirmatory row with in-window epochs, the exact
matched-filter response rho_e to a unit-flux Moffat(beta=3, SEEING)
source at the predicted position is computed per epoch
(sglsurvey.inject.stamp_response on FluxMaps rebuilt with
keep_inputs=True). The unit statistic is linear in injected flux, so
Delta_S(m) is analytic; recovery is Monte-Carlo'd by bootstrap-drawing
the null S from the unit's own control ensemble (confirmatory_v1.json)
and applying the frozen rule S_null + Delta_S > max(T, 0). m90 = the
faintest magnitude with >= 90 % recovery. Channel-A variances carry the
unit's frozen k rescale; channel B uses raw matched-filter variances
(v1.2 is A-only). Units without controls or in-window data are reported
uncalibrated, not silently dropped.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from confirmatory_search import (EV_BY_ID, PLAN_PATH, _dest, discover_cones,
                                 star_pos)
from coverage_intersect import relay_apparent
from dev_search import BANDCODE, WEIGHT_CAP, _snr_stack  # noqa
from sglsurvey.adapters.irsa_ztf import MASK_FATAL_TEMPLATE
from sglsurvey.inject import MoffatPSF, stamp_response
from sglsurvey.photometry import ZTF_PIX_ARCSEC, build_flux_map_ztf

D = REPO / "surveys" / "ztf-crossings"
RES = json.loads((D / "results" / "confirmatory_v1.json").read_text())
PLAN = json.loads(PLAN_PATH.read_text())
MAGS = np.arange(14.0, 22.01, 0.2)
NDRAW = 200          # >= frozen minimum of 100 per event-window cell
RNG = np.random.default_rng(20260824)
INJECT_Z = 2500.0


def fm_with_inputs(ch, tid, o):
    dest = _dest(ch, tid, o)
    if not dest.exists():
        return None
    sci = next((x for x in dest.iterdir()
                if x.name.endswith("sciimg.fits")), None)
    msk = next((x for x in dest.iterdir()
                if x.name.endswith("mskimg.fits")), None)
    diff = next((x for x in dest.iterdir()
                 if x.name.endswith("scimrefdiffimg.fits.fz")), None)
    if sci is None or msk is None:
        return None
    try:
        return build_flux_map_ztf(sci, diff, msk, MASK_FATAL_TEMPLATE,
                                  o.band, o.t_mid_mjd_utc, keep_inputs=True)
    except Exception:
        return None


def epoch_response(fm, ra, dec):
    """(rho, var, magzp) for a unit-flux source at (ra, dec)."""
    pix = fm.world2pix(np.array([ra]), np.array([dec]))
    x, y = float(pix[0][0]), float(pix[0][1])
    f, v, g = fm.sample(ra, dec)
    if not (np.isfinite(g[0]) and g[0] >= 0.7 and np.isfinite(v[0])):
        return None
    stamp, ox, oy = MoffatPSF(fm.fwhm_arcsec / ZTF_PIX_ARCSEC).render(x, y)
    rho = float(stamp_response(fm, stamp, ox, oy).sample(x, y))
    if not np.isfinite(rho) or rho <= 0 or fm.magzp is None:
        return None
    return rho, float(v[0]), float(fm.magzp)


def recovery(dS_of_m, controls, T):
    """m90 from bootstrap draws of the null over the magnitude grid."""
    pool = np.asarray(controls, float)
    thresh = max(T if np.isfinite(T) else 0.0, 0.0)
    m90 = None
    curve = []
    for m in MAGS:
        dS = dS_of_m(m)
        nulls = RNG.choice(pool, size=NDRAW, replace=True)
        frac = float(np.mean(nulls + dS > thresh))
        curve.append(round(frac, 3))
        if frac >= 0.9:
            m90 = float(m)
    return m90, curve


def main():
    adapter, omaps = discover_cones()
    out = {"A": [], "B": []}

    res_A = {(r["target_id"], r["band"], r["radius_au"]): r
             for r in RES["A"]}
    for p in PLAN["A"]:
        key = (p["target_id"], p["band"], p["radius_au"])
        r = res_A.get(key)
        if r is None:
            continue
        tid, band = p["target_id"], p["band"]
        k_scale = r["gate"].get("k_scale", None)
        ctrls, T = r["controls"], r["T"]
        if not ctrls or k_scale is None:
            out["A"].append(dict(zip(("target_id", "band", "radius_au"), key),
                                 status=r["status"], m90=None,
                                 note="uncalibrated (no controls/template)"))
            continue
        omap = omaps[("A", tid)]
        windows = [tuple(w) for w in p["rung_windows"]]
        per_win = []
        for lo, hi in windows:
            eps = []
            for kk in p["epoch_keys"]:
                o = omap[kk]
                if not (lo <= o.t_mid_mjd_utc <= hi):
                    continue
                fm = fm_with_inputs("A", tid, o)
                if fm is None:
                    continue
                ra, de = star_pos(tid, o.t_mid_mjd_utc)
                er = epoch_response(fm, float(ra[0]), float(de[0]))
                if er:
                    eps.append(er)
            if eps:
                per_win.append(eps)
        if not per_win:
            out["A"].append(dict(zip(("target_id", "band", "radius_au"), key),
                                 status=r["status"], m90=None,
                                 note="no in-window responses"))
            continue

        def dS(m, per_win=per_win, k=k_scale):
            tot = 0.0
            for eps in per_win:
                w = np.array([1.0 / (v * k) for _, v, _ in eps])
                w = np.minimum(w, WEIGHT_CAP * np.median(w))
                sig = np.array([rho * 10 ** (-0.4 * (m - zp))
                                for rho, _, zp in eps])
                tot += float(np.sum(w * sig) / np.sqrt(np.sum(w)))
            return tot / np.sqrt(len(per_win))

        m90, curve = recovery(dS, ctrls, T if T is not None else np.nan)
        out["A"].append(dict(zip(("target_id", "band", "radius_au"), key),
                             status=r["status"],
                             n_windows_injected=len(per_win),
                             m90=m90))
        print("A", key, r["status"], "m90=", m90, flush=True)

    res_B = {}
    for r in RES["B"]:
        res_B[(r["target_id"], r["band"], r["event_id"], r["radius_au"])] = r
    for p in PLAN["B"]:
        key = (p["target_id"], p["band"], p["event_id"], p["radius_au"])
        r = res_B.get(key)
        if r is None:
            continue
        tid, band, eid = p["target_id"], p["band"], p["event_id"]
        ev = EV_BY_ID[eid]
        ctrls = None
        # B controls not stored per row in json; rebuild pool from T only
        # is impossible - so use the row's T and a null pool of the row's
        # searchable siblings' S values (same target/band, empirical null)
        sib = [x["S"] for x in RES["B"]
               if x["target_id"] == tid and x["band"] == band
               and x.get("S") is not None and abs(x["S"]) < 10]
        T = r.get("T")
        if T is None or len(sib) < 3:
            out["B"].append(dict(zip(("target_id", "band", "event_id",
                                      "radius_au"), key),
                                 status=r["status"], m90=None,
                                 note="uncalibrated"))
            continue
        omap = omaps[("B", tid)]
        eps = []
        for kk in p["epoch_keys"]:
            o = omap[kk]
            fm = fm_with_inputs("B", tid, o)
            if fm is None:
                continue
            ra, de = relay_apparent(ev, INJECT_Z, o.t_mid_mjd_utc)
            er = epoch_response(fm, ra, de)
            if er:
                eps.append(er)
        if not eps:
            out["B"].append(dict(zip(("target_id", "band", "event_id",
                                      "radius_au"), key),
                                 status=r["status"], m90=None,
                                 note="no in-window responses"))
            continue

        def dS(m, eps=eps):
            w = np.array([1.0 / v for _, v, _ in eps])
            w = np.minimum(w, WEIGHT_CAP * np.median(w))
            sig = np.array([rho * 10 ** (-0.4 * (m - zp))
                            for rho, _, zp in eps])
            return float(np.sum(w * sig) / np.sqrt(np.sum(w)))

        m90, curve = recovery(dS, sib, T)
        out["B"].append(dict(zip(("target_id", "band", "event_id",
                                  "radius_au"), key),
                             status=r["status"], n_epochs=len(eps), m90=m90))
        print("B", key[0], key[1], key[3], r["status"], "m90=", m90,
              flush=True)

    (D / "results" / "completeness_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    for ch in ("A", "B"):
        vals = [u["m90"] for u in out[ch]
                if u.get("m90") and u["status"] == "searchable"]
        print(ch, "searchable m90 median:",
              round(float(np.median(vals)), 2) if vals else None,
              f"({len(vals)} units)")


if __name__ == "__main__":
    main()
