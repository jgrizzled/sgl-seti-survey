"""The gj-1276 recurrence test (joint-crossings plan.md §2, frozen
before this script ran).

Tests the PS1 retained-ambiguous exceedance (evt-344c12d32a8d,
i_AB ~ 22.8 night-consistent at the z = 550 AU track node, 2014) for
recurrence in ZTF's covered gj-1276 wide-rung windows (2019-2026),
using the retained ZTF confirmatory cutouts and the frozen ZTF
channel-B construction with z FIXED at 550 AU (one pre-registered
joint trial).

Arm R1: per (event, band in g/r) S_550 with the 8 designated ring
controls; joint J = sum S / sqrt(N); joint controls pair ring index
across event-bands; exceedance J > max(T_J, 0). Sensitivity: exact
stamp-response injection of a flat-Fnu source (g = r = 22.8 AB, the
PS1 amplitude) through the identical chain; bootstrap null from the
per-unit ring ensembles (the ZTF completeness convention); joint
recovery at 22.8 and joint m90.

Arm R2 (i band) and the static arm are decided by arithmetic /
citation in plan.md - no computation here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "surveys" / "ztf-crossings" / "scripts"))

from sglsurvey.adapters.base import ConeRegion, MjdRange
from sglsurvey.adapters.irsa_ztf import ZtfSciAdapter
from sglsurvey.snapshots import SnapshotStore

from coverage_intersect import relay_apparent, load_events  # ztf's
from dev_search import _snr_stack  # ztf's frozen stack (WEIGHT_CAP)
import completeness_v1 as ztfc      # fm_with_inputs, epoch_response

D = REPO / "surveys" / "joint-crossings"
ZD = REPO / "surveys" / "ztf-crossings"
RUN = REPO / "runs" / "joint-crossings"
Z_FIXED = 550.0
BANDCODE = {"g": "zg", "r": "zr"}
RING = [(20.0, 0.0), (-20.0, 0.0), (30.0, 0.0), (-30.0, 0.0),
        (40.0, 0.0), (-40.0, 0.0), (0.0, 30.0), (0.0, -30.0)]
ERA = MjdRange(58178.0, 61275.0)
GOOD = 0.7
EPOCH_CLIP_SNR = 20.0     # amendment v1.1: ZTF single-epoch clip
M_REF = 22.8
MAGS = np.arange(20.0, 24.01, 0.2)
NDRAW = 200
RNG = np.random.default_rng(20260826)


def nkey(o):
    return json.dumps(o.native_key, sort_keys=True)


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    (D / "results").mkdir(parents=True, exist_ok=True)
    plan = json.loads((ZD / "configs" / "conf_plan_v1.json").read_text())
    rows = [r for r in plan["B"] if r["target_id"] == "gj-1276"
            and float(r.get("radius_au", 0.1)) == 0.1
            and r.get("epoch_keys")]
    print(f"{len(rows)} gj-1276 wide-rung plan rows with epochs")

    evs = load_events()["B"]
    sub = evs[np.asarray([str(x) == "gj-1276"
                          for x in evs["target_id"]])]
    ev_by_id = {str(e["event_id"]): e for e in sub}

    # one discovery cone (ZTF coverage_intersect convention)
    ras = np.array([float(e["relay_icrs_ra_deg"]) for e in sub])
    des = np.array([float(e["relay_icrs_dec_deg"]) for e in sub])
    ra0, de0 = float(np.median(ras)), float(np.median(des))
    cosd = max(np.cos(np.radians(de0)), 0.05)
    radius = max(np.ptp(ras) * cosd, np.ptp(des)) / 2.0 + 0.15
    adapter = ZtfSciAdapter()
    store = SnapshotStore(RUN)
    omap = {nkey(o): o for o in adapter.discover(
        ConeRegion(ra0, de0, radius), ERA, store)}
    print(f"rediscovered {len(omap)} quadrant-exposures")

    units = []
    for r in rows:
        band, eid = r["band"], r["event_id"]
        ev = ev_by_id[eid]
        eps = []
        for k in r["epoch_keys"]:
            o = omap.get(k)
            if o is None or o.band != BANDCODE.get(band):
                continue
            fm = ztfc.fm_with_inputs("B", "gj-1276", o)
            if fm is None:
                continue
            eps.append((o, fm))
        if not eps:
            continue

        clipped = [0]

        def S_at(dx, dy):
            sm = []
            for o, fm in eps:
                ra, dec = relay_apparent(ev, Z_FIXED, o.t_mid_mjd_utc)
                ra += dx / 3600.0 / max(np.cos(np.radians(dec)), .05)
                dec += dy / 3600.0
                f, v, g = fm.sample(ra, dec)
                if not (np.isfinite(g[0]) and g[0] >= GOOD):
                    continue
                if (v[0] > 0 and abs(f[0]) / np.sqrt(v[0])
                        > EPOCH_CLIP_SNR):
                    clipped[0] += 1
                    continue
                sm.append((f[0], v[0]))
            if not sm:
                return np.nan
            return _snr_stack(sm)[0]

        S = S_at(0.0, 0.0)
        ctrls = [S_at(dx, dy) for dx, dy in RING]
        # injection responses at the z=550 node (same clip)
        resp = []
        for o, fm in eps:
            ra, dec = relay_apparent(ev, Z_FIXED, o.t_mid_mjd_utc)
            f, v, g = fm.sample(ra, dec)
            if (np.isfinite(g[0]) and g[0] >= GOOD and v[0] > 0
                    and abs(f[0]) / np.sqrt(v[0]) <= EPOCH_CLIP_SNR):
                r_ = ztfc.epoch_response(fm, ra, dec)
                if r_ is not None:
                    resp.append(r_)
        units.append({"event_id": eid, "band": band,
                      "t_ca_mjd": float(ev["t_ca_mjd"]),
                      "n_epochs": len(eps),
                      "S550": None if np.isnan(S) else round(S, 3),
                      "controls": [None if np.isnan(c) else round(c, 3)
                                   for c in ctrls],
                      "n_clipped_samples": clipped[0],
                      "_resp": resp})

    live = [u for u in units if u["S550"] is not None
            and all(c is not None for c in u["controls"])]
    N = len(live)
    J = float(sum(u["S550"] for u in live) / np.sqrt(N))
    Jc = [float(sum(u["controls"][c] for u in live) / np.sqrt(N))
          for c in range(8)]
    T_J = max(Jc)
    exceed = bool(J > max(T_J, 0.0))

    # joint injection: dJ(m) analytic; null J bootstrapped from the
    # per-unit ring ensembles (ZTF completeness convention)
    def dJ(m):
        tot = 0.0
        for u in live:
            sig, var = [], []
            for rho, v, zp in u["_resp"]:
                sig.append(rho * 10 ** (-0.4 * (m - zp)))
                var.append(v)
            if not sig:
                continue
            w = 1.0 / np.asarray(var)
            w = np.minimum(w, 20.0 * np.median(w))
            tot += float(np.sum(w * np.asarray(sig))
                         / np.sqrt(np.sum(w)))
        return tot / np.sqrt(N)

    thresh = max(T_J, 0.0)
    rec = np.zeros(len(MAGS))
    for _ in range(NDRAW):
        Jn = float(sum(u["controls"][int(RNG.integers(8))]
                       for u in live) / np.sqrt(N))
        for mi, m in enumerate(MAGS):
            if Jn + dJ(m) > thresh:
                rec[mi] += 1
    rec /= NDRAW
    rec_ref = float(np.interp(M_REF, MAGS, rec))
    m90 = None
    below = np.where(rec < 0.9)[0]
    if rec[0] >= 0.9 and below.size:
        i = below[0]
        fr = (rec[i - 1] - 0.9) / max(rec[i - 1] - rec[i], 1e-9)
        m90 = float(MAGS[i - 1] + fr * (MAGS[i] - MAGS[i - 1]))
    elif rec[0] >= 0.9:
        m90 = float(MAGS[-1])

    if exceed:
        verdict = "joint_retained_anomaly"
    elif rec_ref >= 0.9:
        verdict = "flat_SED_recurrence_refuted_90pct"
    else:
        verdict = "inconclusive_depth"
    out = {
        "plan": "surveys/joint-crossings/plan.md section 2 "
                "(frozen; amendment v1.1 single-epoch clip applied)",
        "pre_amendment_v1_0": {"J": -13.45, "T_J": 1.116,
                               "exceedance": False,
                               "joint_recovery_at_22p8": 0.535,
                               "note": "bad-epoch anomaly units "
                                       "included; see plan amendment"},
        "z_fixed_au": Z_FIXED, "reference_mag_ab": M_REF,
        "n_event_bands": N,
        "units": [{k: u[k] for k in u if k != "_resp"} for u in units],
        "J": round(J, 3), "J_controls": [round(c, 3) for c in Jc],
        "T_J": round(T_J, 3), "exceedance": exceed,
        "joint_recovery_at_22p8": round(rec_ref, 3),
        "joint_m90": None if m90 is None else round(m90, 2),
        "recovery_curve": {"mags": [float(m) for m in MAGS],
                           "recovery": [round(float(x), 3)
                                        for x in rec]},
        "verdict": verdict,
        "arm_R2_i_band": "untestable: 2 in-window i epochs archive-"
                         "wide, 2-epoch depth ~21.1 vs required 22.8 "
                         "(plan.md, decided by arithmetic)",
        "static_arm": "settled by the PS1 DR2 stack catalog: no "
                      "counterpart within 6 arcsec to ~23+ depth",
    }
    (D / "results" / "gj1276_recurrence_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("n_event_bands", "J", "T_J", "exceedance",
                       "joint_recovery_at_22p8", "joint_m90",
                       "verdict")}, indent=1))


if __name__ == "__main__":
    main()
