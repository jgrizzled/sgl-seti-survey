"""Dev stage v1 (hypotheses.md D7/D8 + amendment v1.1): the dev units
(gj-908 A/B 0.1 AU; A gj-1276 1.2 Rsun demoted), the D1 forced-dev window
(B van-maanen 2020-04-02 at three rungs as pseudo-units), the two positive
controls run through the identical chain, and the machinery checks.
Writes results/dev_v1.json and results/dev_v1.md."""
from __future__ import annotations

import json

import numpy as np
from astropy.io import fits
from astropy.time import Time

import hlib as H
import recon_scan as R

CONTROLS = {
    "ctrl-grb130427a": {"ra": 173.136, "dec": 27.699, "t0": "2013-04-27T07:47:06", "half_d": 0.32, "rung_like": "1.2Rsun",
                        "expect": ["S_event_L", "S_event_H", "S_burst"]},
    "ctrl-3c454.3": {"ra": 343.491, "dec": 16.148, "t0": "2010-11-19T12:00:00", "half_d": 2.5, "rung_like": "0.1AU",
                     "expect": ["S_event_L"]},
}


def fake_windows(t0_utc, half_d):
    t = Time(t0_utc).mjd
    return [{"event_id": "ctrl", "t_ca_utc": t0_utc, "t_ca_mjd": t, "mjd0": t - half_d, "mjd1": t + half_d,
             "dur_d": 2 * half_d, "b_rsun": 0.0}]


def summarize(a):
    return {"stats": a["stats"], "n_searchable": a["n_searchable"],
            "diag": {k: v for k, v in a["diag"].items() if not k.startswith("stack_ens")},
            "windows": [{k: w[k] for k in ("t_ca_utc", "searchable", "live_s", "moon_excluded_s", "n_pw", "expo_L", "expo_H",
                                          "n_L", "lam_L", "n_H", "lam_H", "n_U", "rate_U")} for w in a["windows"]]}


def main():
    events = R.load_events()
    T = H.FREEZE["T"]
    out = {"utc": Time.now().isot, "T": T, "dev_units": {}, "forced_dev_window": {}, "controls": {}, "checks": {}}
    md = ["# Dev stage v1 — 2026-09-07", "", f"T = {T:.3f} (145 trials, FWER 0.05).", ""]

    # 1. dev units
    for u in H.all_units(events):
        if not u["dev"]:
            continue
        ph = H.photons(u["key"], u["ra"], u["dec"])
        ws = H.unit_windows(u["channel"], u["target"], u["rung"], events)
        a = H.unit_analysis(u["channel"], u["target"], u["rung"], ph, u["ra"], u["dec"], ws)
        name = f"{u['channel']} {u['target']} {u['rung']}"
        out["dev_units"][name] = summarize(a)
        md.append(f"- dev unit **{name}**: " + ", ".join(f"{k} {v:.2f}" for k, v in a["stats"].items())
                  + f" ({a['n_searchable']} searchable windows; exceed T: {[k for k, v in a['stats'].items() if v > T]})")
        print(f"[dev] {name}: {a['stats']}", flush=True)

    # 2. the D1 forced-dev window as pseudo-units
    pos = R.positions(events)
    ra, de, _ = pos[("B", "van-maanen")]
    ph = H.photons("B-van-maanen", ra, de)
    for rung in H.FREEZE["rungs"]:
        ws = [w for w in H.unit_windows("B", "van-maanen", rung, events) if w["t_ca_utc"].startswith("2020-04-02")]
        a = H.unit_analysis("B", "van-maanen", rung, ph, ra, de, ws)
        out["forced_dev_window"][rung] = summarize(a)
        md.append(f"- D1 window B van-maanen 2020-04-02 {rung}: " + ", ".join(f"{k} {v:.2f}" for k, v in a["stats"].items())
                  + f"; window n_L {a['windows'][0]['n_L']} vs λ {a['windows'][0]['lam_L']:.2f}, n_H {a['windows'][0]['n_H']} vs λ {a['windows'][0]['lam_H']:.3f}")
        print(f"[dev] D1 {rung}: {a['stats']}", flush=True)

    # 3. positive controls
    for key, c in CONTROLS.items():
        ph = H.photons(key, c["ra"], c["dec"])
        ws = fake_windows(c["t0"], c["half_d"])
        a = H.unit_analysis("C", key, c["rung_like"], ph, c["ra"], c["dec"], ws)
        passed = all(a["stats"][s] > T for s in c["expect"])
        out["controls"][key] = dict(summarize(a), expect=c["expect"], passed=passed)
        md.append(f"- control **{key}**: " + ", ".join(f"{k} {v:.2f}" for k, v in a["stats"].items())
                  + f"; expected {c['expect']} > T → **{'PASS' if passed else 'FAIL'}**; window n_L {a['windows'][0]['n_L']} vs λ {a['windows'][0]['lam_L']:.2f}")
        print(f"[dev] control {key}: {a['stats']} pass={passed}", flush=True)

    # 4. machinery checks
    # (a) gate census over all positions fetched so far
    census = {}
    for (ch, tid), (r0, d0, n) in sorted(pos.items()):
        key = f"{ch}-{tid}"
        try:
            p = H.photons(key, r0, d0)
        except FileNotFoundError:
            continue
        census[key] = {"raw": int(p["n_raw"][0]), "source_class_zenith": int(p["n_class"][0]), "gated": int(p["n_gated"][0])}
    out["checks"]["gate_census"] = census
    # (b) photon identity vs the recon data-server pull (van-maanen B, 2020-04-01 -> 04, 5 deg)
    recon = H.REPO / "runs" / "highenergy-crossings" / "recon" / "lat_query_L260907093219AD225ECD61_PH00.fits"
    if recon.exists():
        with fits.open(recon) as h:
            d = h["EVENTS"].data
            rt = np.sort(np.asarray(d["TIME"]))
        p = H.photons("B-van-maanen", ra, de)
        m0, m1 = Time("2020-04-01").mjd, Time("2020-04-04").mjd
        sel = (p["mjd"] >= m0) & (p["mjd"] < m1)
        ours = p["TIME"][sel]
        found = np.isin(np.round(ours, 3), np.round(rt, 3))
        out["checks"]["recon_identity"] = {"ours_gated_in_span": int(sel.sum()), "found_in_recon_file": int(found.sum()),
                                           "recon_file_photons": int(len(rt))}
    # (c) Moon exclusion effect
    cov = json.loads((H.RES / "coverage_v1.json").read_text())
    moon_s = sum(w["moon_excluded_s"] for u in cov["units"] for w in u["windows"])
    live_s = sum(w["live_s"] for u in cov["units"] for w in u["windows"])
    n_moon_win = sum(1 for u in cov["units"] for w in u["windows"] if w["moon_excluded_s"] > 0)
    out["checks"]["moon"] = {"excluded_s_all_windows": moon_s, "gated_live_s_all_windows": live_s, "windows_affected": n_moon_win,
                             "windows_total": sum(len(u["windows"]) for u in cov["units"])}
    # (d) lane containment / effective area for the record
    out["checks"]["irf"] = {lane: {"containment_cos0.8": H.lane_containment(lane), "containment_cos0.5": H.lane_containment(lane, 0.5),
                                   "aeff_onaxis_m2": float(H.aeff_lane(lane)[1][-1]), "mean_E_MeV": float(H.lane_mean_energy_mev(lane))}
                            for lane in ("L", "H", "U")}
    md += ["", "## Machinery checks", "",
           f"- Gate census (raw → SOURCE+zenith → interval-gated): " + "; ".join(f"{k} {v['raw']}/{v['source_class_zenith']}/{v['gated']}" for k, v in census.items()),
           f"- Recon-file identity: {out['checks'].get('recon_identity')}",
           f"- Moon/Sun exclusion: {moon_s/1e3:.1f} ks excluded of {live_s/1e3:.0f} ks gated over {n_moon_win}/{out['checks']['moon']['windows_total']} windows",
           f"- IRF: " + "; ".join(f"{k} A_eff on-axis {v['aeff_onaxis_m2']:.3f} m², containment {v['containment_cos0.8']:.2f} (cosθ 0.8) / {v['containment_cos0.5']:.2f} (0.5), ⟨E⟩ {v['mean_E_MeV']:.0f} MeV" for k, v in out["checks"]["irf"].items())]
    (H.RES / "dev_v1.json").write_text(json.dumps(H.clean(out), indent=1))
    (H.RES / "dev_v1.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
