"""Completeness supplement v1: (a) gj-1276 B NUV 2010-03-03 event
(the 1,637 s deep visit) — S_burst / S_period injections vs the
unit's locked thresholds (its S_rate lane is D1a forced_dev: a
labeled reference grid is run, no threshold claim); (b) wolf-359 A
NUV extended persistent grid (the variability-inflated threshold
exceeded the v1 grid). Merged into completeness_v1.json."""
import json, sys
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(Path(__file__).resolve().parent))
import galexlib as gl
from confirmatory_search import Conf, SEED, OUT

conf = Conf()
res = json.load(open(OUT / "confirmatory_v1.json"))["results"]
comp = json.load(open(OUT / "completeness_v1.json"))
rng = np.random.default_rng(SEED + 99)

def setup(field, date, band):
    ev = conf.event_row(field, date)
    frags = conf.inwindow_visits(field, ev)
    ra, dec = conf.unit_positions(field, ev, frags[0][0])[0]
    usable = conf.usable_at(frags, ra, dec)
    live = sum(1 for b in usable.values() if band in b)
    stamps = np.array(sorted(s * 1000 + 995 for s, b in usable.items()
                             if band in b), dtype=np.int64)
    span = (stamps[-1] + 1000 - stamps[0]) / 1000.0
    ph = conf.photons_at(frags, band, ra, dec, usable)
    return live, stamps, span, ph

# (a) gj-1276 B NUV 2010 event
live, stamps, span, ph_u = setup("gj-1276_antipode", "2010-03-03", "NUV")
T = res["gj-1276_B_NUV"]["thresholds"]
rows = {"event": "2010-03-03", "live_s": live, "persistent_reference_d1a": [],
        "pulse": [], "train": []}
def rec(times, stat, thr, joff=0):
    m = gl.merge_series(ph_u, times)
    if stat == "S_rate": return gl.s_rate(m, live) > max(thr, 0)
    if stat == "S_burst": return gl.s_burst(m, live) > max(thr, 0)
    sp = gl.s_period(m, span, jitter_rng=np.random.default_rng(SEED + 6000 + joff))
    return (sp or 0.0) > max(thr, 0)
for rate in (0.02, 0.04, 0.07, 0.12, 0.2):
    n = sum(bool(rec(gl.inject_persistent(rng, rate, stamps), "S_rate", T["S_rate"])) for _ in range(25))
    rows["persistent_reference_d1a"].append({"rate_cts_s": rate, "recovered": n, "of": 25})
    print("2010 persist(ref)", rate, n, flush=True)
for w in (0.05, 0.5):
    for n_ph in (3, 4, 5, 6, 8, 12):
        n = sum(bool(rec(gl.inject_pulse(rng, n_ph, w, stamps), "S_burst", T["S_burst"])) for _ in range(25))
        rows["pulse"].append({"width_s": w, "n_photons": n_ph, "recovered": n, "of": 25})
    print("2010 pulse w", w, "done", flush=True)
for p in (0.05, 0.5, 5.0, 50.0):
    for n_tot in (15, 30, 60, 120):
        frac = n_tot * p / span
        n = 0
        for j in range(5):
            drift = rng.uniform(-gl.ORBIT_DRIFT_RATE, gl.ORBIT_DRIFT_RATE)
            n += bool(rec(gl.inject_train(rng, p, frac, 0.1, stamps, drift_rate=drift), "S_period", T["S_period"], j))
        rows["train"].append({"period_s": p, "n_injected": n_tot, "recovered": n, "of": 5})
    print("2010 train P", p, "done", flush=True)
comp["gj-1276_B_NUV_2010_supplement"] = rows

# (b) wolf-359 extended persistent grid
live, stamps, span, ph_u = setup("wolf-359_star", "2007-03-03", "NUV")
T = res["wolf-359_A_NUV"]["thresholds"]
ext = []
for rate in (1.5, 2.0, 2.5, 3.0, 4.0):
    n = sum(bool(rec(gl.inject_persistent(rng, rate, stamps), "S_rate", T["S_rate"])) for _ in range(25))
    ext.append({"rate_cts_s": rate, "recovered": n, "of": 25})
    print("wolf persist", rate, n, flush=True)
comp["wolf-359_A_NUV"]["persistent"].extend(ext)
(OUT / "completeness_v1.json").write_text(json.dumps(comp, indent=2, default=float) + "\n")
print("supplement merged")
