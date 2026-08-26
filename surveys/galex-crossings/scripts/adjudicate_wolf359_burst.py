"""Adjudication of the single confirmatory exceedance: wolf-359 A NUV
S_burst 4.679 > T 4.127 (2007-03-03 event, 97 s visit). Frozen veto
ladder order: (1) stellar-flare veto = FRED morphology + the two-band
discriminator (FUV live on this unit); (3) strict re-run (boresight
<= 25') and detector-fixed clustering. Writes
results/adjudication_wolf359_v1.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import galexlib as gl
from confirmatory_search import Conf, FIELDS
from sglsurvey.adapters.mast_gphoton import galex_ms_to_iso

OUT = REPO / "surveys" / "galex-crossings" / "results"


def burst_window(times_ms, t_live_s):
    """Locate the (width, start, count, z) that maximizes S_burst."""
    t = np.asarray(times_ms, float) / 1000.0
    lam = len(t) / t_live_s
    best = None
    for w in gl.BURST_WIDTHS_S:
        starts = np.arange(t.min() - w, t.max() + w / 2.0, w / 2.0)
        left = np.searchsorted(t, starts)
        right = np.searchsorted(t, starts + w)
        counts = right - left
        i = int(np.argmax(counts))
        mu = lam * w
        z = (counts[i] - mu) / np.sqrt(max(mu, 0.5))
        if best is None or z > best["z"]:
            best = {"width_s": w, "start_s": float(starts[i]),
                    "count": int(counts[i]), "mu": float(mu),
                    "z": float(z)}
    return best


def main():
    conf = Conf()
    field = "wolf-359_star"
    ev = conf.event_row(field, "2007-03-03")
    frags = conf.inwindow_visits(field, ev)
    ra, dec = conf.unit_positions(field, ev, frags[0][0])[0]
    usable = conf.usable_at(frags, ra, dec)

    res = {"unit": "wolf-359_A_NUV", "event": "2007-03-03",
           "s_burst": 4.679, "threshold": 4.127}

    # NUV series and the maximizing burst window
    live_n = sum(1 for b in usable.values() if "NUV" in b)
    ph_n = conf.photons_at(frags, "NUV", ra, dec, usable)
    bw = burst_window(ph_n, live_n)
    res["burst"] = bw
    t0 = ph_n.min()
    res["visit_utc"] = galex_ms_to_iso(frags[0][0])

    # light curves, 5 s bins, both bands
    live_f = sum(1 for b in usable.values() if "FUV" in b)
    ph_f = conf.photons_at(frags, "FUV", ra, dec, usable)
    tn = (np.asarray(ph_n, float) - t0) / 1000.0
    tf = (np.asarray(ph_f, float) - t0) / 1000.0
    nb = int(max(tn.max(), tf.max() if len(tf) else 0) // 5) + 1
    lc_n, _ = np.histogram(tn, bins=nb, range=(0, nb * 5))
    lc_f, _ = np.histogram(tf, bins=nb, range=(0, nb * 5))
    res["lc_5s"] = {"nuv": lc_n.tolist(), "fuv": lc_f.tolist()}

    # two-band discriminator: FUV counts inside the NUV burst window
    w0 = bw["start_s"] - (ph_n.min() - t0) / 1000.0 + (t0 - t0) / 1e3
    ws = bw["start_s"] * 1000.0  # absolute ms
    in_f = int(np.sum((np.asarray(ph_f, float) >= ws)
                      & (np.asarray(ph_f, float) < ws
                         + bw["width_s"] * 1000.0)))
    in_n = bw["count"]
    base_f = len(ph_f) / max(live_f, 1) * bw["width_s"]
    from scipy.stats import poisson
    p_f = float(poisson.sf(in_f - 1, base_f)) if in_f > 0 else 1.0
    res["two_band"] = {
        "fuv_in_window": in_f, "fuv_expected": float(base_f),
        "fuv_excess_p": p_f, "fuv_total": int(len(ph_f)),
        "fuv_live_s": live_f}

    # morphology around the burst: rise/decay in 5 s bins
    ipk = int((bw["start_s"] - (t0 - t0) / 1e3) // 5)
    ipk = int(np.argmax(lc_n))
    quiesc = float(np.median(lc_n) / 5.0)
    peak = float(lc_n[ipk] / 5.0)
    res["morphology"] = {
        "quiescent_cts_s": quiesc, "peak_bin_cts_s": peak,
        "contrast": peak / quiesc if quiesc > 0 else None,
        "peak_bin_index": ipk, "n_bins": nb,
        "dev_flare_templates": "rise <=10-20 s, decay e-fold "
                               "30-40 s, contrast 3-21.6x"}

    # strict re-run: boresight <= 25'
    strict = {}
    for t_ms, ra0, de0, band, flag in conf.aspect_frags(frags):
        s = int(gl.sec_id(t_ms))
        if s in strict or flag % 2 != 0:
            continue
        if gl.ang_arcmin(ra0, de0, ra, dec) > 25.0:
            continue
        strict[s] = band
    live_s = sum(1 for b in strict.values() if "NUV" in b)
    ph_s = conf.photons_at(frags, "NUV", ra, dec, strict)
    res["strict_rerun"] = {
        "live_s": live_s,
        "s_burst": gl.s_burst(ph_s, live_s) if live_s else None}

    # detector-fixed check: xi/eta of burst-window photons
    xi_rows = []
    for v0, v1 in frags:
        sql = (f"select time, xi, eta from NUVPhotonsV where time >= "
               f"{int(v0)} and time < {int(v1) + 1000} and ra between "
               f"{ra - 11/3600} and {ra + 11/3600} and dec between "
               f"{dec - 11/3600} and {dec + 11/3600}")
        _, rows = conf.client.query(sql, conf.store)
        xi_rows.extend((int(r[0]), float(r[1]), float(r[2]))
                       for r in rows)
    inb = [(x, e) for t, x, e in xi_rows
           if ws <= t < ws + bw["width_s"] * 1000.0]
    outb = [(x, e) for t, x, e in xi_rows
            if not (ws <= t < ws + bw["width_s"] * 1000.0)]
    if inb and outb:
        mi = np.mean(np.array(inb), axis=0)
        mo = np.mean(np.array(outb), axis=0)
        res["detector_fixed"] = {
            "burst_xi_eta_mean": [float(mi[0]), float(mi[1])],
            "other_xi_eta_mean": [float(mo[0]), float(mo[1])],
            "n_in": len(inb), "n_out": len(outb)}

    (OUT / "adjudication_wolf359_v1.json").write_text(
        json.dumps(res, indent=2, default=float) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "lc_5s"},
                     indent=2, default=float))
    print("NUV lc (5s bins):", lc_n.tolist())
    print("FUV lc (5s bins):", lc_f.tolist())


if __name__ == "__main__":
    main()
