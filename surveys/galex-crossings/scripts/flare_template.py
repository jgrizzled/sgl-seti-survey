"""Dev step iv addendum: FRED morphology of the wolf-359 off-window
flares (2009-03-25 / 2010-02-23) for the frozen A-channel flare veto —
peak epoch, rise time, decay e-fold from 10 s binned light curves, and
the flare/quiescent contrast. Merged into dev_v1.json as
iv_flare_template."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from dev_stage import Dev, DEV_JSON, merge
from sglsurvey.adapters.mast_gphoton import galex_ms_to_iso
import galexlib as gl

BIN_S = 10.0


def analyse(dev, field, v0, v1):
    pos = dev.field_pos(field, (v0 + v1) / 2)
    usable = dev.usable_at([(v0, v1)], pos[0], pos[1])
    ph = dev.photons_at([(v0, v1)], "NUV", pos[0], pos[1], usable)
    t = (np.asarray(ph, float) - ph.min()) / 1000.0
    nbins = int(t.max() // BIN_S) + 1
    counts, edges = np.histogram(t, bins=nbins,
                                 range=(0, nbins * BIN_S))
    rate = counts / BIN_S
    quiesc = float(np.median(rate))
    ipk = int(np.argmax(rate))
    peak = float(rate[ipk])
    # rise: last bin before peak at < quiesc + 20% of (peak-quiesc)
    thr = quiesc + 0.2 * (peak - quiesc)
    i = ipk
    while i > 0 and rate[i - 1] > thr:
        i -= 1
    rise_s = (ipk - i + 1) * BIN_S
    # decay e-fold: first bin after peak below quiesc + (peak-q)/e
    thr_e = quiesc + (peak - quiesc) / np.e
    j = ipk
    while j < nbins - 1 and rate[j + 1] > thr_e:
        j += 1
    decay_s = (j - ipk + 1) * BIN_S
    return {"visit_utc": galex_ms_to_iso(v0),
            "n_photons": int(len(ph)),
            "quiescent_cts_s": quiesc, "peak_cts_s": peak,
            "contrast": peak / quiesc if quiesc > 0 else None,
            "peak_offset_s": ipk * BIN_S,
            "rise_s_le": rise_s, "decay_efold_s": decay_s,
            "fred_like": bool(decay_s > rise_s)}


def main():
    dev = Dev()
    field = "wolf-359_star"
    rows = []
    for v0, v1 in dev.visits(field):
        if v1 - v0 < 900_000:
            continue
        r = analyse(dev, field, v0, v1)
        rows.append(r)
        print(f"[iv-b] {r['visit_utc'][:16]}: quiesc "
              f"{r['quiescent_cts_s']:.2f}, peak {r['peak_cts_s']:.2f} "
              f"cts/s, rise ≤{r['rise_s_le']:.0f}s, decay e-fold "
              f"{r['decay_efold_s']:.0f}s, FRED-like "
              f"{r['fred_like']}", flush=True)
    merge("iv_flare_template", {"bin_s": BIN_S, "flares": rows})


if __name__ == "__main__":
    main()
