"""Channel-A per-target saturation cut (hypotheses.md freeze v1.0).

WISE single-frame saturation W1 ~ 8.1 / W2 ~ 7.0 / W3 ~ 3.8 /
W4 ~ -0.4 Vega. Frozen rule per band:

    excluded  m_est < E_b        E = {W1 8.6, W2 7.5, W3 4.3, W4 0.1}
    marginal  E_b <= m_est < E_b + 0.5   (searchable, flagged)
    ok        m_est >= E_b + 0.5

m_est: census W1mag (CatWISE-derived CNS5 column) where present, else
Ks - 0.1 (M-dwarf sequence, +/-0.3 declared); W2 = W1 - 0.2,
W3 = W2 - 0.1, W4 = W3 - 0.05. Non-sequence objects (T/Y dwarfs,
unresolved multiples without census rows) from a hand table declared
+/-0.5. Unlike the optical surveys, "optically dark" objects are
IR-bright ordinary targets here.

Applies to every elongation-gate target; joins
results/coverage_v1_events.ecsv for surviving covered-event counts.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import yaml
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "surveys" / "wise-crossings" / "results"
CENSUS = REPO / "targets" / "census" / "cns5_10pc_2026-08-20.csv"
REGISTRY = REPO / "registries" / "pilot_wise_2026.yaml"

BANDS = ("W1", "W2", "W3", "W4")
EXCLUDE = {"W1": 8.6, "W2": 7.5, "W3": 4.3, "W4": 0.1}
MARGIN = 0.5
COLOR = {"W2": -0.2, "W3": -0.3, "W4": -0.35}   # offsets from W1

# T/Y dwarfs and objects without usable census photometry (+/-0.5):
HAND = {  # tid: (W1, W2, W3, W4)
    "wise-0855": (16.6, 13.9, 12.0, 11.0),      # Y0, very red
    "luhman16-a": (9.8, 9.4, 8.5, 8.2),         # L7.5, resolved half
    "luhman16-b": (9.9, 9.5, 8.6, 8.3),         # T0.5
    "ez-aqr": (5.4, 5.2, 5.0, 4.9),             # GJ 866 ABC combined
    "gj-11068": (9.3, 9.1, 8.8, 8.7),           # ~M6.5 @6.8 pc
    # census rows without W1/Ks — all far below every threshold:
    "barnard-star": (4.2, 3.6, 3.3, 3.2),       # M4V, Ks 4.5
    "fomalhaut": (0.9, 0.9, 0.9, 0.2),          # A4V + debris disc
    "gj-783": (3.5, 3.4, 3.3, 3.2),             # K2.5V, Ks 3.7
    "gj65-b": (5.0, 4.7, 4.4, 4.3),             # UV Cet, blended pair
    "procyon-a": (0.0, -0.1, -0.2, -0.2),       # F5IV
    "procyon-b": (0.0, -0.1, -0.2, -0.2),       # inside Procyon A's PSF
}


def status(band, m):
    if m is None or (isinstance(m, float) and np.isnan(m)):
        return "no_photometry"
    e = EXCLUDE[band]
    return ("excluded" if m < e
            else "marginal" if m < e + MARGIN else "ok")


def main():
    reg = yaml.safe_load(open(REGISTRY))["targets"]
    rows_by_gaia, rows_by_cns5 = {}, {}
    for r in csv.DictReader(open(CENSUS)):
        if r.get("GaiaDR3"):
            rows_by_gaia[r["GaiaDR3"]] = r
        if r.get("CNS5"):
            rows_by_cns5[r["CNS5"].strip()] = r

    cov = Table.read(OUT / "coverage_v1_events.ecsv")
    q = cov[np.asarray([s == "queried" for s in cov["status"]])]
    era_targets = sorted(set(str(t) for t in q["target_id"]))

    out = []
    for tid in era_targets:
        t = reg[tid]
        ids = {str(i.get("catalog", "")): str(i["id"])
               for i in t.get("identifiers", [])}
        gaia = next((v for k, v in ids.items() if k.startswith("Gaia")),
                    None)
        row = rows_by_gaia.get(gaia)
        source = "gaia_join" if row else None
        if row is None:
            cns5 = ids.get("CNS5")
            if cns5:
                row = (rows_by_cns5.get(cns5)
                       or rows_by_cns5.get(cns5.split("-")[-1]))
                source = "cns5_join" if row else None
        w1 = None
        if tid in HAND:
            mags = dict(zip(BANDS, HAND[tid]))
            source = "hand_table"
        else:
            if row and row.get("W1mag"):
                w1 = float(row["W1mag"])
                source = (source or "join") + "_w1"
            elif row and row.get("Ksmag"):
                w1 = float(row["Ksmag"]) - 0.1
                source = (source or "join") + "_ks"
            else:
                source = source or "missing"
            mags = ({b: (w1 + COLOR.get(b, 0.0)) if w1 is not None
                     else None for b in BANDS})
        st = {b: status(b, mags[b]) for b in BANDS}
        sub = q[np.asarray([str(x) == tid for x in q["target_id"]])]
        rec = {"target_id": tid, "source": source}
        for b in BANDS:
            m = mags[b]
            rec[f"{b}_est"] = np.nan if m is None else round(m, 2)
            rec[f"status_{b}"] = st[b]
            n_cov = int((np.asarray(sub[f"n_{b}"]) > 0).sum())
            rec[f"cov_rows_{b}"] = (n_cov if st[b] in ("ok", "marginal")
                                    else 0)
        out.append(rec)

    tab = Table(rows=out)
    tab.write(OUT / "saturation_cut_v1.ecsv", format="ascii.ecsv",
              overwrite=True)
    summary = {"rule": {"excluded_below": EXCLUDE,
                        "marginal_width": MARGIN,
                        "transform_accuracy_mag": 0.3},
               "n_targets": len(out)}
    for b in BANDS:
        col = [r[f"status_{b}"] for r in out]
        summary[b] = {s: col.count(s) for s in
                      ("excluded", "marginal", "ok", "no_photometry")}
        summary[b]["surviving_covered_event_rows"] = int(
            sum(r[f"cov_rows_{b}"] for r in out))
    (OUT / "saturation_cut_v1.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    for r in out:
        print(f"{r['target_id']:16s} "
              + " ".join(f"{b}={r[f'{b}_est']:>6}" for b in BANDS)
              + " [" + "/".join(r[f"status_{b}"][:3] for b in BANDS)
              + f"] src={r['source']}")


if __name__ == "__main__":
    main()
