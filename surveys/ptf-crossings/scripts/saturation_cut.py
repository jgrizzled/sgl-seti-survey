"""Channel-A per-target saturation cut (hypotheses.md freeze v1.0 §6).

Port of surveys/ps1-crossings/scripts/saturation_cut.py to PTF g/R.
Frozen per-band rule (hypotheses v1.0: 60 s P48 point-source
saturation estimates R 14.0 / g 14.5 plus a 0.5 mag conservative
margin; to be verified against a bright-star dmask bit-8 extent at
dev, amendment if off by > 0.5 mag):

    excluded  m_est < E_b        E = {R 14.5, g 15.0}
    marginal  E_b <= m_est < E_b + 0.5   (searchable, flagged)
    ok        m_est >= E_b + 0.5

m_est from Gaia DR3 G and G-RP via the PS1-cut piecewise-linear
dwarf-sequence table (declared accuracy +/-0.3 mag, absorbed in the
margin): PTF g = the SDSS-like g estimate directly; PTF Mould R via
the Jordi et al. (2006) Cousins transform R = r - 0.153(r-i) - 0.117
applied to the sequence r/i estimates. Hand-table targets reuse the
PS1 cut's per-band grizy values through the same R transform.
Optically dark targets (late-T/Y) cannot saturate; 'ok' but flagged.

Applies to every channel-A era target; also reports the surviving
covered-event counts by joining results/coverage_v1_events.ecsv.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import yaml
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "surveys" / "ptf-crossings" / "results"
CENSUS = REPO / "targets" / "census" / "cns5_10pc_2026-08-20.csv"
REGISTRY = REPO / "registries" / "pilot_wise_2026.yaml"

BANDS = ("g", "R")
EXCLUDE = {"R": 14.5, "g": 15.0}
MARGIN = 0.5

# dwarf-sequence colors vs G-RP (approx. F8..M7), +/-0.3 mag — the
# PS1/ZTF cut tables, verbatim
GRP = np.array([0.25, 0.45, 0.65, 0.85, 1.05, 1.25, 1.50])
G_R = np.array([0.20, 0.40, 0.80, 1.30, 1.45, 1.50, 1.70])   # g-r
R_G = np.array([0.00, 0.05, 0.15, 0.30, 0.45, 0.60, 0.80])   # r-G
R_I = np.array([0.00, 0.10, 0.25, 0.55, 0.95, 1.50, 2.20])   # r-i

# Non-dwarf-sequence objects: the PS1 cut's hand table (grizy),
# reduced here to (g, r, i) for the g / Mould-R estimates.
HAND = {  # tid: (g, r, i)
    "alpha-cen-a": (0.4, -0.2, -0.4),
    "alpha-cen-b": (2.0, 1.0, 0.7),
    "sirius-a": (-1.3, -1.4, -1.4),
    "sirius-b": (8.3, 8.6, 8.8),
    "procyon-a": (0.7, 0.2, 0.1),
    "procyon-b": (11.2, 10.8, 10.7),
    "fomalhaut": (1.2, 1.2, 1.3),
    "van-maanen": (12.4, 12.3, 12.4),
    "ez-aqr": (13.4, 11.1, 9.4),
    "gj-11068": (16.9, 15.4, 13.2),
    "teegarden": (16.6, 14.2, 11.9),
}
OPTICALLY_DARK = {"wise-0855", "luhman16-a", "luhman16-b"}


def mould_r(r: float, i: float) -> float:
    """Jordi et al. (2006) Cousins R from SDSS-like r, i (declared
    +/-0.3, absorbed in the margin; Mould R ~ Cousins R)."""
    return r - 0.153 * (r - i) - 0.117


def est_mags(G: float, RP: float):
    c = np.clip(G - RP, GRP[0], GRP[-1])
    r = G + np.interp(c, GRP, R_G)
    i = r - np.interp(c, GRP, R_I)
    return {"g": float(r + np.interp(c, GRP, G_R)),
            "R": float(mould_r(r, i))}


def status(band, m):
    if m is None or np.isnan(m):
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
    A = cov[cov["channel"] == "A"]
    era_targets = sorted(set(str(t) for t in A["target_id"]))

    out = []
    for tid in era_targets:
        t = reg[tid]
        ids = {str(i.get("catalog", "")): str(i["id"])
               for i in t.get("identifiers", [])}
        gaia = next((v for k, v in ids.items() if k.startswith("Gaia")), None)
        row = rows_by_gaia.get(gaia)
        source = "gaia_join" if row else None
        if row is None:
            cns5 = ids.get("CNS5")
            if cns5:
                row = (rows_by_cns5.get(cns5)
                       or rows_by_cns5.get(cns5.split("-")[-1]))
                source = "cns5_join" if row else None
        mags = {b: None for b in BANDS}
        grp_clipped = False
        if tid in HAND:
            g, r, i = HAND[tid]
            mags = {"g": g, "R": mould_r(r, i)}
            source = "hand_table"
        elif row and row.get("Gmag") and row.get("RPmag"):
            mags = est_mags(float(row["Gmag"]), float(row["RPmag"]))
            grp_clipped = (float(row["Gmag"]) - float(row["RPmag"])
                           > GRP[-1])
        elif tid in OPTICALLY_DARK:
            source = "optically_dark"
        else:
            source = source or "missing"
        dark = tid in OPTICALLY_DARK
        st = {b: ("ok" if dark else status(b, mags[b])) for b in BANDS}
        sub = A[np.asarray([str(x) == tid for x in A["target_id"]])]
        rec = {"target_id": tid, "source": source, "optically_dark": dark,
               "grp_clipped": grp_clipped}
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
    summary = {
        "rule": {"excluded_below": EXCLUDE, "marginal_width": MARGIN,
                 "transform_accuracy_mag": 0.3,
                 "R_transform": "Jordi 2006: R = r - 0.153(r-i) - 0.117"},
        "n_targets": len(out),
    }
    for b in BANDS:
        col = [row[f"status_{b}"] for row in out]
        summary[b] = {s: col.count(s) for s in
                      ("excluded", "marginal", "ok", "no_photometry")}
        summary[b]["surviving_covered_event_rows"] = int(
            sum(row[f"cov_rows_{b}"] for row in out))
    (OUT / "saturation_cut_v1.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    for row in out:
        print(f"{row['target_id']:16s} "
              + " ".join(f"{b}={row[f'{b}_est']:>6}" for b in BANDS)
              + " [" + "/".join(row[f"status_{b}"][:3] for b in BANDS)
              + f"] src={row['source']}")


if __name__ == "__main__":
    main()
