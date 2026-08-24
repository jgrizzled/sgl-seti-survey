"""Channel-A per-target saturation cut (hypotheses.md freeze v1.0).

Port of surveys/ztf-crossings/scripts/saturation_cut.py to PS1 grizy.
Frozen per-band rule (hypotheses v1.0: derived warp point-source
saturation g 13.3 / r 13.3 / i 13.4 / z 12.5 / y 11.5 plus ~0.5 mag
conservative margin; to be verified against CELL.SATURATION at dev):

    excluded  m_est < E_b        E = {g 14.0, r 14.0, i 14.0, z 13.0, y 12.0}
    marginal  E_b <= m_est < E_b + 0.5   (searchable, flagged)
    ok        m_est >= E_b + 0.5

m_est per PS1 band from Gaia DR3 G and G-RP via a piecewise-linear
dwarf-sequence color table (declared accuracy +/-0.3 mag, absorbed in
the margin). Sources: CNS5 census snapshot (Gaia DR3 join, CNS5-id
fallback); Gaia-saturated classics and non-dwarf-sequence objects
(white dwarfs, unresolved late-M multiples) from a hand table. Objects
with no optical photometry because they are optically dark (late-T/Y
dwarfs) are 'ok' for saturation but flagged `optically_dark` - they
cannot saturate, and the search sensitivity statement for them is the
empty-field limit, not a contrast limit.

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
OUT = REPO / "surveys" / "ps1-crossings" / "results"
CENSUS = REPO / "targets" / "census" / "cns5_10pc_2026-08-20.csv"
REGISTRY = REPO / "registries" / "pilot_wise_2026.yaml"

BANDS = ("g", "r", "i", "z", "y")
EXCLUDE = {"g": 14.0, "r": 14.0, "i": 14.0, "z": 13.0, "y": 12.0}
MARGIN = 0.5

# dwarf-sequence colors vs G-RP (approx. F8..M7), +/-0.3 mag; g-r/r-G/r-i
# reused from the ZTF cut, r-z/r-y extended for PS1 from the same sequence
GRP = np.array([0.25, 0.45, 0.65, 0.85, 1.05, 1.25, 1.50])
G_R = np.array([0.20, 0.40, 0.80, 1.30, 1.45, 1.50, 1.70])   # g-r
R_G = np.array([0.00, 0.05, 0.15, 0.30, 0.45, 0.60, 0.80])   # r-G
R_I = np.array([0.00, 0.10, 0.25, 0.55, 0.95, 1.50, 2.20])   # r-i
R_Z = np.array([0.05, 0.13, 0.30, 0.65, 1.15, 1.80, 2.60])   # r-z
R_Y = np.array([0.08, 0.20, 0.45, 0.90, 1.45, 2.20, 3.10])   # r-y

# Gaia-saturated classics and non-dwarf-sequence objects: per-band
# estimates by hand (grizy). Classics are all far below every exclusion
# threshold, so +/-0.5 mag is irrelevant. V sources: Hipparcos / Bond
# et al. orbits as in registry; WD/late-M literature colors.
HAND = {  # tid: (g, r, i, z, y)
    "alpha-cen-a": (0.4, -0.2, -0.4, -0.5, -0.5),   # V=0.01 G2V
    "alpha-cen-b": (2.0, 1.0, 0.7, 0.6, 0.5),       # V=1.33 K1V
    "sirius-a": (-1.3, -1.4, -1.4, -1.3, -1.2),     # V=-1.46 A1V
    "sirius-b": (8.3, 8.6, 8.8, 9.0, 9.2),          # V=8.44 DA2 (blue)
    "procyon-a": (0.7, 0.2, 0.1, 0.1, 0.1),         # V=0.37 F5IV
    "procyon-b": (11.2, 10.8, 10.7, 10.7, 10.8),    # V=10.94 DQZ
    "fomalhaut": (1.2, 1.2, 1.3, 1.4, 1.4),         # V=1.16 A4V
    # WD: flat AB SED, dwarf-sequence transform invalid; PS1 DR2 mean
    # photometry (declared +/-0.5):
    "van-maanen": (12.4, 12.3, 12.4, 12.5, 12.5),   # DZ8 @4.3 pc
    # no census/Gaia photometry; literature-based (declared +/-0.5):
    "ez-aqr": (13.4, 11.1, 9.4, 9.2, 8.9),          # GJ 866 ABC, M5V colors
    "gj-11068": (16.9, 15.4, 13.2, 12.7, 12.3),     # CNS5 1819, ~M6.5 colors
    # G-RP = 1.59, beyond the table's declared range (clip biases r
    # bright by > 1 mag at M7); literature M7V colors (declared +/-0.5):
    "teegarden": (16.6, 14.2, 11.9, 10.9, 10.3),
}
OPTICALLY_DARK = {"wise-0855", "luhman16-a", "luhman16-b"}  # Y0 / L-T


def est_mags(G: float, RP: float):
    c = np.clip(G - RP, GRP[0], GRP[-1])
    r = G + np.interp(c, GRP, R_G)
    return {"g": float(r + np.interp(c, GRP, G_R)), "r": float(r),
            "i": float(r - np.interp(c, GRP, R_I)),
            "z": float(r - np.interp(c, GRP, R_Z)),
            "y": float(r - np.interp(c, GRP, R_Y))}


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
            mags = dict(zip(BANDS, HAND[tid]))
            source = "hand_table"
        elif row and row.get("Gmag") and row.get("RPmag"):
            mags = est_mags(float(row["Gmag"]), float(row["RPmag"]))
            # beyond the table range the transform is bright-biased,
            # i.e. conservative (may exclude a searchable target)
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
            n_cov = int((np.asarray(sub[f"n_{b}_pri"]) > 0).sum())
            rec[f"cov_rows_{b}"] = (n_cov if st[b] in ("ok", "marginal")
                                    else 0)
        out.append(rec)

    tab = Table(rows=out)
    tab.write(OUT / "saturation_cut_v1.ecsv", format="ascii.ecsv",
              overwrite=True)
    summary = {
        "rule": {"excluded_below": EXCLUDE, "marginal_width": MARGIN,
                 "transform_accuracy_mag": 0.3},
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
