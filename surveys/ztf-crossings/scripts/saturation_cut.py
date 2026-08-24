"""Channel-A per-target saturation cut (hypotheses.md freeze v1.0).

ZTF single-epoch images saturate at ~12.5 mag (ZSDS); within ~0.5 mag
of that the blended-channel photometric statistic (window-locked excess
at the ~1-2 % level) is unreliable anyway. Frozen rule, per band:

    excluded  m_est < 13.0     (saturation + 0.5 mag margin)
    marginal  13.0 <= m_est < 13.5   (searchable, flagged)
    ok        m_est >= 13.5

m_est per ZTF band from Gaia DR3 G and G-RP via a piecewise-linear
dwarf-sequence color table (declared accuracy +/-0.3 mag, absorbed in
the margin). Sources: CNS5 census snapshot (Gaia DR3 join, CNS5-id
fallback); seven Gaia-saturated classics from a hand table (Hipparcos /
standard photometry). Objects with no optical photometry because they
are optically dark (late-T/Y dwarfs) are 'ok' for saturation but flagged
`optically_dark` - they cannot saturate, and the search sensitivity
statement for them is the empty-field limit, not a contrast limit.

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
OUT = REPO / "surveys" / "ztf-crossings" / "results"
CENSUS = REPO / "targets" / "census" / "cns5_10pc_2026-08-20.csv"
REGISTRY = REPO / "registries" / "pilot_wise_2026.yaml"

EXCLUDE_MAG, MARGINAL_MAG = 13.0, 13.5

# dwarf-sequence colors vs G-RP (approx. F8..M7 + WD/blue end), +/-0.3 mag
GRP = np.array([0.25, 0.45, 0.65, 0.85, 1.05, 1.25, 1.50])
G_R = np.array([0.20, 0.40, 0.80, 1.30, 1.45, 1.50, 1.70])   # g-r
R_G = np.array([0.00, 0.05, 0.15, 0.30, 0.45, 0.60, 0.80])   # r-G
R_I = np.array([0.00, 0.10, 0.25, 0.55, 0.95, 1.50, 2.20])   # r-i

# Gaia-saturated classics: V and V->ZTF-band offsets folded into per-band
# estimates by hand (all far below the exclusion threshold, so +/-0.5 mag
# is irrelevant). V sources: Hipparcos / Bond et al. orbits as in registry.
HAND = {  # tid: (g, r, i)
    "alpha-cen-a": (0.4, -0.2, -0.4),   # V=0.01 G2V
    "alpha-cen-b": (2.0, 1.0, 0.7),     # V=1.33 K1V
    "sirius-a": (-1.3, -1.4, -1.4),     # V=-1.46 A1V
    "sirius-b": (8.3, 8.6, 8.8),        # V=8.44 DA2 (blue)
    "procyon-a": (0.7, 0.2, 0.1),       # V=0.37 F5IV
    "procyon-b": (11.2, 10.8, 10.7),    # V=10.94 DQZ
    "fomalhaut": (1.2, 1.2, 1.3),       # V=1.16 A4V
    # no census/Gaia photometry; literature-based (declared +/-0.5):
    "ez-aqr": (13.4, 11.1, 9.4),        # GJ 866 ABC combined V~12.2, M5V colors
    "gj-11068": (16.9, 15.4, 13.2),     # CNS5 1819: 2MASS J=10.63 @6.8pc, ~M6.5 colors
}
OPTICALLY_DARK = {"wise-0855", "luhman16-a", "luhman16-b"}  # Y0 / L-T


def est_mags(G: float, RP: float):
    c = np.clip(G - RP, GRP[0], GRP[-1])
    r = G + np.interp(c, GRP, R_G)
    g = r + np.interp(c, GRP, G_R)
    i = r - np.interp(c, GRP, R_I)
    return float(g), float(r), float(i)


def status(m):
    if m is None or np.isnan(m):
        return "no_photometry"
    return ("excluded" if m < EXCLUDE_MAG
            else "marginal" if m < MARGINAL_MAG else "ok")


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
                key = cns5.replace("gj-", "").replace("wise-0855", "")
                row = (rows_by_cns5.get(cns5)
                       or rows_by_cns5.get(cns5.split("-")[-1]))
                source = "cns5_join" if row else None
        g = r = i = None
        if tid in HAND:
            g, r, i = HAND[tid]
            source = "hand_table"
        elif row and row.get("Gmag") and row.get("RPmag"):
            g, r, i = est_mags(float(row["Gmag"]), float(row["RPmag"]))
        elif tid in OPTICALLY_DARK:
            source = "optically_dark"
        else:
            source = source or "missing"
        dark = tid in OPTICALLY_DARK
        st = {b: ("ok" if dark else status(m))
              for b, m in (("g", g), ("r", r), ("i", i))}
        sub = A[np.asarray([str(x) == tid for x in A["target_id"]])]
        surv = {}
        for band, col in (("g", "n_g_pri"), ("r", "n_r_pri"), ("i", "n_i_pri")):
            n_cov = int((np.asarray(sub[col]) > 0).sum())
            surv[band] = n_cov if st[band] in ("ok", "marginal") else 0
        out.append({
            "target_id": tid, "source": source,
            "g_est": np.nan if g is None else round(g, 2),
            "r_est": np.nan if r is None else round(r, 2),
            "i_est": np.nan if i is None else round(i, 2),
            "status_g": st["g"], "status_r": st["r"], "status_i": st["i"],
            "optically_dark": dark,
            "cov_rows_g": surv["g"], "cov_rows_r": surv["r"],
            "cov_rows_i": surv["i"],
        })

    tab = Table(rows=out)
    tab.write(OUT / "saturation_cut_v1.ecsv", format="ascii.ecsv",
              overwrite=True)
    summary = {
        "rule": {"excluded_below": EXCLUDE_MAG, "marginal_below": MARGINAL_MAG,
                 "transform_accuracy_mag": 0.3},
        "n_targets": len(out),
    }
    for b in ("g", "r", "i"):
        col = [row[f"status_{b}"] for row in out]
        summary[b] = {s: col.count(s) for s in
                      ("excluded", "marginal", "ok", "no_photometry")}
        summary[b]["surviving_covered_event_rows"] = int(
            sum(row[f"cov_rows_{b}"] for row in out))
    (OUT / "saturation_cut_v1.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    for row in out:
        print(f"{row['target_id']:16s} g={row['g_est']:>6} r={row['r_est']:>6} "
              f"[{row['status_g'][:3]}/{row['status_r'][:3]}/{row['status_i'][:3]}]"
              f" src={row['source']}")


if __name__ == "__main__":
    main()
