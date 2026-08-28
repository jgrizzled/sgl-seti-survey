"""DASCH positive controls (thresholds.md v1.0, dev stage).

A-chain: RY Cnc (Algol EB, P = 1.092943 d, the DASCH tutorial
object; ICRS 129.97759 +19.82191, SIMBAD/Gaia) through the identical
querycat -> lightcurve -> usable-row chain; success = the folded
lightcurve shows the eclipse (phase-binned faint excursions at the
eclipse phase, depth >= 0.5 mag) and the out-of-eclipse scatter is
~0.15 mag-class.

B-chain: (run separately once sb_ident yields an in-field asteroid;
see results/sbident notes.)
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import daschlib as dl        # noqa: E402
import search_lib as sl      # noqa: E402

D = Path(__file__).resolve().parents[1]
SNAP = dl.RUNS / "dev" / "control"
RA, DEC = 129.97759, 19.82191
PERIOD = 1.092943  # days (DASCH tutorial)


def main() -> int:
    resp = dl.post("dasch/dr7/querycat",
                   {"ra_deg": RA, "dec_deg": DEC, "radius_arcsec": 15,
                    "refcat": "apass"}, SNAP / "querycat_rycnc.json")
    ent = dl.rows_of(resp)
    if not ent:
        print("no querycat entry for RY Cnc"); return 1
    e = min(ent, key=lambda r: float(r.get("stdmag") or 99))
    lc = dl.post("dasch/dr7/lightcurve",
                 {"gsc_bin_index": int(e["gsc_bin_index"]),
                  "ref_number": int(e["ref_number"]),
                  "refcat": "apass"}, SNAP / "lightcurve_rycnc.json")
    rows = [r for r in dl.rows_of(lc)
            if r.get("magcal_magdep") not in ("", "99.0")
            and r.get("date_jd") and r.get("reject_flag") in ("", "0")]
    t = np.array([float(r["date_jd"]) for r in rows])
    m = np.array([float(r["magcal_magdep"]) for r in rows])
    ph = ((t - 2400000.0) / PERIOD) % 1.0
    med = np.median(m)
    mad = np.median(np.abs(m - med)) * 1.4826
    # An Algol eclipse occupies only a small phase fraction, so bin
    # medians dilute it. Test instead that the faint excursions
    # (m > med + 0.4) concentrate at a common phase and are deep.
    faint = m > med + 0.4
    phf = ph[faint]
    # modal phase via circular mean of the faint excursions
    ang = 2 * np.pi * phf
    mode = (np.arctan2(np.sin(ang).sum(), np.cos(ang).sum())
            / (2 * np.pi)) % 1.0
    dph = np.abs(((phf - mode + 0.5) % 1.0) - 0.5)
    concentration = float((dph <= 0.15).mean()) if len(phf) else 0.0
    in_ecl = faint & (np.abs(((ph - mode + 0.5) % 1.0) - 0.5) <= 0.15)
    out_mask = np.abs(((ph - mode + 0.5) % 1.0) - 0.5) > 0.2
    depth = float(np.median(m[in_ecl]) - np.median(m[out_mask])) \
        if in_ecl.any() else 0.0
    scatter = float(np.median(np.abs(m[out_mask]
                                     - np.median(m[out_mask]))) * 1.4826)
    # null concentration for uniform phases = 0.3 (window width)
    passed = bool(len(phf) >= 20 and concentration >= 0.5
                  and depth >= 0.5 and scatter <= 0.3)
    result = dict(
        generated_utc=datetime.now(timezone.utc).isoformat(
            timespec="seconds"),
        object="RY Cnc", period_days=PERIOD,
        querycat_stdmag=float(e["stdmag"]), n_usable_rows=len(rows),
        global_median_B=float(med), global_mad=float(mad),
        n_faint_excursions=int(faint.sum()),
        eclipse_phase_mode=float(mode),
        faint_phase_concentration=concentration,
        concentration_uniform_null=0.3,
        eclipse_depth_mag=depth,
        out_of_eclipse_scatter_mag=scatter,
        passed=passed,
        snapshots={p.name: dl.sha256_file(p)
                   for p in SNAP.glob("*rycnc.json")})
    (D / "results" / "control_rycnc_v1.json").write_text(
        json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: result[k] for k in
                      ("querycat_stdmag", "n_usable_rows",
                       "n_faint_excursions", "faint_phase_concentration",
                       "eclipse_depth_mag", "out_of_eclipse_scatter_mag",
                       "passed")}, indent=1))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    sys.exit(main())
