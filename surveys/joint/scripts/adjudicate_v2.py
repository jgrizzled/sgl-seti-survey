"""Stage-7 adjudication of cells retained by joint_ps1_ztf_mugrid.py:
(1) direct per-band ZTF forced photometry along the fitted (z, mu, T0)
track for BOTH roles (Rx/Tx loci coincide at the antipode: a real source
must appear in both); (2) catalogued sources (PS1 DR2 mean, with PM from
Gaia DR3 via VizieR) within 5" of the track at the epoch quantiles;
(3) the control-maximum margin. Writes runs/joint/ps1_ztf_v2/adjudication.json
and candidate_adjudicated.jsonl.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import requests

from sglseti import load_target_registry

from sglsurvey.geometry import GeometryContext
from sglsurvey.records import Candidate, append_records, read_records
from sglsurvey.vetting import load_ps1_mean, track_position

sys.path.insert(0, str(Path(__file__).resolve().parent))
from marginal_ztf_test import run_marginal_tests  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "runs" / "joint" / "ps1_ztf_v2"
sys.path.insert(0, str(REPO / "surveys" / "panstarrs" / "scripts"))
from ps1_corridors import CORRIDOR_OF  # noqa: E402


def gaia_near(ra, dec, r_arcsec=5.0):
    r = requests.get("https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
                     params={"-source": "I/355/gaiadr3", "-c": f"{ra:.6f} {dec:.6f}",
                             "-c.rs": f"{r_arcsec}", "-out": "RA_ICRS,DE_ICRS,Gmag,pmRA,pmDE,Plx",
                             "-out.add": "_r"}, timeout=60)
    rows = [l.split("\t") for l in r.text.splitlines() if l and not l.startswith("#")]
    return [dict(zip(rows[0], [x.strip() for x in v])) for v in rows[3:] if len(v) == len(rows[0])]


def main():
    registry = load_target_registry(REPO / "registries" / "pilot_wise_2026.yaml")
    ctx = GeometryContext.ps1_default()
    tr = json.load(open(OUT / "threshold_report.json"))
    run_id = read_records(OUT / "records" / "analysis_run.jsonl")[0]["analysis_run_id"]
    cands = [c for c in read_records(OUT / "records" / "candidate.jsonl") if c["status"] == "retained"]
    t0 = None
    cells = {}
    for c in cands:
        key = f"{c['endpoint_id']}/{c['role']}/{c['extra']['band']}"
        t0 = c["fitted_residual_motion"]["t0_mjd"]
        mc = c["model_comparison"]
        cells[key] = {"cell": {"z_au": c["fitted_z_au"], "mu": c["fitted_residual_motion"]["mu_arcsec_yr"]},
                      "S_real": mc["real_max_S"], "T": mc["threshold_8_controls"]}
        # other role, same cell
        e, role, band = key.split("/")
        other = "tx" if role == "rx" else "rx"
        cells[f"{e}/{other}/{band}"] = dict(cells[key], other_role_of=key)
    run_marginal_tests(OUT, run_id, cells=cells, t0_mjd=t0, out_name="direct_ztf_tests.json")
    direct = json.load(open(OUT / "direct_ztf_tests.json"))["cells"]
    report, sup = {}, []
    for c in cands:
        key = f"{c['endpoint_id']}/{c['role']}/{c['extra']['band']}"
        e, role, band = key.split("/")
        r = tr[key]
        z, mu = r["real_max_z"], r["real_max_mu"]
        cat = load_ps1_mean(REPO / "runs" / "panstarrs" / "screen_v1", CORRIDOR_OF[e], band)
        mjd = np.array(sorted([t for t in (55300, 56500, 58500, 59800, 61000)]))
        near = []
        for t in mjd:
            ra, dec = track_position(ctx, registry[e], role, z, mu, t0, float(t))
            n = cat.nearest(ra, dec, min_ndet=1)
            g = gaia_near(ra, dec)
            near.append({"mjd": float(t), "ps1_nearest": n,
                         "gaia_within_5arcsec": [{k: v.get(k) for k in ("_r", "Gmag", "pmRA", "pmDE", "Plx")} for v in g]})
        other = "tx" if role == "rx" else "rx"
        d_self, d_other = direct[key]["ztf_by_band"], direct[f"{e}/{other}/{band}"]["ztf_by_band"]
        margin = (r["real_max_S"] - r["T"]) / r["T"]
        reasons = []
        a_o = d_other["all"]
        if a_o["S_track"] < 2.0:
            reasons.append(f"other role ({other}) absent along the same track: ZTF all-band S={a_o['S_track']:.2f} "
                           f"(T={a_o['T_controls']:.2f}, n={a_o['n_valid']}); Rx/Tx loci coincide so a real source must appear in both")
        if min(r["control_maxima"]) > 0.8 * r["T"] and r["T"] > 10:
            reasons.append(f"field-wide systematic: all 8 control trajectories at S {min(r['control_maxima']):.1f}-{max(r['control_maxima']):.1f}; "
                           f"real-track margin {100 * margin:.1f} % is not point-source evidence")
        gaia_hits = [h for n in near for h in n["gaia_within_5arcsec"]]
        if gaia_hits and any(float(h["_r"] or 9) < 2.5 for h in gaia_hits):
            reasons.append("Gaia DR3 source within 2.5\" of the track at a test epoch")
        status = "vetoed" if reasons else "marginal"
        report[key] = {"cell": {"z_au": z, "mu": mu, "t0": t0}, "S_real": r["real_max_S"], "T": r["T"],
                       "control_maxima": r["control_maxima"], "direct_ztf_self": d_self,
                       "direct_ztf_other_role": d_other, "catalog_along_track": near,
                       "status": status, "reasons": reasons}
        print(f"{key}: S={r['real_max_S']:.2f}/T={r['T']:.2f} other-role all S={a_o['S_track']:.2f}/T={a_o['T_controls']:.2f} -> {status.upper()}")
        for rs in reasons:
            print("   -", rs)
        sup.append(Candidate(**{**c, "status": status, "veto_reason": "; ".join(reasons) or None,
                                "extra": {**c["extra"], "adjudication": "stage-7 v2 (adjudicate_v2.py)",
                                          "supersedes_status": "retained"}}))
    (OUT / "adjudication.json").write_text(json.dumps(report, indent=2, default=float))
    p = OUT / "records" / "candidate_adjudicated.jsonl"
    if p.exists():
        p.unlink()
    append_records(p, sup)


if __name__ == "__main__":
    main()
