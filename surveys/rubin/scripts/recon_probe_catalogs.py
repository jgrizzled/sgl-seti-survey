"""Rubin DP2 recon: catalog-side probes (depth, density, spatial ADQL).

- per-band VisitDetector depth/seeing aggregates (whole DP2)
- DiaSource / Object / ForcedSource probes in a 6 arcmin cone at the
  ross-128 antipode (the best-covered corridor + the 1.86 Rsun
  crossing event position)
- verifies CONTAINS(POINT,CIRCLE) spatial ADQL on Qserv

Output: surveys/rubin/results/recon_catalog_probes_v0.json
"""
import json
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from recon_fetch_visits import get_service  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[3]
OUT = REPO / "surveys" / "rubin" / "results" / "recon_catalog_probes_v0.json"

RA, DEC, R = 356.943, -0.788, 0.1  # ross-128 antipode, 6 arcmin
CONE = f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {RA}, {DEC}, {R})) = 1"


def run(svc, label, query, out, mode="sync"):
    t0 = time.time()
    try:
        if mode == "async":
            job = svc.submit_job(query)
            job.run()
            job.wait(phases=["COMPLETED", "ERROR"])
            assert job.phase == "COMPLETED", job.phase
            tab = job.fetch_result().to_table()
            job.delete()
        else:
            tab = svc.search(query).to_table()
        rows = [dict(zip(tab.colnames, (None if hasattr(x, "mask") and x.mask
                                        else (x.item() if hasattr(x, "item") else x)
                                        for x in r))) for r in tab]
        out[label] = {"query": query, "mode": mode,
                      "elapsed_s": round(time.time() - t0, 1), "rows": rows}
        print(f"{label}: {len(rows)} rows in {out[label]['elapsed_s']}s")
    except Exception as exc:  # keep probing on individual failures
        out[label] = {"query": query, "mode": mode, "error": str(exc)[:500]}
        print(f"{label}: ERROR {str(exc)[:200]}")


def main():
    svc = get_service()
    out = {"cone": {"ra": RA, "dec": DEC, "radius_deg": R,
                    "note": "ross-128 antipode"}}

    run(svc, "visitdetector_depth",
        "SELECT band, COUNT(*) AS n, AVG(magLim) AS maglim_avg, "
        "MIN(magLim) AS maglim_min, MAX(magLim) AS maglim_max, "
        "AVG(seeing) AS seeing_avg, AVG(zeroPoint) AS zp_avg "
        "FROM dp2.VisitDetector GROUP BY band", out, mode="async")

    run(svc, "diasource_cone",
        f"SELECT band, COUNT(*) AS n, AVG(reliability) AS rel_avg, "
        f"MIN(midpointMjdTai) AS mjd_min, MAX(midpointMjdTai) AS mjd_max "
        f"FROM dp2.DiaSource WHERE {CONE} GROUP BY band", out)

    run(svc, "diaobject_cone",
        f"SELECT COUNT(*) AS n FROM dp2.DiaObject WHERE {CONE}", out)

    run(svc, "object_cone",
        "SELECT COUNT(*) AS n FROM dp2.Object WHERE "
        "CONTAINS(POINT('ICRS', coord_ra, coord_dec), "
        f"CIRCLE('ICRS', {RA}, {DEC}, {R})) = 1", out)

    run(svc, "forcedsource_sample",
        "SELECT TOP 5 o.objectId, o.coord_ra, o.coord_dec "
        "FROM dp2.Object AS o WHERE CONTAINS(POINT('ICRS', o.coord_ra, "
        f"o.coord_dec), CIRCLE('ICRS', {RA}, {DEC}, 0.01)) = 1", out)

    # forced photometry row count for one object (linkage check)
    oid = None
    rows = out.get("forcedsource_sample", {}).get("rows") or []
    if rows:
        oid = rows[0]["objectId"]
        run(svc, "forcedsource_lightcurve",
            "SELECT band, COUNT(*) AS n, COUNT(psfDiffFlux) AS n_diff "
            f"FROM dp2.ForcedSource WHERE objectId = {oid} GROUP BY band",
            out)

    run(svc, "ssobject_total",
        "SELECT COUNT(*) AS n FROM dp2.SSObject", out)

    OUT.write_text(json.dumps(out, indent=1, default=str))
    print("wrote", OUT.relative_to(REPO))


if __name__ == "__main__":
    main()
