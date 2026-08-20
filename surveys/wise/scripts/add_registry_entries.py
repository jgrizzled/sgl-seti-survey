"""Generate linear-astrometry registry entries for a batch of universal-
list systems and append them to registries/pilot_wise_2026.yaml.

Earlier batches were curated by hand; this script makes the routine
part reproducible. Per component with a Gaia DR3 id: pull the 5p
solution + RV from VizieR I/355/gaiadr3 (epoch 2016.0 TCB). Gaia-
saturated stars (no DR3 id, e.g. Fomalhaut) fall back to Hipparcos-2
(I/311/hip2, epoch 1991.25) with RV from CNS5/SIMBAD. RVs missing from
both sources get a 0 ± 50 km/s allocation (locus impact < 0.1",
van-Maanen precedent). Uncertainties are inflated x3 when RUWE > 1.4
(tau Cet precedent). Every entry carries a note naming the basket(s)
that motivated it and any multiplicity approximation.

Usage:
    uv run python surveys/wise/scripts/add_registry_entries.py \
        --batch batch-5 --tag universal-v2 [--dry-run] SYSTEM [SYSTEM ...]
SYSTEM names are `system` values from targets/universal_v2.json; append
`:A` to restrict to one component (e.g. "GJ 229 AB:A").
Prints the CORRIDOR_OF lines to paste into wise_corridors.py.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import time
from pathlib import Path

import requests
import yaml

REPO = Path(__file__).resolve().parents[3]
REGISTRY = REPO / "registries" / "pilot_wise_2026.yaml"
ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
TODAY = "2026-08-20"

RV_OVERRIDES = {  # eid: (rv, err, source) for stars lacking catalog RV
    "fomalhaut": (6.5, 0.5, "Gontcharov 2006 (PCRV) via SIMBAD"),
}
SLUG_OVERRIDES = {"82 Eri": "82-eri", "61 Vir": "61-vir",
                  "sigma Dra": "sigma-dra", "HD 219134": "hd-219134",
                  "Fomalhaut": "fomalhaut", "Wolf 437": "wolf-437",
                  "Wolf 1069": "wolf-1069", "LHS 1723": "lhs-1723"}


def slug(name, comp=""):
    base = SLUG_OVERRIDES.get(name)
    if base is None:
        n = re.sub(r"\(.*?\)", "", name).strip()
        n = re.sub(r"\s+(AB|ABC|ABCD)$", "", n)
        base = re.sub(r"[^a-z0-9]+", "-", n.lower()).strip("-")
    return base + (f"-{comp.lower()}" if comp else "")


def asu(source, out, retries=5, **constraints):
    p = {"-source": source, "-out.max": "50", "-out": ",".join(out)}
    p.update(constraints)
    for k in range(retries):
        r = requests.get(ASU, params=p, timeout=120)
        if r.status_code == 200 and "Postgres connect error" not in r.text:
            lines = [l for l in r.text.splitlines()
                     if l and not l.startswith("#")]
            if not lines:
                return []
            hdr = lines[0].split("\t")
            return [dict(zip(hdr, [v.strip() for v in l.split("\t")]))
                    for l in lines[3:] if len(l.split("\t")) == len(hdr)]
        time.sleep(2 * (k + 1))
    raise RuntimeError(f"ASU failed for {source}")


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def gaia_entry(gid):
    rows = asu("I/355/gaiadr3",
               ["Source", "RA_ICRS", "DE_ICRS", "e_RA_ICRS", "e_DE_ICRS",
                "Plx", "e_Plx", "pmRA", "e_pmRA", "pmDE", "e_pmDE",
                "RV", "e_RV", "RUWE"], Source="=" + gid)
    if not rows:
        return None
    r = rows[0]
    return {"ra": f(r["RA_ICRS"]), "dec": f(r["DE_ICRS"]),
            "e_ra": f(r["e_RA_ICRS"]), "e_dec": f(r["e_DE_ICRS"]),
            "plx": f(r["Plx"]), "e_plx": f(r["e_Plx"]),
            "pmra": f(r["pmRA"]), "e_pmra": f(r["e_pmRA"]),
            "pmde": f(r["pmDE"]), "e_pmde": f(r["e_pmDE"]),
            "rv": f(r["RV"]), "e_rv": f(r["e_RV"]), "ruwe": f(r["RUWE"]),
            "epoch": 2016.0, "scale": "tcb",
            "src": f"Gaia DR3 {gid} via VizieR I/355/gaiadr3 (retrieved {TODAY})",
            "catalog": "Gaia DR3", "id": gid}


def hip2_entry(hip):
    rows = asu("I/311/hip2",
               ["HIP", "RArad", "DErad", "e_RArad", "e_DErad", "Plx",
                "e_Plx", "pmRA", "e_pmRA", "pmDE", "e_pmDE"], HIP="=" + hip)
    if not rows:
        return None
    r = rows[0]
    ra, de = f(r["RArad"]), f(r["DErad"])
    # VizieR's hip2 RArad/DErad arrive scaled by 180/pi relative to
    # degrees (observed 2026-08-20: Fomalhaut 19733.34 -> 344.41 deg).
    import math
    while abs(ra) >= 360.0:
        ra, de = ra / (180 / math.pi), de / (180 / math.pi)
    return {"ra": ra, "dec": de,
            "e_ra": f(r["e_RArad"]), "e_dec": f(r["e_DErad"]),
            "plx": f(r["Plx"]), "e_plx": f(r["e_Plx"]),
            "pmra": f(r["pmRA"]), "e_pmra": f(r["e_pmRA"]),
            "pmde": f(r["pmDE"]), "e_pmde": f(r["e_pmDE"]),
            "rv": None, "e_rv": None, "ruwe": None,
            "epoch": 1991.25, "scale": "tt",
            "src": f"Hipparcos-2 (van Leeuwen 2007) HIP {hip} via VizieR I/311/hip2 (retrieved {TODAY})",
            "catalog": "Hipparcos", "id": hip}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("systems", nargs="+")
    ap.add_argument("--batch", required=True)
    ap.add_argument("--tag", default="universal-v2")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    uni = json.load(open(REPO / "targets" / "universal_v2.json"))
    by_name = {s["system"]: s for s in uni["systems"]}
    census = {r["CNS5"]: r for r in csv.DictReader(
        open(REPO / "targets" / uni["config"]["census"]["snapshot"]))}
    # Map Gaia id -> CNS5 row for RV fallback and HIP lookup
    by_gaia = {r["GaiaDR3"]: r for r in census.values() if r["GaiaDR3"]}
    by_sysname = {}
    reg = yaml.safe_load(REGISTRY.read_text())
    existing = set(reg["targets"])

    entries, corridors = {}, []
    for spec in a.systems:
        name, _, only = spec.partition(":")
        s = by_name[name]
        label = ", ".join(b.replace("_", " ").title() for b in s["baskets"])
        for c in s["components"]:
            if only and c["comp"] != only:
                continue
            comp = c["comp"] if s["n_comp"] > 1 or len(c["comp"]) > 1 else ""
            if len(comp) > 1:     # merged CNS5 row (e.g. "ABC"): primary
                comp = comp[0]
            eid = slug(name, comp)
            if eid in existing:
                print(f"  skip {eid}: already in registry")
                continue
            ast = gaia_entry(c["gaia_dr3"]) if c["gaia_dr3"] else None
            cns = by_gaia.get(c["gaia_dr3"]) if c["gaia_dr3"] else None
            if ast is None:
                # Gaia-saturated: find the CNS5 row by name to get HIP
                hip = next((r["HIP"] for r in census.values()
                            if r["GJ"].strip() == s.get("gj", "") and r["HIP"]), None)
                if not hip:
                    print(f"  !! {eid}: no Gaia id and no HIP — skipped")
                    continue
                ast = hip2_entry(hip)
                cns = next(r for r in census.values() if r["HIP"] == hip)
            if ast is None:
                print(f"  !! {eid}: no astrometry — skipped")
                continue
            notes = [f"Universal v2 baskets: {label}."]
            rv_src = ast["src"]
            if ast["rv"] is None:
                if eid in RV_OVERRIDES:
                    ast["rv"], ast["e_rv"], rv_src = RV_OVERRIDES[eid]
                elif cns and cns["RV"]:
                    ast["rv"], ast["e_rv"] = f(cns["RV"]), f(cns["e_RV"]) or 1.0
                    rv_src = f"CNS5 compilation RV (retrieved {TODAY})"
                else:
                    ast["rv"], ast["e_rv"] = 0.0, 50.0
                    rv_src = "no catalog RV; 0 +/- 50 km/s allocation (locus impact < 0.1\")"
                    notes.append("RV unavailable: +/-50 km/s allocation.")
            infl = 1.0
            if ast["ruwe"] is not None and ast["ruwe"] > 1.4:
                infl = 3.0
                notes.append(f"Gaia RUWE {ast['ruwe']:.2f}: uncertainties inflated x3 (tau Cet precedent).")
            elif ast["ruwe"] is not None:
                notes.append(f"RUWE {ast['ruwe']:.2f}.")
            if s["multiplicity"] != "single":
                notes.append("Multiple system: component endpoint on its own Gaia "
                             "solution; orbital curvature over the WISE baseline "
                             "negligible at the declared separation "
                             f"({s['min_sep_au']} AU) — Struve 2398 precedent.")
            if ast["catalog"] == "Hipparcos":
                notes.append("Gaia-saturated: Hipparcos-2 solution adopted "
                             "(Sirius/Procyon precedent; single star, no "
                             "barycenter correction needed).")
            disp = name + (f" {comp}" if comp else "")
            if s["engineering"].get("sptype"):
                disp += f" ({s['engineering']['sptype']})"

            def prov(src, unc, unit):
                return {"source": src, "uncertainty": round(unc * infl, 4), "unit": unit}
            entries[eid] = {
                "display_name": disp,
                "endpoint_kind": "component" if comp else "star",
                "priority": 0.5,
                "tags": [a.batch, a.tag],
                "notes": " ".join(notes),
                "identifiers": [{"catalog": ast["catalog"], "id": ast["id"],
                                 "version": ast["src"]}],
                "quality": ("Gaia DR3 five-parameter solution" if ast["catalog"] == "Gaia DR3"
                            else "Hipparcos-2 solution") + ("; see notes" if len(notes) > 1 else ""),
                "model_rationale": ("Linear astrometry family; orbital/companion effects "
                                    "either negligible at WISE-baseline precision or "
                                    "absorbed into inflated uncertainties and the "
                                    "residual-motion search cell (see notes)."),
                "state": {
                    "provider": "linear_astrometry_v1",
                    "astrometry": {
                        "frame": "icrs", "ra_deg": ast["ra"], "dec_deg": ast["dec"],
                        "parallax_mas": ast["plx"],
                        "pm_ra_cosdec_mas_per_yr": ast["pmra"],
                        "pm_dec_mas_per_yr": ast["pmde"],
                        "radial_velocity_km_s": ast["rv"],
                        "reference_epoch_jyear": ast["epoch"],
                        "reference_epoch_scale": ast["scale"],
                        "source": ast["src"] + "; RV per provenance"},
                    "provenance": {
                        "ra_deg": prov(ast["src"], ast["e_ra"], "mas"),
                        "dec_deg": prov(ast["src"], ast["e_dec"], "mas"),
                        "parallax_mas": prov(ast["src"], ast["e_plx"], "mas"),
                        "pm_ra_cosdec_mas_per_yr": prov(ast["src"], ast["e_pmra"], "mas/yr"),
                        "pm_dec_mas_per_yr": prov(ast["src"], ast["e_pmde"], "mas/yr"),
                        "radial_velocity_km_s": prov(rv_src, ast["e_rv"], "km/s"),
                    }}}
            corridors.append((eid, slug(name).replace("-", "")))
            print(f"  {eid:16s} {ast['catalog']:10s} plx={ast['plx']:.2f} "
                  f"ruwe={ast['ruwe']} rv={ast['rv']}")

    print("\n# paste into wise_corridors.py CORRIDOR_OF:")
    for eid, cor in corridors:
        print(f'    "{eid}": "{cor}",')
    if a.dry_run:
        return
    text = REGISTRY.read_text().rstrip("\n") + "\n\n"
    text += f"  # ---- {a.batch} ({TODAY}): universal v2 picky-network expansion ----\n"
    body = yaml.safe_dump({"targets": entries}, sort_keys=False,
                          allow_unicode=True, width=100)
    body = body.split("\n", 1)[1]  # drop "targets:" line; keep 2-space indent
    REGISTRY.write_text(text + body)
    print(f"\nappended {len(entries)} entries to {REGISTRY.name}")


if __name__ == "__main__":
    main()
