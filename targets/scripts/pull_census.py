"""Snapshot the CNS5 census and the engineering-column joins.

Produces two files in targets/census/ (dated, never overwritten in
place — bump the date and point config.yaml at the new snapshot):

  cns5_<horizon>pc_<date>.csv        CNS5 rows with plx >= 1000/horizon
  engineering_<date>.csv              per-CNS5-row joins used by the
                                      desirability classifier and the
                                      engineering basket:
      TIC (IV/39/tic82)        Teff, Rad, Mass, Lum   (uniform; Mann+2019
                               relations for cool dwarfs, Torres for hot)
      Gaia DR3 main (I/355)    RUWE, IPDfmp
      Gaia DR3 paramp (I/355)  FLAME mass/radius/age/evolstage, ESP-ELS
                               Halpha EW + active-M-dwarf probability,
                               ESP-CS Ca IRT activity index
      Kervella+2022 (J/A+A/657/A7)  Hip2-EDR3 proper-motion anomaly S/N,
                               binary flag, SIMBAD SpType
      HGCA (J/ApJS/254/42)     chi2 of the acceleration fit
      Boro Saikia+2018 (J/A+A/616/A108)  log R'HK (FGK activity),
                               matched by position (15")
      Jeffers+2018 (J/A+A/614/A76)  CARMENES Halpha pEW + vsini (M
                               dwarfs), matched by position (15")
      2RXS (J/A+A/588/A103)    ROSAT count rate, cone 40" at the
                               PM-propagated 1990.5 position

Sources: VizieR ASU (TAP was throttled on 2026-08-20). Query lists are
chunked; each chunk is retried because CDS intermittently drops the
Postgres connection.

Usage: uv run python targets/scripts/pull_census.py [--horizon-pc 10]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
HERE = Path(__file__).resolve().parent
CENSUS = HERE.parent / "census"

CNS5_COLS = ["CNS5", "GJ", "Comp", "NComp", "P?", "GJp", "GaiaDR3", "HIP",
             "RAJ2000", "DEJ2000", "Epoch", "plx", "e_plx", "pmRA", "pmDE", "RV",
             "e_RV", "Gmag", "RPmag", "Jmag", "Hmag", "Ksmag", "W1mag",
             "SimbadName"]


def asu_query(source, out, retries=6, **constraints):
    params = {"-source": source, "-out.max": "unlimited",
              "-out": ",".join(out)}
    params.update(constraints)
    for attempt in range(retries):
        try:
            r = requests.get(ASU, params=params, timeout=180)
            text = r.text
        except requests.RequestException as e:
            text = f"#ERR {e}"
        if r.status_code == 200 and "Postgres connect error" not in text \
                and "#INFO\tError" not in text:
            break
        time.sleep(3 * (attempt + 1))
    else:
        raise RuntimeError(f"ASU failed for {source}: {text[:300]}")
    lines = [l for l in text.splitlines() if l and not l.startswith("#")]
    if not lines:
        return []
    header = lines[0].split("\t")
    # lines[1] units, lines[2] dashes
    rows = []
    for l in lines[3:]:
        vals = [v.strip() for v in l.split("\t")]
        if len(vals) != len(header):
            continue
        rows.append(dict(zip(header, vals)))
    return rows


def chunked(ids, n=40):
    ids = [i for i in ids if i]
    for k in range(0, len(ids), n):
        yield ids[k:k + n]


def list_query(source, key, ids, out):
    result = {}
    for chunk in chunked(ids):
        for r in asu_query(source, out, **{key: "=" + ",".join(chunk)}):
            result.setdefault(r[key], r)
        print(f"  {source} {key}: {len(result)}/{len(ids)}", end="\r")
    print()
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--horizon-pc", type=float, default=10.0)
    ap.add_argument("--date", default=dt.date.today().isoformat())
    a = ap.parse_args()
    plx_min = 1000.0 / a.horizon_pc

    print(f"CNS5 plx >= {plx_min:.1f} mas")
    rows = asu_query("J/A+A/670/A19/cns5", CNS5_COLS, plx=f">={plx_min:.2f}")
    rows.sort(key=lambda r: -float(r["plx"]))
    if len(rows) < 100:
        raise RuntimeError(f"CNS5 returned only {len(rows)} rows — "
                           "transient CDS failure; nothing written")
    print(f"  {len(rows)} rows")

    gaia = [r["GaiaDR3"] for r in rows]
    hip = [r["HIP"] for r in rows]

    print("engineering joins")
    tic = list_query("IV/39/tic82", "GAIA", gaia,
                     ["GAIA", "TIC", "Teff", "Rad", "Mass", "Lum", "Lclass"])
    gmain = list_query("I/355/gaiadr3", "Source", gaia,
                       ["Source", "RUWE", "IPDfmp"])
    gpar = list_query("I/355/paramp", "Source", gaia,
                      ["Source", "Mass-Flame", "Rad-Flame", "Age-Flame",
                       "Evol", "EWHa", "f_EWHa", "PdMactive", "CA-CS",
                       "Teff", "logg"])
    kerv = list_query("J/A+A/657/A7/tablea1", "HIP", hip,
                      ["HIP", "Name", "SpType", "snrPMaH2EG3b",
                       "BinH2EG3b", "RUWE", "DMS", "W"])
    hgca = list_query("J/ApJS/254/42/catalog", "HIP", hip,
                      ["HIP", "chi2"])

    print("activity joins")
    bs = asu_query("J/A+A/616/A108/catalog",
                   ["Name", "logRpHK", "_RA", "_DE"], Plx=">=90")
    jf = asu_query("J/A+A/614/A76/tablea2",
                   ["Karmn", "pEWHa", "vsini", "_RA", "_DE"])
    print(f"  Boro Saikia rows {len(bs)}, Jeffers rows {len(jf)}")

    def nearest(cat, ra, dec, rs_arcsec):
        best, bd = None, rs_arcsec
        cosd = math.cos(math.radians(dec))
        for r in cat:
            try:
                dra = (float(r["_RA"]) - ra) * cosd
                dde = float(r["_DE"]) - dec
            except ValueError:
                continue
            d = math.hypot(dra, dde) * 3600.0
            if d < bd:
                best, bd = r, d
        return best

    def one_activity(r):
        try:
            ra, dec = float(r["RAJ2000"]), float(r["DEJ2000"])
        except ValueError:
            return None
        pmra = float(r["pmRA"] or 0) / 3.6e6
        pmde = float(r["pmDE"] or 0) / 3.6e6
        cosd = math.cos(math.radians(dec))
        # CNS5 positions are at a per-row epoch (1991-2018, `Epoch`);
        # propagate to J2000 for the SIMBAD-positioned catalogs and to
        # the RASS epoch (~1990.5) for ROSAT.
        ep = float(r["Epoch"] or 2000.0)

        def at(epoch):
            return (ra + pmra * (epoch - ep) / cosd,
                    dec + pmde * (epoch - ep))
        ra00, de00 = at(2000.0)
        b = nearest(bs, ra00, de00, 20.0)
        j = nearest(jf, ra00, de00, 20.0)
        ra90, de90 = at(1990.5)
        rx = asu_query("J/A+A/588/A103/cat2rxs",
                       ["_r", "2RXS", "CRate", "e_CRate"],
                       **{"-c": f"{ra90:.5f} {de90:+.5f}", "-c.rs": "60"})
        rx.sort(key=lambda x: float(x["_r"]))
        return {
            "logRpHK": b["logRpHK"] if b else "",
            "carm_pEWHa": j["pEWHa"] if j else "",
            "carm_vsini": j["vsini"] if j else "",
            "rosat_crate": rx[0]["CRate"] if rx else "",
            "rosat_sep": rx[0]["_r"] if rx else "",
        }

    act = {}
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(one_activity, r): r["CNS5"] for r in rows}
        for i, f in enumerate(as_completed(futs)):
            res = f.result()
            if res:
                act[futs[f]] = res
            if i % 20 == 0:
                print(f"  activity {i + 1}/{len(rows)}", flush=True)

    out = CENSUS / f"cns5_{a.horizon_pc:g}pc_{a.date}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CNS5_COLS)
        w.writeheader()
        w.writerows(rows)
    print(f"  {len(rows)} rows -> {out.name}")

    eng_cols = ["CNS5", "GaiaDR3", "HIP",
                "tic_id", "tic_teff", "tic_rad", "tic_mass", "tic_lum",
                "tic_lclass",
                "ruwe", "ipd_frac_multi_peak",
                "flame_mass", "flame_rad", "flame_age", "flame_evol",
                "ew_halpha", "ew_halpha_flag", "p_active_mdwarf",
                "ca_irt_activity", "gspphot_teff", "gspphot_logg",
                "kerv_name", "sptype", "pma_snr", "pma_binary",
                "hip_dms", "in_wds", "hgca_chi2",
                "logRpHK", "carm_pEWHa", "carm_vsini", "rosat_crate",
                "rosat_sep"]
    out2 = CENSUS / f"engineering_{a.date}.csv"
    with open(out2, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=eng_cols)
        w.writeheader()
        for r in rows:
            t = tic.get(r["GaiaDR3"], {})
            gm = gmain.get(r["GaiaDR3"], {})
            gp = gpar.get(r["GaiaDR3"], {})
            k = kerv.get(r["HIP"], {})
            h = hgca.get(r["HIP"], {})
            w.writerow({
                "CNS5": r["CNS5"], "GaiaDR3": r["GaiaDR3"], "HIP": r["HIP"],
                "tic_id": t.get("TIC", ""), "tic_teff": t.get("Teff", ""),
                "tic_rad": t.get("Rad", ""), "tic_mass": t.get("Mass", ""),
                "tic_lum": t.get("Lum", ""), "tic_lclass": t.get("Lclass", ""),
                "ruwe": gm.get("RUWE", ""),
                "ipd_frac_multi_peak": gm.get("IPDfmp", ""),
                "flame_mass": gp.get("Mass-Flame", ""),
                "flame_rad": gp.get("Rad-Flame", ""),
                "flame_age": gp.get("Age-Flame", ""),
                "flame_evol": gp.get("Evol", ""),
                "ew_halpha": gp.get("EWHa", ""),
                "ew_halpha_flag": gp.get("f_EWHa", ""),
                "p_active_mdwarf": gp.get("PdMactive", ""),
                "ca_irt_activity": gp.get("CA-CS", ""),
                "gspphot_teff": gp.get("Teff", ""),
                "gspphot_logg": gp.get("logg", ""),
                "kerv_name": k.get("Name", ""), "sptype": k.get("SpType", ""),
                "pma_snr": k.get("snrPMaH2EG3b", ""),
                "pma_binary": k.get("BinH2EG3b", ""),
                "hip_dms": k.get("DMS", ""), "in_wds": k.get("W", ""),
                "hgca_chi2": h.get("chi2", ""),
                **act.get(r["CNS5"], {}),
            })
    print(f"  -> {out2.name}")


if __name__ == "__main__":
    sys.exit(main())
