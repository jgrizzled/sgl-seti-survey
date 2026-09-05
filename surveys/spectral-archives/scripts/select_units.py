"""Fix the v1 unit population and null ensembles from the recon rows
(hypotheses D2/D3, thresholds §1) before any download.

Output: results/units_v1.json
  combos: {target|instrument: {n_pool, null_draw:[...], in_window:[...]}}
  units:  [{unit_id, target, instrument, event_id, t_ca, b_rsun, rungs,
            spectra:[...], cells}]
Each spectrum record carries the download route (url / zip member).
"""
from __future__ import annotations

import json
import random
from collections import defaultdict
from pathlib import Path

from astropy.time import Time

HERE = Path(__file__).resolve().parents[1]
CFG = json.load(open(HERE / "configs" / "threshold_freeze_v1_1.json"))
rows = json.load(open(HERE / "results" / "recon_rows_v0.json"))
scan = json.load(open(HERE / "results" / "recon_scan_v0.json"))
ONSTAR = scan["onstar_arcsec"]
SEED = CFG["seed"]
INSTR = CFG["instruments"]
ESO_FILE = "https://dataportal.eso.org/dataportal_new/file/{dp_id}"
CADC_FILE = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/CFHT/{pid}.fits"
CARM_ZIP = "http://carmenes.cab.inta-csic.es/gto/getDR1DataPublic.action?id={karmn}_VIS.zip"  # karmn percent-encoded below
KARMN = {"teegarden": "J02530+168", "wolf-359": "J10564+070", "ross-128": "J11477+008", "ross-154": "J18498-238"}


def route(r):
    """Return (instrument_key, record) for eligible rows, else None."""
    a = r["archive"]; inst = r["instrument"]
    if a == "ESO-phase3" and inst in ("HARPS", "ESPRESSO", "NIRPS", "XSHOOTER"):
        if inst == "XSHOOTER":                      # per-arm ensembles (hypotheses §2)
            lo = (r.get("band_nm") or [0])[0]
            inst = "XSHOOTER-UVB" if lo < 400 else ("XSHOOTER-VIS" if lo < 900 else "XSHOOTER-NIR")
        return inst, {"url": ESO_FILE.format(dp_id=r["ident"]), "file": f"{r['ident']}.fits"}
    if a == "CADC/CFHT" and inst == "SPIRou":
        pid = r["ident"].split("/")[-1]
        if not pid.endswith("t"):
            return None
        return inst, {"url": CADC_FILE.format(pid=pid), "file": f"{pid}.fits"}
    if a == "CARMENES-DR1" and inst == "CARMENES-VIS":
        return inst, {"zip": CARM_ZIP.format(karmn=KARMN[r["target_id"]].replace("+", "%2B")), "member": r["ident"], "file": r["ident"]}
    return None


def main():
    combos = defaultdict(lambda: {"in_window": [], "pool": []})
    seen = set()
    for r in rows:
        if r["sep_arcsec"] > ONSTAR or r.get("t_mid_mjd") is None:
            continue
        if str(r.get("public")) != "True":
            continue
        rt = route(r)
        if rt is None:
            continue
        inst, rec = rt
        key = f"{r['target_id']}|{inst}"
        dk = (key, rec["file"])
        if dk in seen:
            continue
        seen.add(dk)
        rec.update(target=r["target_id"], instrument=inst, t_mid_mjd=r["t_mid_mjd"],
                   utc=Time(r["t_mid_mjd"], format="mjd").utc.isot[:19], exptime_s=r.get("exptime_s"),
                   in_window=r.get("in_window") or {}, nearest_event=r.get("nearest_event"),
                   nearest_dt_days=r.get("nearest_event_dt_days"), band_nm=r.get("band_nm"))
        if rec["in_window"]:
            combos[key]["in_window"].append(rec)
        else:
            combos[key]["pool"].append(rec)

    rng = random.Random(SEED)
    units = []
    out_combos = {}
    evinfo = {}
    for tid, td in scan["targets"].items():
        for e in td["events"]:
            evinfo[e["event_id"]] = e
    for key in sorted(combos):
        c = combos[key]
        tid, inst = key.split("|")
        n_pool = len(c["pool"])
        eligible = n_pool >= CFG["min_null_pool"] and len(c["in_window"]) > 0
        # stratified-by-year seeded draw
        draw = []
        if eligible:
            byyear = defaultdict(list)
            for p in sorted(c["pool"], key=lambda p: p["t_mid_mjd"]):
                byyear[p["utc"][:4]].append(p)
            years = sorted(byyear)
            need = min(CFG["n_null"], n_pool)
            # proportional allocation with at least 1 per year where possible
            alloc = {y: 1 for y in years} if need >= len(years) else {y: 0 for y in years}
            rem = need - sum(alloc.values())
            weights = {y: len(byyear[y]) - alloc[y] for y in years}
            while rem > 0:
                tot = sum(max(w, 0) for w in weights.values())
                if tot == 0:
                    break
                for y in years:
                    if rem == 0:
                        break
                    share = round(need * len(byyear[y]) / n_pool)
                    if alloc[y] < min(share, len(byyear[y])):
                        alloc[y] += 1; weights[y] -= 1; rem -= 1
                else:
                    # fill leftovers from years with capacity
                    for y in years:
                        if rem == 0:
                            break
                        if alloc[y] < len(byyear[y]):
                            alloc[y] += 1; rem -= 1
            for y in years:
                pool_y = list(byyear[y])
                rng.shuffle(pool_y)
                draw += pool_y[:alloc[y]]
            draw = sorted(draw, key=lambda p: p["t_mid_mjd"])
        out_combos[key] = {"target": tid, "instrument": inst, "n_pool": n_pool, "n_in_window": len(c["in_window"]),
                           "eligible": eligible, "reason": None if eligible else ("no in-window spectra" if not c["in_window"] else f"null pool {n_pool} < {CFG['min_null_pool']}"),
                           "null_draw": draw if eligible else [], "in_window": c["in_window"]}
        if not eligible:
            continue
        byev = defaultdict(list)
        for s in c["in_window"]:
            eid = s["in_window"].get("0.1AU") or next(iter(s["in_window"].values()))
            byev[eid].append(s)
        band = INSTR[inst]["band_nm"] if inst in INSTR else {"XSHOOTER-UVB": [300, 560], "XSHOOTER-VIS": [530, 1020], "XSHOOTER-NIR": [990, 2480]}[inst]
        cells = [cn for cn, iv in CFG["cells_nm"].items() if iv is None or (iv[0] >= band[0] and iv[1] <= band[1])]
        for eid, specs in sorted(byev.items(), key=lambda kv: evinfo[kv[0]]["t_ca_mjd"]):
            e = evinfo[eid]
            rungs = sorted(set(r for s in specs for r in s["in_window"]))
            units.append({"unit_id": f"{tid}|{inst}|{e['t_ca_utc'][:10]}", "target": tid, "instrument": inst,
                          "event_id": eid, "t_ca_utc": e["t_ca_utc"], "t_ca_mjd": e["t_ca_mjd"], "b_rsun": e["b_rsun"],
                          "rungs": rungs, "family": "dev" if tid in CFG["split"]["dev"] else "confirmatory",
                          "cells": cells, "n_spectra": len(specs),
                          "spectra": sorted(specs, key=lambda s: s["t_mid_mjd"])})
    out = {"frozen": CFG["frozen_utc"], "seed": SEED, "combos": out_combos, "units": units}
    (HERE / "results" / "units_v1.json").write_text(json.dumps(out, indent=1))
    print(f"{len(units)} units")
    for u in units:
        print(f"  {u['family']:12s} {u['unit_id']:36s} n={u['n_spectra']:2d} rungs={'+'.join(u['rungs']):22s} cells={u['cells']}")
    for k, c in out_combos.items():
        print(f"  combo {k:26s} pool {c['n_pool']:4d} in-window {c['n_in_window']:3d} eligible={c['eligible']} draw={len(c['null_draw'])} {c['reason'] or ''}")


if __name__ == "__main__":
    main()
