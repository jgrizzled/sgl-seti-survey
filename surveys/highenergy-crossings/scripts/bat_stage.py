"""Swift/BAT second family (hypotheses.md §8/D10, constraint-only): for
every grazing-family window (<= 2.5 Rsun, both channels, Swift era) the
BAT survey pointings whose boresight is within 20 deg of the channel
position, processed with HEASoft batsurvey (chbrandt/heasoft container,
local Swift BAT CALDB) with all 14 channel positions + the Crab as the
input catalogue, so every pointing yields an 8-band rate +/- error at the
positions. Per window: exposure-weighted mean rate over the overlapping
pointings, 3-sigma upper limit in batsurvey units, converted to mCrab by
the Crab rate measured in the same chain (pointings that contain the
Crab), then to 14-195 keV flux (Crab = 2.31e-8 erg/cm2/s) and to power
through the rung cone. No threshold, no trial.

Stages: query (swiftmastr per window), fetch (obs trees), run (batsurvey
per obsid, N workers), collect. Writes results/bat_v1.json / .md."""
from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

import hlib as H
import recon_scan as R

BAT = H.RUN / "bat"
OBS = BAT / "obs"
OUT = BAT / "out"
CATALOG = BAT / "positions.cat"
IMAGE = "chbrandt/heasoft:latest"
HEADAS = "/usr/local/heasoft/x86_64-pc-linux-gnu-libc2.12"
OFFSET_MAX_DEG = 20.0
SWIFT_ERA = ("2004-11-20", "2026-09-07")
CRAB = (83.633, 22.0145)
CRAB_FLUX_14_195 = 2.31e-8       # erg cm^-2 s^-1
MJDREF_SWIFT = 51910.0 + 7.428703703703703e-4
N_WORKERS = 8


def positions_catalog():
    pos = R.positions(R.load_events())
    names, ras, decs = [], [], []
    for (ch, tid), (ra, de, n) in sorted(pos.items()):
        names.append(f"{ch}-{tid}"); ras.append(ra); decs.append(de)
    names.append("crab"); ras.append(CRAB[0]); decs.append(CRAB[1])
    t = Table({"NAME": names, "RA_OBJ": ras, "DEC_OBJ": decs, "CATNUM": list(range(1, len(names) + 1))})
    BAT.mkdir(parents=True, exist_ok=True)
    t.write(CATALOG, format="fits", overwrite=True)
    return {n: (r, d) for n, r, d in zip(names, ras, decs)}


def stage_query(store):
    """Per grazing window (2.5 Rsun), swiftmastr pointings within 20 deg with
    BAT survey exposure overlapping the window."""
    events = R.load_events()
    era = [Time(x).mjd for x in SWIFT_ERA]
    out = []
    for ch, tab in events.items():
        for ev in tab:
            if float(ev["b_min_au"]) >= 2.5 * H.RSUN_AU or not (era[0] <= float(ev["mjd"]) <= era[1]):
                continue
            w = R.windows_for(ev)
            w25 = w["2.5Rsun"]; w12 = w["1.2Rsun"]
            ra, de = R.source_radec(ch, ev)
            adql = ("SELECT obsid, name, ra, dec, start_time, stop_time, bat_expo_sv FROM swiftmastr WHERE "
                    f"start_time <= {w25[1]:.5f} AND stop_time >= {w25[0]:.5f} AND bat_expo_sv > 0 AND "
                    f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.5f}, {de:.5f}, {OFFSET_MAX_DEG}))=1")
            rows, err = R.heasarc(store, adql, f"bat-q {ch} {ev['target_id']} {ev['t_ca_utc'][:10]}")
            out.append({"channel": ch, "target": str(ev["target_id"]), "key": f"{ch}-{ev['target_id']}", "event_id": str(ev["event_id"]),
                        "t_ca_utc": str(ev["t_ca_utc"]), "b_rsun": float(ev["b_min_au"]) / H.RSUN_AU, "ra": ra, "dec": de,
                        "w25": [w25[0], w25[1]], "w12": [w12[0], w12[1]] if w12 else None,
                        "pointings": [{"obsid": r["obsid"], "start": R.ffloat(r["start_time"]), "stop": R.ffloat(r["stop_time"]),
                                       "offset_deg": R.sep_deg(R.ffloat(r["ra"]), R.ffloat(r["dec"]), ra, de), "bat_expo_sv": R.ffloat(r["bat_expo_sv"])}
                                      for r in rows or []], "error": err})
            print(f"[bat-q] {ch} {ev['target_id']:11s} {ev['t_ca_utc'][:10]} -> {len(rows or [])} pointings", flush=True)
    (BAT / "windows_query.json").write_text(json.dumps(out, indent=1))
    return out


def fetch_obs(obsid, mjd):
    d = OBS / obsid
    if (d / "fetched.flag").exists():
        return True
    ym = Time(mjd, format="mjd").datetime.strftime("%Y_%m")
    base = f"https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{ym}/{obsid}/"
    ok = True
    for sub in ("bat/survey", "bat/hk", "auxil"):
        (d / sub).mkdir(parents=True, exist_ok=True)
        try:
            idx = requests.get(base + sub + "/", timeout=120)
            if idx.status_code != 200:
                ok = False; continue
            for f in sorted(set(re.findall(r'href="(sw[^"]+)"', idx.text))):
                p = d / sub / f
                if not p.exists():
                    r = requests.get(base + sub + "/" + f, timeout=600)
                    if r.status_code == 200:
                        p.write_bytes(r.content)
        except Exception as exc:
            print(f"[bat-fetch] {obsid} {sub}: {exc}", file=sys.stderr, flush=True); ok = False
    if ok and any((d / "bat" / "survey").glob("*.dph*")):
        (d / "fetched.flag").write_text(Time.now().isot)
        return True
    return False


def run_batsurvey(obsid):
    outd = OUT / obsid
    if (outd / "stats_point.fits").exists() or (outd / "done.flag").exists():
        return True
    cmd = (f"export HEADAS={HEADAS}; . $HEADAS/headas-init.sh; export HEADASNOQUERY=1 HEADASPROMPT=/dev/null; "
           "export CALDB=/work/bat/caldb/local CALDBCONFIG=/work/bat/caldb/local/caldb.config CALDBALIAS=/work/bat/caldb/local/alias_config.fits; "
           f"export PFILES=\"/tmp/pfiles_{obsid};$HEADAS/syspfiles\"; mkdir -p /tmp/pfiles_{obsid}; cd /work/bat; "
           f"batsurvey obs/{obsid} out/{obsid} detthresh=6 detthresh2=6 incatalog=/work/bat/positions.cat clobber=yes > out/{obsid}.log 2>&1; "
           f"touch out/{obsid}/done.flag")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / obsid).mkdir(exist_ok=True)
    subprocess.run(["docker", "run", "--rm", "--user", "1000:1000", "-e", "HOME=/tmp", "--entrypoint", "/bin/bash", "-v", f"{BAT}:/work/bat:z", IMAGE, "-c", cmd],
                   capture_output=True, timeout=3600)
    return (OUT / obsid / "done.flag").exists()


def collect_obsid(obsid, names):
    """Per pointing: rates at the catalogue positions from point_*_2.cat."""
    rows = []
    for cat in sorted((OUT / obsid).glob("point_*/point_*_2.cat")):
        try:
            with fits.open(cat) as h:
                d = h[1].data; hd = h[1].header
                for r in d:
                    nm = str(r["NAME"]).strip()
                    if nm not in names:
                        continue
                    t0 = float(r["TIME"]); t1 = float(r["TIME_STOP"])
                    rate = np.asarray(r["RATE"], dtype=float); err = np.asarray(r["RATE_ERR"], dtype=float)
                    rows.append({"obsid": obsid, "point": cat.parent.name, "name": nm, "met0": t0, "met1": t1,
                                 "mjd0": MJDREF_SWIFT + t0 / 86400.0, "mjd1": MJDREF_SWIFT + t1 / 86400.0,
                                 "exposure": float(r["EXPOSURE"]) if "EXPOSURE" in d.columns.names else float(hd.get("EXPOSURE", 0)),
                                 "pcodefr": float(r["PCODEFR"]), "rate_tot": float(np.nansum(rate)), "err_tot": float(np.sqrt(np.nansum(err ** 2))),
                                 "rate_bands": rate.tolist(), "err_bands": err.tolist(), "snr": float(r["SNR"]) if "SNR" in d.columns.names else None})
        except Exception as exc:
            print(f"[bat-collect] {cat}: {exc}", file=sys.stderr, flush=True)
    return rows


def main(stages):
    from sglsurvey.snapshots import SnapshotStore
    store = SnapshotStore(BAT / "query")
    names = positions_catalog()
    if "query" in stages or not (BAT / "windows_query.json").exists():
        wins = stage_query(store)
    else:
        wins = json.loads((BAT / "windows_query.json").read_text())
    todo = {}
    for w in wins:
        for p in w["pointings"]:
            if p["offset_deg"] <= OFFSET_MAX_DEG:
                todo[p["obsid"]] = p["start"]
    print(f"[bat] {len(wins)} windows, {len(todo)} distinct pointings", flush=True)
    if "fetch" in stages:
        with ThreadPoolExecutor(6) as ex:
            res = list(ex.map(lambda kv: fetch_obs(*kv), sorted(todo.items())))
        print(f"[bat] fetched {sum(res)}/{len(res)}", flush=True)
    if "run" in stages:
        ok = [o for o in sorted(todo) if (OBS / o / "fetched.flag").exists()]
        with ThreadPoolExecutor(N_WORKERS) as ex:
            res = list(ex.map(run_batsurvey, ok))
        print(f"[bat] batsurvey done {sum(res)}/{len(ok)}", flush=True)
    if "collect" in stages:
        rates = []
        for o in sorted(todo):
            if (OUT / o / "done.flag").exists():
                rates += collect_obsid(o, set(names))
        (BAT / "rates_all.json").write_text(json.dumps(rates, indent=1))
        # Crab calibration
        crab = [r for r in rates if r["name"] == "crab" and r["pcodefr"] >= 0.3 and r["exposure"] > 100]
        crab_rate = float(np.median([r["rate_tot"] for r in crab])) if crab else None
        out = {"utc": Time.now().isot, "crab_rate_units": crab_rate, "n_crab_pointings": len(crab), "crab_flux_14_195": CRAB_FLUX_14_195,
               "offset_max_deg": OFFSET_MAX_DEG, "windows": []}
        for w in wins:
            for rung, span in (("2.5Rsun", w["w25"]), ("1.2Rsun", w["w12"])):
                if span is None:
                    continue
                rr = [r for r in rates if r["name"] == w["key"] and r["mjd1"] > span[0] and r["mjd0"] < span[1] and r["pcodefr"] >= 0.1 and r["err_tot"] > 0]
                rec = {"channel": w["channel"], "target": w["target"], "rung": rung, "t_ca_utc": w["t_ca_utc"], "b_rsun": w["b_rsun"],
                       "n_pointings": len(rr), "exposure_s": float(sum(r["exposure"] for r in rr))}
                if rr:
                    wgt = np.array([1 / r["err_tot"] ** 2 for r in rr]); rt = np.array([r["rate_tot"] for r in rr])
                    mean = float((wgt * rt).sum() / wgt.sum()); err = float(1 / np.sqrt(wgt.sum()))
                    rec.update({"rate": mean, "rate_err": err, "snr": mean / err, "ul3_rate": max(mean, 0.0) + 3 * err,
                                "max_pointing_snr": float(max(r["rate_tot"] / r["err_tot"] for r in rr))})
                    if crab_rate:
                        rec["ul3_mcrab"] = rec["ul3_rate"] / crab_rate * 1e3
                        rec["ul3_flux_erg"] = rec["ul3_rate"] / crab_rate * CRAB_FLUX_14_195
                        b_cm = H.RUNG_R_AU[rung] * H.CM_PER_AU
                        rec["ul3_power_W"] = rec["ul3_flux_erg"] * np.pi * b_cm ** 2 * 1e-7
                out["windows"].append(rec)
        (H.RES / "bat_v1.json").write_text(json.dumps(H.clean(out), indent=1))
        cov = [r for r in out["windows"] if r["n_pointings"] > 0]
        md = ["# Swift/BAT second family v1 (constraint-only) — 2026-09-07", "",
              f"Crab rate in the chain: {crab_rate} (batsurvey units, {len(crab)} pointings); Crab 14–195 keV = {CRAB_FLUX_14_195:.2e} erg/cm²/s.",
              f"Windows with ≥ 1 usable pointing (offset ≤ {OFFSET_MAX_DEG}°, pcode ≥ 0.1): {len(cov)} of {len(out['windows'])}.", "",
              "| Channel | Target | Rung | Windows | covered | median 3σ UL (mCrab) | min UL (mCrab) | median 3σ power (W) | max pointing SNR |", "|---|---|---|---|---|---|---|---|---|"]
        keys = sorted(set((r["channel"], r["target"], r["rung"]) for r in out["windows"]))
        for k in keys:
            rs = [r for r in out["windows"] if (r["channel"], r["target"], r["rung"]) == k]
            cs = [r for r in rs if r.get("ul3_mcrab")]
            if cs:
                md.append(f"| {k[0]} | {k[1]} | {k[2]} | {len(rs)} | {len(cs)} | {np.median([r['ul3_mcrab'] for r in cs]):.0f} | {min(r['ul3_mcrab'] for r in cs):.0f} | "
                          f"{np.median([r['ul3_power_W'] for r in cs]):.2e} | {max(r['max_pointing_snr'] for r in cs):.1f} |")
            else:
                md.append(f"| {k[0]} | {k[1]} | {k[2]} | {len(rs)} | 0 | — | — | — | — |")
        hi = sorted([r for r in cov if r.get("snr", 0) > 3], key=lambda r: -r["snr"])
        md += ["", f"Windows with combined SNR > 3 (constraint-only, recorded): {len(hi)}"] + [f"- {r['channel']} {r['target']} {r['rung']} {r['t_ca_utc'][:10]}: rate {r['rate']:.2e} ± {r['rate_err']:.2e} (SNR {r['snr']:.1f}, {r['n_pointings']} pointings)" for r in hi[:20]]
        (H.RES / "bat_v1.md").write_text("\n".join(md))
        print("\n".join(md[:8]))


if __name__ == "__main__":
    main(set(sys.argv[1:]) or {"query", "fetch", "run", "collect"})
