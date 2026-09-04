"""Gattini-IR DR1 reachability recon scan (plan §5.8 item 6; 2026-09-04).

Geometry-only era scope of the universal crossing list against the
PGIR DR1 era, then a per-target-channel *coverage proxy* pull from
NOIRLab Data Lab (``pgir_dr1`` tables, anonymous sync queries,
snapshotted verbatim under runs/gattini-crossings/recon).

Why a proxy: DR1 is a light-curve catalog of 2MASS point sources, not
an image archive. The epoch list of any 2MASS source in a 1.24 deg
sub-quadrant *is* the visit list of that sub-quadrant (non-detections
are kept as rows), so the union of epochs of the nearest few 2MASS
sources to a channel position is the coverage at that position.

Two DR1 defects are quantified per event:
- ``obsjd`` is float32 at the source (paper Table: obsjd float32),
  quantised to 0.25 d. Each epoch is treated as the interval
  [q - 0.125, q + 0.125] d; in-window counts are reported as strict
  (interval inside the flat-chord window), nominal (bin centre
  inside) and loose (interval overlaps).
- Photometry is forced at 2MASS-epoch (1998-2001) positions; the
  proper-motion offset of each target star at the PGIR era is
  measured as the separation between the events-table star position
  and the nearest 2MASS source.

No signal statistic is formed; nothing here is a search.
"""
from __future__ import annotations

import csv
import io
import json
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

RUN = REPO / "runs" / "gattini-crossings" / "recon"
OUT = Path(__file__).resolve().parents[1] / "results" / "recon_scan_v0.json"
EVENTS = REPO / "crossings" / "universal_v1" / "events.ecsv"

# measured 2026-09-04: SELECT MIN(obsjd), MAX(obsjd) FROM pgir_dr1.exposures
ERA_JD = (2458407.0, 2459868.0)
ERA_MJD = (ERA_JD[0] - 2400000.5, ERA_JD[1] - 2400000.5)
DEC_MIN = -28.5           # DR1 footprint floor (Data Lab dataset page)
QUANT_D = 0.25            # float32 obsjd resolution at JD 2.458e6
KM_PER_AU = 1.495978707e8
RUNGS = [("B", "1.2Rsun", 1.2 * 0.00465047), ("B", "2.5Rsun", 2.5 * 0.00465047),
         ("B", "0.1AU", 0.1), ("A", "0.1AU", 0.1)]
CONE_DEG = 0.1            # 6 arcmin: sources cone for the coverage proxy
N_PROXY = 3               # epochs unioned over the N nearest sources


def q(store, sql, timeout=600):
    from dl import queryClient as qc
    request_utc = datetime.now(timezone.utc).isoformat()
    for attempt in range(3):
        try:
            text = qc.query(sql=sql, fmt="csv", timeout=timeout)
            break
        except Exception as exc:          # transient service errors
            if attempt == 2:
                raise
            print(f"  retry after: {exc}", file=sys.stderr)
            time.sleep(5)
    rows = list(csv.DictReader(io.StringIO(text)))
    store.store(service_url="datalab:queryClient/query", query=sql,
                request_utc=request_utc, response_bytes=text.encode(),
                row_count=len(rows), http_status=200)
    return rows


def load_events():
    t = Table.read(EVENTS)
    t["t_ca_mjd"] = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    mjd = np.asarray(t["t_ca_mjd"])
    era = (mjd >= ERA_MJD[0]) & (mjd <= ERA_MJD[1])
    side = np.asarray(t["axis_distance_au"])
    is_a = np.asarray(t["link_direction"]) == "inbound"
    not_sunward = np.where(is_a, side > 0, side < 0)
    return t[era & not_sunward], t


def window(ev, r_au):
    b = float(ev["b_min_au"])
    half_d = np.sqrt(max(r_au * r_au - b * b, 0.0)) * KM_PER_AU / float(ev["v_perp_km_s"]) / 86400.0
    return float(ev["t_ca_mjd"]) - half_d, float(ev["t_ca_mjd"]) + half_d


def count_modes(epochs_mjd, lo, hi):
    e = np.asarray(epochs_mjd, dtype=float)
    if e.size == 0:
        return {"strict": 0, "nominal": 0, "loose": 0}
    h = QUANT_D / 2
    return {"strict": int(((e - h >= lo) & (e + h <= hi)).sum()),
            "nominal": int(((e >= lo) & (e <= hi)).sum()),
            "loose": int(((e + h >= lo) & (e - h <= hi)).sum())}


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    store = SnapshotStore(RUN)
    ev, t_all = load_events()
    is_a = np.asarray(ev["link_direction"]) == "inbound"
    b_au = np.asarray(ev["b_min_au"])
    dec_pos = np.where(is_a, ev["star_icrs_dec_deg"], ev["relay_icrs_dec_deg"])
    ra_pos = np.where(is_a, ev["star_icrs_ra_deg"], ev["relay_icrs_ra_deg"])

    # ---- 1. geometry-only era scope --------------------------------
    scope = {}
    for ch, rung, lim in RUNGS:
        m = (is_a if ch == "A" else ~is_a) & (b_au <= lim)
        north = m & (dec_pos > DEC_MIN)
        per = defaultdict(int)
        for tid in ev["target_id"][north]:
            per[str(tid)] += 1
        scope[f"{ch}_{rung}"] = {"events": int(m.sum()), "events_in_footprint": int(north.sum()),
                                 "targets_in_footprint": len(per), "per_target": dict(sorted(per.items()))}
    print("era scope:", {k: (v["events"], v["events_in_footprint"]) for k, v in scope.items()})

    # ---- 2. coverage proxy per target-channel ----------------------
    units = {}
    sel = (b_au <= 0.1) & (dec_pos > DEC_MIN)
    for i in np.where(sel)[0]:
        key = (str(ev["target_id"][i]), "A" if is_a[i] else "B")
        units.setdefault(key, []).append(int(i))
    print(f"{len(units)} target-channels with in-footprint events at <= 0.1 AU")

    results = {}
    for (tid, ch), idx in sorted(units.items()):
        ra = float(np.mean(ra_pos[idx])); dec = float(np.mean(dec_pos[idx]))
        srcs = q(store, f"SELECT pts_key,tmcra,tmcdec,tmcjmag,tmchmag,psfcontam,meanmag,rmsmag,"
                        f"q3c_dist(tmcra,tmcdec,{ra:.6f},{dec:.6f})*3600 AS sep_arcsec "
                        f"FROM pgir_dr1.sources WHERE q3c_radial_query(tmcra,tmcdec,{ra:.6f},{dec:.6f},{CONE_DEG}) "
                        f"ORDER BY sep_arcsec LIMIT 20")
        unit = {"channel": ra and ch, "ra_deg": ra, "dec_deg": dec, "n_events": len(idx),
                "n_sources_6arcmin": len(srcs),
                "nearest_sources": [{k: (float(r[k]) if k != "pts_key" else int(r[k])) for k in
                                     ("pts_key", "sep_arcsec", "tmcjmag", "meanmag", "rmsmag", "psfcontam")}
                                    for r in srcs[:3]]}
        epochs = {}
        for r in srcs[:N_PROXY]:
            rows = q(store, f"SELECT p.obsjd,p.stackquadid,p.magpsf,p.magpsferr,p.magpsflim,p.flags,"
                            f"e.limmag,e.nightid,e.exptime,e.filename FROM pgir_dr1.photometry p "
                            f"JOIN pgir_dr1.exposures e ON p.stackquadid=e.stackquadid "
                            f"WHERE p.pts_key={int(r['pts_key'])}")
            for x in rows:
                epochs.setdefault(int(x["stackquadid"]), x)
        mjd = np.array([float(x["obsjd"]) - 2400000.5 for x in epochs.values()])
        lim = np.array([float(x["limmag"]) for x in epochs.values()])
        nights = sorted({int(x["nightid"]) for x in epochs.values()})
        unit.update({"n_epochs": int(len(mjd)), "n_nights": len(nights),
                     "mjd_first": float(mjd.min()) if mjd.size else None,
                     "mjd_last": float(mjd.max()) if mjd.size else None,
                     "limmag_median": float(np.nanmedian(lim)) if lim.size else None,
                     "limmag_p10_p90": [float(np.nanpercentile(lim, 10)), float(np.nanpercentile(lim, 90))] if lim.size else None,
                     "distinct_obsjd_values": int(len(set(mjd.tolist())))})
        # per-event in-window counts per rung
        evs = []
        for i in idx:
            e = ev[i]
            row = {"event_id": str(e["event_id"]), "t_ca_utc": str(e["t_ca_utc"]),
                   "b_rsun": round(float(e["b_min_solar_radii"]), 3), "v_perp": round(float(e["v_perp_km_s"]), 1),
                   "rungs": {}}
            for rch, rung, r_au in RUNGS:
                if rch != ch or float(e["b_min_au"]) > r_au:
                    continue
                lo, hi = window(e, r_au)
                row["rungs"][rung] = {"half_width_d": round((hi - lo) / 2, 3), **count_modes(mjd, lo, hi)}
            closest = float(np.min(np.abs(mjd - float(e["t_ca_mjd"])))) if mjd.size else None
            row["closest_epoch_d"] = round(closest, 3) if closest is not None else None
            evs.append(row)
        unit["events"] = evs
        results[f"{tid}:{ch}"] = unit
        tot = {rung: sum(x["rungs"].get(rung, {}).get("nominal", 0) for x in evs) for _, rung, _ in RUNGS}
        print(f"  {tid}:{ch}  epochs={len(mjd)} nights={len(nights)} limmag~{unit['limmag_median']} "
              f"nearest={unit['nearest_sources'][0]['sep_arcsec']:.1f}\" in-window(nominal)={tot}")

    out = {"generated_utc": datetime.now(timezone.utc).isoformat(), "era_jd": ERA_JD, "era_mjd": ERA_MJD,
           "dec_min": DEC_MIN, "quantisation_d": QUANT_D, "cone_deg": CONE_DEG, "n_proxy_sources": N_PROXY,
           "events_input": str(EVENTS.relative_to(REPO)), "n_events_total": len(t_all),
           "n_events_in_era_not_sunward": len(ev), "era_scope": scope, "units": results}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print("wrote", OUT.relative_to(REPO))


if __name__ == "__main__":
    main()
