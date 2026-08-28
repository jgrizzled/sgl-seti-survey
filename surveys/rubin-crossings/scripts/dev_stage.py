"""Rubin DP2 crossings dev stage (pseudo-units only, freeze D8).

Runs the identical frozen chain (rubinlib) on every off-window visit
at the ross-128 antipode field — epochs where no relay is predicted —
to measure the association-rate null, exercise the control rule, run
the SkyBoT known-object positive control, and census the archive
injection population. The in-window unit visit 2025091300612 is
NEVER queried here (blindness guard asserts it).

Outputs: results/dev_v1.json (+ prints); snapshots under
runs/rubin-crossings/dev_v1/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import rubinlib as rl  # noqa: E402

from sglsurvey.adapters.rsp_tap import RspTapClient  # noqa: E402
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

REPO = rl.REPO
OUT = REPO / "surveys" / "rubin-crossings" / "results" / "dev_v1.json"
RUN = REPO / "runs" / "rubin-crossings" / "dev_v1"

EVENT_ID = "evt-5197f21cd67b"
UNIT_VISIT = 2025091300612          # blind — must never be queried at dev
WINDOW_MJD = (60932.10, 60943.78)   # the 0.1 AU rung window (coverage v1)
SKYBOT = "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"
TAI_MINUS_UTC_S = 37.0


class BlindGuardClient(RspTapClient):
    def query(self, adql, store=None):
        assert str(UNIT_VISIT) not in adql, "blindness violation"
        return super().query(adql, store)


def main():
    ev = None
    t = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
    ev = t[np.array(t["event_id"]) == EVENT_ID][0]
    client = BlindGuardClient()
    store = SnapshotStore(RUN)
    pivot = rl.axis_point(ev)
    out = {"event_id": EVENT_ID, "pivot": pivot}

    # -- pseudo-unit family: off-window visits at the field ------------
    _, vrows = client.query(
        "SELECT visit, band, expMidptMJD FROM dp2.Visit WHERE "
        "CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', "
        f"{pivot[0]:.6f}, {pivot[1]:.6f}, 2.0)) = 1", store)
    visits = [(int(r[0]), str(r[1]), float(r[2])) for r in vrows]
    off = [(v, b, m) for v, b, m in visits
           if not (WINDOW_MJD[0] <= m <= WINDOW_MJD[1])]
    inw = [v for v, b, m in visits
           if WINDOW_MJD[0] <= m <= WINDOW_MJD[1]]
    print(f"{len(visits)} field visits; {len(off)} off-window pseudo-"
          f"units; in-window (excluded, blind): {inw}")
    assert inw == [UNIT_VISIT]

    objects = rl.fetch_static_objects(client, store, pivot,
                                      rl.CONE_ARCSEC + 30.0)
    out["n_static_objects_90as"] = len(objects)

    records = []
    for v, b, m in sorted(off, key=lambda x: x[2]):
        rec = rl.run_chain(ev, v, m, b, client, store, objects=objects)
        rec["band"] = b
        records.append(rec)
        if rec["status"] == "searched":
            print(f"  {v} {b} mjd {m:.3f}: S={rec['s_det']:.2f} "
                  f"T={rec['threshold']:.2f} "
                  f"exc={rec['exceedance']} ndia={rec['n_diasources_cone']}")
        else:
            print(f"  {v} {b}: {rec['status']}")
    searched = [r for r in records if r["status"] == "searched"]
    n_exc = sum(r["exceedance"] for r in searched)
    out["pseudo_units"] = records
    out["null_summary"] = {
        "n_pseudo_units": len(records),
        "n_searched": len(searched),
        "n_exceedances": n_exc,
        "expected": len(searched) / 9.0,
        "assoc_rate_s_gt0": (sum(r["s_det"] > 0 for r in searched)
                             / max(len(searched), 1)),
        "extra_rotations": sum(c["extra_rotation_deg"] > 0
                               for r in searched
                               for c in r["control_census"]),
    }
    print("null:", json.dumps(out["null_summary"], indent=1))

    # -- positive control: SkyBoT-predicted known object ---------------
    # find ssObject-linked DiaSources on off-window r visits (field-wide)
    pc = {"attempts": []}
    for v, b, m in sorted(off, key=lambda x: x[2]):
        if b != "r" or pc["attempts"] and pc["attempts"][-1]["pass"]:
            continue
        _, rows = client.query(
            "SELECT TOP 5 diaSourceId, ssObjectId, ra, dec, psfFlux, "
            "psfFluxErr FROM dp2.DiaSource WHERE "
            f"visit = {v} AND ssObjectId IS NOT NULL AND ssObjectId != 0 "
            "AND CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', "
            f"{pivot[0]:.6f}, {pivot[1]:.6f}, 1.0)) = 1", store)
        if not rows:
            continue
        ds = rows[0]
        ra_d, de_d = float(ds[2]), float(ds[3])
        jd_utc = m - TAI_MINUS_UTC_S / 86400.0 + 2400000.5
        resp = requests.get(SKYBOT, params={
            "-ep": f"{jd_utc:.8f}", "-ra": f"{ra_d:.6f}",
            "-dec": f"{de_d:.6f}", "-rd": "0.05", "-mime": "text",
            "-output": "object", "-loc": "X05"}, timeout=60)
        store.store(service_url=SKYBOT,
                    query=f"ep={jd_utc:.8f} ra={ra_d:.6f} dec={de_d:.6f}",
                    request_utc="2026-08-26",
                    response_bytes=resp.content, row_count=None,
                    http_status=resp.status_code)
        lines = [ln for ln in resp.text.splitlines()
                 if ln and not ln.startswith(("#", "-"))]
        best = None
        for ln in lines:
            parts = [p.strip() for p in ln.split("|")]
            if len(parts) < 4:
                continue
            try:
                from astropy.coordinates import Angle
                sra = Angle(parts[2] + " hours").degree
                sde = Angle(parts[3] + " degrees").degree
            except Exception:
                continue
            d = rl.sep_arcsec(ra_d, de_d, sra, sde)
            if best is None or d < best[1]:
                best = ((parts[1], sra, sde), d)
        att = {"visit": v, "mjd": m, "diaSourceId": ds[0],
               "ssObjectId": ds[1], "dia_pos": [ra_d, de_d],
               "n_skybot": len(lines), "pass": False}
        if best:
            (name, sra, sde), doff = best
            dias = rl.fetch_diasources(client, store, v, (sra, sde),
                                       30.0)
            s, hit = rl.s_det([(sra, sde)], dias)
            att.update({"skybot_name": name,
                        "skybot_offset_arcsec": round(doff, 3),
                        "s_det_at_prediction": s,
                        "associated": bool(hit),
                        "pass": bool(hit)
                        and hit["diaSourceId"] == ds[0]})
        pc["attempts"].append(att)
        print("positive control:", json.dumps(att, default=str))
        if att["pass"]:
            break
    pc["pass"] = any(a.get("pass") for a in pc["attempts"])
    out["positive_control"] = pc

    # -- archive-injection census (off-window r visits, field cone) ----
    inj = []
    for v, b, m in off:
        if b != "r":
            continue
        try:
            _, rows = client.query(
                "SELECT psfFlux, psfFluxErr FROM dp2.DiaSource WHERE "
                f"visit = {v} AND (pixelFlags_injectedCenter = 1 OR "
                "pixelFlags_injected = 1) AND CONTAINS(POINT('ICRS', "
                f"ra, dec), CIRCLE('ICRS', {pivot[0]:.6f}, "
                f"{pivot[1]:.6f}, 0.5)) = 1", store)
        except Exception as exc:
            out["injection_census_error"] = str(exc)[:300]
            break
        for r in rows:
            if r[0] not in ("", None):
                f = float(r[0])
                inj.append({"visit": v, "psfFlux": f,
                            "snr": (f / float(r[1])
                                    if r[1] not in ("", None) else None),
                            "mag": (31.4 - 2.5 * np.log10(abs(f))
                                    if f else None)})
    mags = sorted(x["mag"] for x in inj if x["mag"])
    out["injection_census"] = {
        "n_injected_diasources_r": len(inj),
        "n_visits_probed": len([1 for v, b, m in off if b == "r"]),
        "mag_range": ([round(mags[0], 2), round(mags[-1], 2)]
                      if mags else None),
        "mag_median": round(float(np.median(mags)), 2) if mags else None,
    }
    print("injections:", json.dumps(out["injection_census"]))

    OUT.write_text(json.dumps(out, indent=1, default=float))
    print("wrote", OUT.relative_to(REPO))


if __name__ == "__main__":
    main()
