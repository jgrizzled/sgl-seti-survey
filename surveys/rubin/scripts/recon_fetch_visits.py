"""Rubin DP2 recon: snapshot dp2.Visit (async TAP) + column inventories.

Writes:
  runs/rubin/recon/dp2_visit_snapshot.ecsv     full visit table
  runs/rubin/recon/dp2_columns.json            column lists for key tables
Auth: RSP_API_TOKEN from .env (x-oauth-basic basic auth per rsp.lsst.io).
"""
from __future__ import annotations

import json
import os
import pathlib

import pyvo as vo
from pyvo.auth import AuthSession

REPO = pathlib.Path(__file__).resolve().parents[3]
OUT = REPO / "runs" / "rubin" / "recon"
TAP_URL = "https://data.lsst.cloud/api/tap"


def get_service() -> vo.dal.TAPService:
    token = None
    for line in (REPO / ".env").read_text().splitlines():
        if line.startswith("RSP_API_TOKEN="):
            token = line.split("=", 1)[1].strip().strip('"')
    assert token, "RSP_API_TOKEN not found in .env"
    session = AuthSession()
    session.credentials.set_password("x-oauth-basic", token)
    session.add_security_method_for_url(TAP_URL, "ivo://ivoa.net/sso#BasicAA")
    session.add_security_method_for_url(TAP_URL + "/sync", "ivo://ivoa.net/sso#BasicAA")
    session.add_security_method_for_url(TAP_URL + "/async", "ivo://ivoa.net/sso#BasicAA")
    session.add_security_method_for_url(TAP_URL + "/tables", "ivo://ivoa.net/sso#BasicAA")
    return vo.dal.TAPService(TAP_URL, session=session)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    svc = get_service()

    # --- column inventories (sync; small) ---
    cols: dict[str, list[str]] = {}
    for table in ("dp2.Visit", "dp2.VisitDetector", "dp2.DiaSource",
                  "dp2.DiaObject", "dp2.ForcedSource", "dp2.Object",
                  "dp2.SSSource", "dp2.SSObject"):
        res = svc.search(
            "SELECT column_name FROM tap_schema.columns "
            f"WHERE table_name = '{table}'")
        cols[table] = sorted(str(r["column_name"]) for r in res)
        print(f"{table}: {len(cols[table])} columns")
    (OUT / "dp2_columns.json").write_text(json.dumps(cols, indent=1))

    # --- full visit table (async; probes the UWS route) ---
    query = ("SELECT visit, ra, dec, band, expMidptMJD, expTime, "
             "airmass, skyRotation FROM dp2.Visit")
    job = svc.submit_job(query)
    job.run()
    job.wait(phases=["COMPLETED", "ERROR"])
    print("async phase:", job.phase)
    tab = job.fetch_result().to_table()
    print("visits:", len(tab))
    tab.meta["recon"] = {"query": query, "tap_url": TAP_URL,
                         "retrieved_utc": "2026-08-26"}
    tab.write(OUT / "dp2_visit_snapshot.ecsv", format="ascii.ecsv",
              overwrite=True)
    job.delete()


if __name__ == "__main__":
    main()
