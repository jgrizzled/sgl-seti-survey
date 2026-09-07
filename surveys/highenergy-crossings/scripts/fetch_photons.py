"""Era-long Fermi-LAT photon pulls per channel position (hypotheses.md
freeze v1.0 §5/D3) through the LAT data server, plus the two positive-
control pulls (D7). Submits, polls, downloads, and records every query
under runs/highenergy-crossings/v1/photons/<key>/ with the submit and
result pages snapshotted. Idempotent: a key with a complete manifest is
skipped."""
from __future__ import annotations

import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from sglsurvey.snapshots import SnapshotStore  # noqa: E402
import recon_scan as R  # noqa: E402

RUN = REPO / "runs" / "highenergy-crossings" / "v1"
PHOT = RUN / "photons"
QUERY_URL = "https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi"
RESULT_URL = "https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/QueryResults.cgi"
ERA = ("2008-08-04 15:43:36", "2026-09-07 00:13:13")
RADIUS_DEG = 3.0
ENERGY = "100,300000"

CONTROLS = {
    "ctrl-grb130427a": (173.136, 27.699, "2012-12-27 00:00:00", "2013-08-27 00:00:00"),
    "ctrl-3c454.3": (343.491, 16.148, "2010-07-20 00:00:00", "2011-03-22 00:00:00"),
}


def now():
    return datetime.now(timezone.utc).isoformat()


def submit(store, ra, dec, t0, t1):
    data = {"coordfield": f"{ra:.4f},{dec:.4f}", "coordsystem": "J2000",
            "shapefield": f"{RADIUS_DEG:.1f}", "timefield": f"{t0},{t1}",
            "timetype": "Gregorian", "energyfield": ENERGY,
            "photonOrExtendedOrNone": "Photon", "spacecraft": "off",
            "destination": "query"}
    r = requests.post(QUERY_URL, data=data, timeout=180)
    store.store(service_url=QUERY_URL, query=json.dumps(data), request_utc=now(),
                response_bytes=r.content, row_count=None, http_status=r.status_code)
    m = re.search(r"QueryResults\.cgi\?id=([A-Z0-9]+)", r.text)
    if not m:
        raise RuntimeError("no query id in response: " + r.text[:300])
    return m.group(1), data


def poll(store, qid, timeout_s=3600):
    t0 = time.monotonic()
    while time.monotonic() - t0 < timeout_s:
        r = requests.get(RESULT_URL, params={"id": qid}, timeout=120)
        files = sorted(set(re.findall(r"queries/(" + qid + r"_(?:PH|SC)\d+\.fits)", r.text)))
        # the results page says "The state of your query is 2 (Query complete)"
        # and lists "Query complete" per server; any other state keeps polling
        done = bool(files) and ("(Query complete)" in r.text) and (r.text.count("Query complete") >= 2)
        if done:
            store.store(service_url=RESULT_URL, query=f"id={qid}", request_utc=now(),
                        response_bytes=r.content, row_count=len(files), http_status=r.status_code)
            return files
        time.sleep(20)
    raise RuntimeError(f"query {qid} not complete after {timeout_s} s")


def fetch(key, ra, dec, t0, t1, store):
    d = PHOT / key
    man = d / "manifest.json"
    if man.exists():
        print(f"[fetch] {key}: present", flush=True)
        return json.loads(man.read_text())
    d.mkdir(parents=True, exist_ok=True)
    qid, data = submit(store, ra, dec, t0, t1)
    print(f"[fetch] {key}: query {qid} submitted", flush=True)
    files = poll(store, qid)
    out = []
    for f in files:
        if "_SC" in f:      # the era-long spacecraft file is GB-scale; the weekly files are the frozen source
            continue
        url = f"https://fermi.gsfc.nasa.gov/FTP/fermi/data/lat/queries/{f}"
        for k in range(3):
            r = requests.get(url, timeout=900)
            if r.status_code == 200 and len(r.content) > 2880:
                (d / f).write_bytes(r.content)
                out.append({"file": f, "bytes": len(r.content),
                            "sha256": __import__("hashlib").sha256(r.content).hexdigest()})
                break
            time.sleep(15)
    m = {"key": key, "query_id": qid, "params": data, "ra": ra, "dec": dec,
         "utc": now(), "files": out}
    man.write_text(json.dumps(m, indent=1))
    print(f"[fetch] {key}: {len(out)} files, {sum(x['bytes'] for x in out)/1e6:.0f} MB", flush=True)
    return m


def main():
    store = SnapshotStore(RUN / "photons")
    events = R.load_events()
    pos = R.positions(events)
    jobs = []
    for (ch, tid), (ra, de, n) in sorted(pos.items()):
        jobs.append((f"{ch}-{tid}", ra, de, ERA[0], ERA[1]))
    for key, (ra, de, t0, t1) in CONTROLS.items():
        jobs.append((key, ra, de, t0, t1))
    # sequential submit/poll keeps the server load at one query at a time
    for key, ra, de, t0, t1 in jobs:
        for attempt in range(3):
            try:
                fetch(key, ra, de, t0, t1, store)
                break
            except Exception as exc:
                print(f"[fetch] {key}: attempt {attempt+1} failed: {exc}", file=sys.stderr, flush=True)
                time.sleep(60)
    (PHOT / "positions.json").write_text(json.dumps(
        {f"{ch}-{tid}": {"ra": ra, "dec": de, "n_events": n} for (ch, tid), (ra, de, n) in pos.items()}
        | {k: {"ra": v[0], "dec": v[1], "t0": v[2], "t1": v[3]} for k, v in CONTROLS.items()}, indent=1))


if __name__ == "__main__":
    main()
