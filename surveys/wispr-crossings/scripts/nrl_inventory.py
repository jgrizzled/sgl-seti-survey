"""Snapshot the NRL WISPR L2 (and L3) trees: every day directory's file
list with sizes and modification stamps -> runs/wispr-crossings/coverage/
l2_inventory.json (and l3_inventory.json). Anonymous plain-HTTP listings;
~1,000 requests per level."""
from __future__ import annotations
import json, re, sys, time
from pathlib import Path
import requests

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "runs" / "wispr-crossings" / "coverage"
BASE = "https://wispr.nrl.navy.mil/data/rel/fits"
PAT = re.compile(r'href="(psp_L\w+_wispr_(\d{8}T\d{6})_V(\d)_(\d{4})\.fits)"[^\n]*?(\d{4}-\d{2}-\d{2} \d{2}:\d{2})\s+</td><td[^>]*>\s*([\d.]+[KMG]?)')


def listing(url):
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.text
        except Exception as e:
            time.sleep(2 * (attempt + 1))
    return None


def parse(html):
    rows = []
    for m in re.finditer(r'href="(psp_L\w+_wispr_(\d{8}T\d{6})_V(\d)_(\d{4})\.fits)"', html):
        rows.append({"file": m.group(1), "stamp": m.group(2), "version": int(m.group(3)),
                     "wxyz": m.group(4)})
    return rows


def main(level="L2"):
    OUT.mkdir(parents=True, exist_ok=True)
    inv = {"level": level, "snapshot_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "days": {}}
    if level == "L2":
        top = listing(f"{BASE}/L2/")
        days = sorted(set(re.findall(r'href="(\d{8})/"', top)))
        for i, d in enumerate(days):
            html = listing(f"{BASE}/L2/{d}/")
            inv["days"][d] = parse(html) if html else None
            if i % 50 == 0:
                print(i, d, len(inv["days"][d] or []), file=sys.stderr, flush=True)
    else:
        top = listing(f"{BASE}/L3/")
        orbits = sorted(set(re.findall(r'href="(orbit\d+)/"', top)))
        for o in orbits:
            html = listing(f"{BASE}/L3/{o}/")
            for d in sorted(set(re.findall(r'href="(\d{8})/"', html or ""))):
                h = listing(f"{BASE}/L3/{o}/{d}/")
                inv["days"][f"{o}/{d}"] = parse(h) if h else None
            print(o, file=sys.stderr, flush=True)
    (OUT / f"{level.lower()}_inventory.json").write_text(json.dumps(inv) + "\n")
    n = sum(len(v or []) for v in inv["days"].values())
    print(level, len(inv["days"]), "days,", n, "files")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "L2")
