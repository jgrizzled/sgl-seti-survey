"""Snapshot the SoloHI L2 (and L1) holdings of the Solar Orbiter Archive
(SOAR TAP, anonymous): every data item's begin/end time, descriptor,
filename, size and item id -> runs/solohi-crossings/coverage/
soar_l2_inventory.json. Paged by month through v_sc_data_item."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import requests

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "runs" / "solohi-crossings" / "coverage"
TAP = "https://soar.esac.esa.int/soar-sl-tap/tap/sync"
COLS = "data_item_id, filename, descriptor, begin_time, end_time, filesize, is_public, generation_time"


def query(sql):
    for attempt in range(4):
        try:
            r = requests.get(TAP, params={"REQUEST": "doQuery", "LANG": "ADQL",
                                          "FORMAT": "json", "QUERY": sql}, timeout=600)
            r.raise_for_status()
            return r.json()["data"]
        except Exception as e:
            print("retry", attempt, e, file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(sql)


def main(level="L2"):
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for y in range(2020, 2027):
        for m in range(1, 13):
            a = f"{y}-{m:02d}-01T00:00:00"
            b = f"{y + (m == 12)}-{(m % 12) + 1:02d}-01T00:00:00"
            r = query(f"SELECT {COLS} FROM v_sc_data_item WHERE instrument='SOLOHI' "
                      f"AND level='{level}' AND begin_time >= '{a}' AND begin_time < '{b}' "
                      f"ORDER BY begin_time")
            rows.extend(r)
            print(a[:7], len(r), file=sys.stderr, flush=True)
    inv = {"level": level, "snapshot_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "columns": [c.strip() for c in COLS.split(",")], "rows": rows}
    (OUT / f"soar_{level.lower()}_inventory.json").write_text(json.dumps(inv) + "\n")
    print(level, len(rows), "items")


if __name__ == "__main__":
    main(*sys.argv[1:])
