"""Rubin Science Platform TAP adapter (snapshot-disciplined).

Sync ADQL against ``data.lsst.cloud/api/tap`` with the RSP_API_TOKEN
Bearer credential from ``.env``; verbatim CSV responses stored via
SnapshotStore (plan §8 reproducibility). The sync endpoint 303s to the
result — the session follows redirects.
"""

from __future__ import annotations

import csv
import io
from datetime import datetime, timezone
from pathlib import Path

import requests

TAP_SYNC_URL = "https://data.lsst.cloud/api/tap/sync"
REPO = Path(__file__).resolve().parents[2]


class RspTapError(RuntimeError):
    pass


def _load_token() -> str:
    for line in (REPO / ".env").read_text().splitlines():
        if line.startswith("RSP_API_TOKEN="):
            return line.split("=", 1)[1].strip().strip('"')
    raise RspTapError("RSP_API_TOKEN not found in .env")


class RspTapClient:
    """Snapshot-disciplined sync ADQL access to the RSP TAP service."""

    def __init__(self, timeout_s: float = 300.0):
        self.timeout_s = timeout_s
        self.session = requests.Session()
        self.session.headers["Authorization"] = f"Bearer {_load_token()}"

    def query(self, adql: str, store=None):
        """Run one sync ADQL query; return (columns, rows) from CSV.

        Snapshots the verbatim CSV when a SnapshotStore is given.
        Raises RspTapError on HTTP failure or a VOTable error payload.
        """
        request_utc = datetime.now(timezone.utc).isoformat()
        resp = self.session.post(
            TAP_SYNC_URL,
            data={"LANG": "ADQL", "REQUEST": "doQuery",
                  "FORMAT": "csv", "RESPONSEFORMAT": "csv",
                  "QUERY": adql},
            timeout=self.timeout_s, allow_redirects=True)
        cols, rows, ok = [], [], False
        if resp.status_code == 200 and resp.text.startswith("<?xml"):
            # the sync /run endpoint may ignore the requested format
            # and serve VOTable; accept it (and its error payloads)
            from astropy.io.votable import parse_single_table
            if 'value="ERROR"' not in resp.text[:2000]:
                tab = parse_single_table(io.BytesIO(resp.content),
                                         verify="ignore").to_table()
                cols = list(tab.colnames)
                rows = [["" if (hasattr(x, "mask") and x.mask) else str(
                    x.item() if hasattr(x, "item") else x) for x in r]
                    for r in tab]
                ok = True
        elif resp.status_code == 200:
            table = list(csv.reader(io.StringIO(resp.text)))
            cols, rows = (table[0], table[1:]) if table else ([], [])
            ok = True
        if store is not None:
            store.store(service_url=TAP_SYNC_URL, query=adql,
                        request_utc=request_utc,
                        response_bytes=resp.content,
                        row_count=len(rows) if ok else None,
                        http_status=resp.status_code)
        if not ok:
            raise RspTapError(
                f"TAP query failed (HTTP {resp.status_code}): "
                f"{resp.text[:300]} :: {adql[:200]}")
        return cols, rows
