"""Verbatim storage of archive query responses (plan §6 reproducibility:
snapshot queries and responses; archives evolve even when survey code
does not)."""

from __future__ import annotations

import hashlib
from pathlib import Path

from sglsurvey.records import QuerySnapshot, append_records


class SnapshotStore:
    """Writes response bytes and QuerySnapshot records under a run dir.

    Layout: ``<run_dir>/snapshots/<snapshot_id>.bin`` plus an
    append-only ``<run_dir>/records/query_snapshot.jsonl``.
    """

    def __init__(self, run_dir: Path):
        self.run_dir = Path(run_dir)
        self.snapshot_dir = self.run_dir / "snapshots"
        self.record_path = self.run_dir / "records" / "query_snapshot.jsonl"

    def store(self, *, service_url: str, query: str, request_utc: str,
              response_bytes: bytes, row_count: int | None,
              http_status: int) -> QuerySnapshot:
        digest = hashlib.sha256(response_bytes).hexdigest()
        snapshot = QuerySnapshot.build(
            service_url=service_url,
            query=query,
            request_utc=request_utc,
            response_sha256=f"sha256:{digest}",
            row_count=row_count,
            http_status=http_status,
        )
        self.snapshot_dir.mkdir(parents=True, exist_ok=True)
        payload_path = self.snapshot_dir / f"{snapshot.snapshot_id}.bin"
        if not payload_path.exists():
            payload_path.write_bytes(response_bytes)
        snapshot = QuerySnapshot(
            **{**snapshot.__dict__,
               "response_path": str(payload_path.relative_to(self.run_dir))},
        )
        append_records(self.record_path, [snapshot])
        return snapshot
