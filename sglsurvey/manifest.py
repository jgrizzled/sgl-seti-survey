"""Content-hash manifests and stale-product detection (WISE v2 plan
§8.1–8.2, review §9).

A manifest lists every input file of a run by sha256 of its bytes plus
the schema version of the stage that wrote it; the manifest's own hash
is what AnalysisRun records carry in ``observation_set_hash`` /
``intersection_set_hash``. Products that derive from a set of inputs
store ``input_hash`` (the hash of the sorted input digests) so that a
rebuild with ``--only-missing`` can refuse a product whose recorded
inputs no longer match the files on disk.

Hashing tens of thousands of cutouts is cached by (path, size, mtime)
in a sidecar JSON so that repeated runs cost seconds, not minutes;
the cache is a convenience, the digests are the truth.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Iterable, Mapping

SCHEMA_VERSION = "manifest-v1"


def file_sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


class HashCache:
    """(path, size, mtime_ns) -> sha256 cache persisted as JSON."""

    def __init__(self, cache_path: Path):
        self.path = cache_path
        self._d: dict = {}
        if cache_path.exists():
            try:
                self._d = json.loads(cache_path.read_text())
            except json.JSONDecodeError:
                self._d = {}
        self._dirty = 0

    def sha256(self, path: Path) -> str:
        st = os.stat(path)
        key = str(path)
        rec = self._d.get(key)
        if rec and rec["size"] == st.st_size and rec["mtime_ns"] == st.st_mtime_ns:
            return rec["sha256"]
        digest = file_sha256(path)
        self._d[key] = {"size": st.st_size, "mtime_ns": st.st_mtime_ns,
                        "sha256": digest}
        self._dirty += 1
        if self._dirty >= 500:
            self.flush()
        return digest

    def flush(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self._d))
        os.replace(tmp, self.path)
        self._dirty = 0


def combined_hash(digests: Iterable[str]) -> str:
    """Order-independent hash of a set of digests."""
    h = hashlib.sha256()
    for d in sorted(digests):
        h.update(d.encode())
    return "sha256:" + h.hexdigest()


def build_manifest(files: Mapping[str, Path], *, stage: str,
                   schema_version: str, cache: HashCache | None = None,
                   extra: Mapping | None = None) -> dict:
    """Manifest dict for a named set of files. ``files`` maps a logical
    name (e.g. an observation id + product kind) to a path."""
    entries = {}
    for name, path in files.items():
        path = Path(path)
        digest = cache.sha256(path) if cache else file_sha256(path)
        entries[name] = {"path": str(path), "sha256": digest,
                         "bytes": os.path.getsize(path)}
    if cache:
        cache.flush()
    man = {"manifest_schema": SCHEMA_VERSION, "stage": stage,
           "schema_version": schema_version, "n_files": len(entries),
           "files": entries, "extra": dict(extra or {})}
    man["manifest_hash"] = combined_hash(
        [e["sha256"] for e in entries.values()] + [schema_version])
    return man


def write_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=1, sort_keys=True))


def load_manifest(path: Path) -> dict:
    return json.loads(Path(path).read_text())


class StaleProductError(RuntimeError):
    pass


def check_product_inputs(recorded_input_hash: str, current_digests: Iterable[str],
                         product: str) -> None:
    """Raise :class:`StaleProductError` if a product's recorded input hash
    differs from the hash of the inputs as they are now on disk."""
    now = combined_hash(current_digests)
    if recorded_input_hash != now:
        raise StaleProductError(
            f"{product}: recorded input hash {recorded_input_hash[:23]}… "
            f"!= current {now[:23]}…; rebuild it (do not --only-missing)")
