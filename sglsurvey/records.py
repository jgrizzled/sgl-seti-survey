"""The five core survey records (plan §3.2) plus query snapshots.

Records are immutable, content-addressed via sglseti's provenance
helpers, and stored append-only as JSON Lines. Identity hashes cover
science-defining fields only; run metadata (wall-clock times, local
paths) is carried on the record but excluded from the ID. Design:
docs/design/records_and_adapters.md.
"""

from __future__ import annotations

import dataclasses
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from sglseti import canonical_json, stable_id

Json = Mapping[str, Any]


def _clean(value: Any) -> Any:
    """Drop None entries from mappings for stable identity payloads."""
    if isinstance(value, Mapping):
        return {k: _clean(v) for k, v in value.items() if v is not None}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    return value


@dataclass(frozen=True)
class QuerySnapshot:
    """Verbatim capture of one archive service request/response."""

    snapshot_id: str
    service_url: str
    query: str
    request_utc: str  # run metadata, excluded from ID
    response_sha256: str
    row_count: int | None
    http_status: int
    response_path: str | None = None  # local convention, excluded from ID

    @classmethod
    def build(cls, *, service_url: str, query: str, request_utc: str,
              response_sha256: str, row_count: int | None,
              http_status: int, response_path: str | None = None,
              ) -> "QuerySnapshot":
        sid = stable_id("snp", {
            "service_url": service_url,
            "query": query,
            "response_sha256": response_sha256,
        })
        return cls(sid, service_url, query, request_utc, response_sha256,
                   row_count, http_status, response_path)


@dataclass(frozen=True)
class Observation:
    """Immutable archive metadata and data-product identity."""

    observation_id: str
    archive_id: str
    collection: str
    release: str
    native_key: Json
    band: str
    wavelength_um: float | None
    t_start_mjd_utc: float
    t_mid_mjd_utc: float
    t_stop_mjd_utc: float
    exptime_s: float
    corners_icrs_deg: tuple[tuple[float, float], ...]
    wcs: Json
    quality_flags: Json = field(default_factory=dict)
    products: Json = field(default_factory=dict)
    extra: Json = field(default_factory=dict)
    snapshot_id: str | None = None
    discovered_utc: str | None = None  # run metadata, excluded from ID

    @classmethod
    def build(cls, *, archive_id: str, collection: str, release: str,
              native_key: Json, **fields: Any) -> "Observation":
        oid = stable_id("obs", {
            "archive_id": archive_id,
            "collection": collection,
            "release": release,
            "native_key": _clean(native_key),
        })
        return cls(oid, archive_id, collection, release, native_key,
                   **fields)


@dataclass(frozen=True)
class IntersectionEvaluation:
    """One observation x one hypothesis cell x one evaluation stage.

    Coarse misses are retained as audit records; a precise hit is still
    not coverage until usability and an analysis run exist.
    """

    intersection_id: str
    observation_id: str
    hypothesis_version: str
    hypothesis_hash: str
    endpoint_id: str
    role: str  # "rx" | "tx"
    registry_source_hash: str
    target_source_hash: str
    model_id: str
    model_version: str
    ephemeris_id: str
    tolerance_arcsec: float
    confidence_level: float
    padding_arcsec: float
    stage: str  # "coarse" | "precise"
    hit: bool
    covered_z_intervals_au: tuple[tuple[float, float], ...] = ()
    envelope_pad_arcsec: float | None = None
    mc_seed: int | None = None
    mc_sample_count: int | None = None
    warnings: tuple[str, ...] = ()
    usable: str = "unknown"  # "unknown" | "usable" | "unusable" | "partial"
    usable_fraction: float | None = None
    usability_notes: str | None = None
    extra: Json = field(default_factory=dict)

    @classmethod
    def build(cls, **fields: Any) -> "IntersectionEvaluation":
        identity_fields = {
            k: v for k, v in fields.items()
            if k not in ("usability_notes", "extra") and v is not None
        }
        iid = stable_id("ixn", _clean(identity_fields))
        return cls(intersection_id=iid, **fields)


@dataclass(frozen=True)
class AnalysisRun:
    """A versioned pipeline actually searched specified data."""

    analysis_run_id: str
    pipeline_id: str
    pipeline_version: str
    config: Json
    observation_set_hash: str
    intersection_set_hash: str
    registry_source_hash: str
    hypothesis_version: str
    environment: Json  # sglseti/sglsurvey commits + dirty flags, versions
    random_seeds: Json = field(default_factory=dict)
    # run metadata, excluded from ID:
    started_utc: str | None = None
    finished_utc: str | None = None
    output_files: Json = field(default_factory=dict)
    warning_summary: tuple[str, ...] = ()


@dataclass(frozen=True)
class Constraint:
    """Injection-calibrated sensitivity or a qualified null result."""

    constraint_id: str
    analysis_run_id: str
    endpoint_id: str
    role: str
    hypothesis_version: str
    z_interval_au: tuple[float, float]
    band: str
    epoch_range_mjd: tuple[float, float]
    duty_cycle_range: tuple[float, float]
    residual_motion_bound_arcsec_per_yr: float
    morphology: str
    kind: str  # "recovery_curve" | "qualified_null" | "not_constrainable"
    recovery_probability: float | None = None
    flux_limit: Json | None = None  # value, unit, physical interpretation
    injection_summary_ref: str | None = None
    false_alarm_rate: float | None = None
    trials_accounting_ref: str | None = None
    extra: Json = field(default_factory=dict)


@dataclass(frozen=True)
class Candidate:
    """A retained event or track with its competing-model tests."""

    candidate_id: str
    analysis_run_id: str
    endpoint_id: str
    role: str
    observation_ids: tuple[str, ...]
    fitted_z_au: float | None = None
    fitted_z_interval_au: tuple[float, float] | None = None
    fitted_residual_motion: Json | None = None
    measurements_ref: str | None = None
    model_comparison: Json = field(default_factory=dict)
    status: str = "new"
    holdout_test_ref: str | None = None
    veto_reason: str | None = None
    extra: Json = field(default_factory=dict)


@dataclass(frozen=True)
class ScreenMatch:
    """A catalog detection within the screen radius of a predicted SGL
    track at its epoch (plan §4.5 layer-1 screening). One record per
    (detection, endpoint, role) pair; catalog absence or presence alone
    is neither a detection nor a null result."""

    screen_match_id: str
    endpoint_id: str
    role: str
    hypothesis_version: str
    registry_source_hash: str
    source_table: str
    source_cntr: str
    source_designation: str | None
    ra_deg: float
    dec_deg: float
    sigra_mas: float | None
    sigdec_mas: float | None
    mjd: float
    scan_id: str | None
    frame_num: int | None
    photometry: Json
    flags: Json
    dist_arcsec: float
    implied_z_au: float
    z_segment_au: tuple[float, float]
    screen_radius_arcsec: float
    locus_mjd: float
    snapshot_id: str | None = None
    extra: Json = field(default_factory=dict)

    @classmethod
    def build(cls, **fields: Any) -> "ScreenMatch":
        sid = stable_id("scr", {
            "endpoint_id": fields["endpoint_id"],
            "role": fields["role"],
            "hypothesis_version": fields["hypothesis_version"],
            "registry_source_hash": fields["registry_source_hash"],
            "source_table": fields["source_table"],
            "source_cntr": fields["source_cntr"],
        })
        return cls(screen_match_id=sid, **fields)


def to_jsonl_line(record: Any) -> str:
    """Serialize a record dataclass to one canonical JSON line."""
    return canonical_json(dataclasses.asdict(record))


def append_records(path: Path, records: Sequence[Any]) -> None:
    """Append records to a JSONL file, creating parents as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        for record in records:
            fh.write(to_jsonl_line(record) + "\n")


def read_records(path: Path) -> list[dict[str, Any]]:
    """Read raw record dicts back from a JSONL file."""
    with path.open("r", encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]
