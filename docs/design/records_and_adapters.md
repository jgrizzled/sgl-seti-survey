---
title: "Core records and archive adapter interfaces"
status: "Draft v0.1 — 2026-08-18"
---

# Core records and archive adapter interfaces

Design for the shared `sglsurvey` package: the five core record types
(plan §3.2, §8) and the adapter boundary that keeps archive-specific
logic out of the pipeline. Grounded in the sglseti v1.1 API
(`../sglseti`, pinned `21f6f3d`).

## Principles

1. **Records are immutable and content-addressed.** We reuse
   `sglseti.provenance` (`stable_hash`, `stable_id`, `canonical_json`,
   `file_sha256`) rather than inventing a second identity scheme —
   `import sglseti` is side-effect free and these helpers are public.
   Same rules apply: IDs never contain filesystem paths; new
   default-valued fields don't move existing IDs.
2. **One record type per state.** Discovery, geometry, analysis, and
   interpretation never share a table (plan §3.2). Records reference each
   other by ID only.
3. **Raw alongside normalized.** Every discovery response is stored
   verbatim as a snapshot; normalized records point at the snapshot they
   were derived from (plan §8 "archive adapters").
4. **The pipeline sees protocols, not archives.** Discovery, product
   download, and valid-pixel logic live behind `ArchiveAdapter`; sglseti
   sees footprints only as an opaque `contains(ra_deg, dec_deg) -> bool`
   callback, which is exactly the signature `covered_z_intervals`
   consumes.

## Storage model

Append-only JSON Lines per record type per run, under
`runs/<archive>/<run-id>/records/<type>.jsonl`, with a `manifest.json`
per run (sglseti-style: science-input hashes separated from run
metadata). Tracked, compact roll-ups (counts, coverage tables,
constraint summaries) are derived into `surveys/<archive>/results/`.
Parquet/ECSV consolidation can come later without changing record
identity, because IDs hash canonical content, not files.

## The five records

Field lists below are the normative core; `extra: dict` on each record
carries archive-specific normalized fields without disturbing identity
(excluded from the ID hash only where marked).

### Observation — `obs-<hash12>`

Immutable identity + metadata of one archive data product. One record
per (archive product row); for WISE L1b, per (mission phase table,
`scan_id`, `frame_num`, `band`, `wrelease`).

- `observation_id` — `stable_id("obs", {archive_id, collection, native_key, release})`
- `archive_id` (e.g. `irsa-wise`), `collection` (e.g.
  `wise.neowiser_merge_p1bm_frm`), `release` (e.g. `wrelease`)
- `native_key` — archive-native identifying mapping
  (`{scan_id, frame_num, band}`)
- `band` — normalized band label (`W1`…`W4`) + `wavelength_um`
- timing: `t_start_mjd_utc`, `t_mid_mjd_utc`, `t_stop_mjd_utc`,
  `exptime_s` (start/stop derived from mid ± exptime/2 for WISE)
- footprint: `corners_icrs_deg` (4×2 nominal polygon) and `wcs`
  (full native WCS mapping — CD matrix, CRVAL/CRPIX, projection type,
  distortion flag) — nominal only; valid-pixel truth requires products
- quality: normalized `quality_flags` mapping (frame quality, moon
  separation, …)
- products: mapping product kind → `{url, checksum_md5|sha256, size}`
  (WISE: `int`, `msk`, `unc`, `art`, `fflag`; IBE publishes `.md5`
  sidecars — record them at discovery time)
- provenance: `snapshot_id` of the discovery query, `discovered_utc`
  (run metadata, excluded from ID)

### IntersectionEvaluation — `ixn-<hash12>`

One observation × one hypothesis cell × one evaluation stage. Coarse
misses are retained as audit records (plan §3.2).

- `intersection_id` — `stable_id("ixn", {...all science fields...})`
- refs: `observation_id`
- hypothesis cell: `hypothesis_version` (+ content hash of the frozen
  hypotheses file), `endpoint_id` (registry target ID), `role`
  (`rx`|`tx`), `registry_source_hash`, `target_source_hash`
- geometry pins: `model_id`, `model_version`, `ephemeris_id`,
  `observer_provider` triple, `tolerance_arcsec`,
  `confidence_level`, `padding_arcsec`, MC `seed`/`sample_count`
  where uncertainty was propagated
- `stage` — `coarse` | `precise`
- result: `hit: bool`; for precise hits `covered_z_intervals_au`
  (list of `[z_min, z_max]` from `covered_z_intervals`, possibly
  disjoint), `envelope_pad_arcsec`, sglseti `warnings` passed through
- usability: `usable` (`unknown` at coarse stage; at precise stage the
  result of mask/quality tests), `usable_fraction`, `usability_notes`

### AnalysisRun — `run-<hash12>`

A versioned pipeline actually searched specified data.

- `analysis_run_id` — `stable_id("run", {science config})`
- `pipeline_id`, `pipeline_version`, full `config` (canonicalized),
  decision thresholds
- inputs: sorted `observation_ids` hash, sorted `intersection_ids` hash,
  `registry_source_hash`, `hypothesis_version`
- environment: `sglseti_commit` + dirty flag, `sglsurvey_commit` + dirty
  flag, package versions, `random_seeds`
- run metadata (excluded from ID): `started_utc`, `finished_utc`,
  `output_files` with checksums, `warning_summary`

### Constraint — `con-<hash12>`

Injection-calibrated sensitivity or a qualified null, always traceable
to an analysis run and usable pixels — never to a bare footprint hit.

- `constraint_id`
- refs: `analysis_run_id`, `endpoint_id`, `role`
- cell: `z_interval_au`, `band`, `epoch_range_mjd`, `duty_cycle_range`,
  `residual_motion_bound_arcsec_per_yr`, `morphology`,
  `hypothesis_version`
- result: `kind` (`recovery_curve` | `qualified_null` |
  `not_constrainable`), `recovery_probability`, `flux_limit` (+ unit,
  band-specific physical interpretation label per hypotheses §6),
  `injection_summary_ref`, `false_alarm_rate`, `trials_accounting_ref`

### Candidate — `cnd-<hash12>`

A retained event or track with its competing-model tests.

- `candidate_id`
- refs: `analysis_run_id`, `endpoint_id`, `role`,
  contributing `observation_ids`
- track: fitted `z_au` (+ interval), fitted residual motion, per-epoch
  measurements ref
- model comparison: scores/likelihoods for SGL track vs.
  inertial-background vs. stellar parallax+PM vs. Keplerian
  solar-system vs. instrumental (plan §3.4)
- vetting: `status` (`new` | `vetoed` | `pending_holdout` |
  `confirmed_prediction_failed` | `retained`), `holdout_test_ref`,
  `veto_reason` (e.g. known-object match via `neowiser_p1ba_mch`)

## Query snapshots — `snp-<hash12>`

Not one of the five, but the audit substrate: `{service_url, exact
query, request_utc, response_sha256, response_path (by convention, not
in ID), row_count, http_status}`. Response bytes stored verbatim under
`runs/<archive>/<run-id>/snapshots/`.

## Adapter interface

```python
class FootprintTest(Protocol):
    def contains(self, ra_deg: float, dec_deg: float) -> bool: ...
    # passed directly as covered_z_intervals(contains=...)

class ArchiveAdapter(Protocol):
    archive_id: str

    def discover(self, region: DiscoveryRegion, time_range: MjdRange,
                 store: SnapshotStore) -> Iterator[Observation]:
        """Query archive metadata; snapshot raw responses; yield
        normalized Observation records."""

    def nominal_footprint(self, obs: Observation) -> FootprintTest:
        """Cheap test from inventory metadata (corner polygon / nominal
        WCS). Suitable for the coarse pass only."""

    def exact_footprint(self, obs: Observation,
                        products: ProductSet) -> FootprintTest:
        """Exact spherical WCS + mask + quality cuts; requires
        downloaded products. Defines 'usable' pixels (plan §3.3)."""

    def fetch(self, obs: Observation, kinds: Sequence[str],
              dest: Path, *, cutout: CutoutSpec | None = None
              ) -> ProductSet:
        """Download products (full or cutout), verify checksums,
        record local paths + hashes."""
```

`DiscoveryRegion` is produced by the geometry layer (`sglsurvey.geometry`)
from sglseti loci: adaptive/swept loci per endpoint × role × epoch chunk,
MC-propagated to the frozen 99% confidence radius, padded +10", then
converted to the adapter's preferred query form (cone/polygon per epoch
chunk now; MOC later, plan §7 "footprint representation"). The geometry
layer records every sglseti identity (`model_id`/`version`,
`ephemeris_id`, provider triples, `target_source_hash`, seeds) for
inclusion in IntersectionEvaluation records.

## WISE adapter specifics

- Discovery: ADQL against `wise.neowiser_merge_p1bm_frm` (verified union
  of all four mission phases) via TAP; `/TAP/async` for production.
- Products: IBE per-frame template (verified — see
  `surveys/wise/notes/irsa_recon.md`), `.md5` sidecars for checksum
  verification, `?center=&size=` cutouts to stay inside disk budget.
- Exact footprint: WCS from `int` header (distortion terms included),
  `msk` bitmask + `fflag`/quality columns define usable pixels; the
  precise pass evaluates the swept locus over `t_mid ± exptime/2`.

## Open items

- Time-resolved unWISE hosting and reject-table names per phase
  (recon TO VERIFY).
- Injection framework design (separate doc once the forced-photometry
  stage exists).
- MOC-based discovery representation once pilot-scale cones prove
  limiting.
