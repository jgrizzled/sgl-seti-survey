# Rubin alert-broker reachability recon — 2026-08-27

Scope: which community brokers can serve **historical** and **live**
LSST alert data for the plan §6 near-term alert project (crossings
screening + corridor track-screen), probed from the dev machine.
Brokers examined: Lasair-LSST, Fink-LSST, Babamul (BOOM/Caltech).
All service facts below verified live today unless marked _(docs)_.

## 0. Programme-level finding: the stream is DOWN

- Both Fink and Babamul nightly statistics show **LSST alert
  production stopped 2026-07-13/14**. Babamul: last non-zero night
  2026-07-13 (473k alerts); zero every night through 2026-08-26.
  Fink: last ingested night 2026-07-14.
- Cause: the record mid-July 2026 winter storm in the Coquimbo
  region; the summit was evacuated 2026-07-14 and recovery is slow
  (snow-blocked summit road, follow-on storms) —
  [community.lsst.org status thread](https://community.lsst.org/t/rubin-observatory-status/12397),
  [LSST Update 2026-07-24](https://community.lsst.org/t/lsst-update-2026-07-24/12332).
- Consequence: there is **no live stream to consume right now**; the
  near-term substrate is the brokers' historical archives
  (science-validation era 2025-10 → operations start 2026-06-30 →
  2026-07-13). Standing watch infrastructure can be built and tested
  against the historical API and switched on when the stream resumes.

Archive span (Fink `/statistics`, anonymous): **132 nights,
2025-10-25 → 2026-07-14, 16,379,145 alerts** — i.e. the broker
archives reach back into commissioning, ~5 months before the DP2
catalog era ends (2025-12-08), and cover ~7.5 months (2025-12-08 →
2026-07-13) accessible **nowhere else** until the Year-1 release.

## 1. Fink-LSST — verified end-to-end, anonymous

Base `https://api.lsst.fink-portal.org` (Fink/LSST object API
v3.7.0; swagger at `/swagger.json`). **No authentication** for the
REST API; all probes below succeeded anonymously.

| Endpoint | Verified | Notes |
| --- | --- | --- |
| `POST /api/v1/conesearch` | ✅ | ra/dec/radius (max 5°!), optional first-detection date window (`startdate`/`stopdate`/`window`, `kind` first/last semantics), column selection, json/csv/parquet/votable. Returns per-alert rows incl. `r:reliability` (real/bogus) and `v:separation_degree`. One query returned alerts from MJD 61095 and 61235 → full historical archive by position. |
| `POST /api/v1/fp` | ✅ | **Forced photometry** by `diaObjectId` (comma-list OK). Probe returned 13 epochs spanning MJD 61219→61235 for an object first alerting 61235 — i.e. the precovery `ForcedSourceOnDiaObject` window (~2 weeks pre-first-alert), not arbitrary positions. |
| `POST /api/v1/cutouts` | ✅ | Per `diaSourceId`, kind Science/Template/Difference/All, **FITS or PNG or array**. Probe fetched a 40 KB FITS difference stamp. |
| `POST /api/v1/sources` | spec | Per-source lightcurve rows by `diaObjectId`. |
| `POST /api/v1/objects` | spec | Aggregated object record. |
| `POST /api/v1/tags` | ✅ | Alerts by Fink science-filter tag with date range. |
| `GET /api/v1/statistics` | ✅ | Nightly ingest stats (used for the span above); per-night flag counts incl. `pixelFlags_*`, `is_sso`, `isDipole`. |
| `/api/v1/skymap`, `/sso`, `/ssoft`, `/resolver`, `/schema`, `/blocks` | spec | SSO endpoints useful for positive controls. |

Caveats: `conesearch` time filters select on an object's **first/last
detection date**, not per-epoch — for in-window screening, query
without date filter (or wide) and cut on `r:midpointMjdTai`
client-side. No published rate limit; be polite. Live Kafka
(fink-client substreams) exists but needs a registration form —
irrelevant while the stream is down; revisit at resumption.

## 2. Lasair-LSST — VERIFIED with account token (2026-08-27)

Base `https://api.lasair.lsst.ac.uk/api/…`; token in `.env` as
`LASAIR_API_TOKEN` (header `Authorization: Token …`). Free account
(email only); **no public demo token** in the LSST docs (unlike
Lasair-ZTF). Limits: 100 calls/hr, 10k rows (standard); 10k
calls/hr, 1M rows (power user on request to
lasair-help@mlist.is.ed.ac.uk).

Verified today at the Fink test position (64.18537, −48.94322):

- `/api/cone/` ✅ — same 7 diaObjects as Fink/Babamul (cross-broker
  ID consistency confirmed).
- `/api/object/` ✅ — record contains `diaSourcesList` (16 rows,
  full per-source fields incl. dipole, detector, flags) **and
  `diaForcedSourcesList`** (13 rows with `psfFlux` *and*
  `scienceFlux` direct-image flux + visit/detector IDs — richer than
  Fink `/fp`, which returns difference fluxes only). No
  non-detection list.
- `/api/query/` (arbitrary SQL over `objects` incl. derived
  features — jump detector, Bazin fits, per-band stats) ✅ — note
  there is **no `nSources` column** (error surfaced the generated
  SQL; use the schema browser for real column names, e.g.
  `latest_psfFlux`, `lastDiaSourceMjdTai` both verified).
- `/api/sherlock/position/` ✅ — arbitrary-position sky context
  (test object classified SN on a 2MASS/DESI galaxy 0.1″ away) —
  useful blend/host annotation, cf. the ross-128 S1 blend lesson.

Unique capabilities for us:

- **Watchmaps: user-uploaded MOC regions** — every ingested alert is
  tagged against the MOC and filterable; **past alerts can be matched
  retroactively on request** (quote the watchmap URL to the helpdesk).
  A corridor MOC (Pipeline A) or crossing-cone MOC set (Pipeline B)
  becomes a server-side standing filter.
- **Watchlists** (position lists with radii) + SQL filters → **Kafka
  topics** (`lasair-lsst-kafka_pub.lsst.ac.uk:9092`, 7-day retention,
  resumable by group_id) or daily e-mail — the prospective-watch
  mechanism, no polling loop needed.
- Python client: `pip install lasair` (v0.1.4).

## 3. Babamul (BOOM, Caltech) — VERIFIED incl. Kafka (2026-08-27)

Base `https://api.kaboom.caltech.edu/babamul` (OpenAPI 3.1 spec
embedded in `/babamul/docs`). Serves **both `lsst` and `ztf`**
through one interface. Auth: JWT in `.env` as `BABAMUL_API_TOKEN`
(header `Authorization: Bearer …`); Kafka pair as
`BABAMUL_KAFKA_USER`/`BABAMUL_KAFKA_PW`. Open endpoints (no auth):
`/stats/nightly` (used for the outage finding), `/stats/kafka`,
`/surveys/{survey}/schemas`.

Verified today:

- `POST /surveys/lsst/objects/cone-search` ✅ and
  `POST /surveys/lsst/alerts/cone-search` ✅ with **JD window**
  (`start_jd`/`end_jd`) + quality filters (`min/max_drb`,
  `min/max_magpsf`, `is_rock`/`is_star`/`is_near_brightstar`/
  `is_stationary`/`is_positive`). **Gotcha: `coordinates` is a
  Kowalski-style map `{"name": [ra, dec]}`** — not `{ra:, dec:}` —
  and is natively **batch** (tested 2 positions in one call; the 88
  crossing positions can go up in one request). Returned the same
  objects as Fink/Lasair; 8 alerts in a 10-day window at the test
  position. Each alert's `candidate` carries **`diffmaglim`
  (per-visit difference-image limiting mag)**, `drb`-class scores,
  forced-flux flags, `glint_trail`, etc.
- `GET /surveys/lsst/cutouts` ✅ — **`which` selects the alert, not
  the image**: enum `First/Last/Brightest/Faintest` (per band via
  `band=`); the response is one JSON with all three stamps
  (`cutoutScience/Template/Difference`, base64). 177 KB for the
  test candid.
- **Kafka live stream ✅** — `pip install babamul` (v0.1.0),
  `AlertConsumer(topics=[…], offset="earliest", group_id=…)`;
  consumed 2 real alerts from `babamul.ztf.lsst-match.hosted`.
  Client env names differ from ours: map
  `BABAMUL_KAFKA_USERNAME=$BABAMUL_KAFKA_USER`,
  `BABAMUL_KAFKA_PASSWORD=$BABAMUL_KAFKA_PW`. Topic census
  (`/stats/kafka`, 14 topics, 7-day retention): **all 8
  `babamul.lsst.*` topics 0 alerts** (stream down) while
  `babamul.ztf.*` topics are live — ZTF flows validate the plumbing
  until LSST resumes. Default consumer offset is `latest` (blocks
  forever on a quiet topic — use `timeout`).

## 4. Gaps common to all three

- **No coverage/non-detection service at arbitrary positions.**
  Forced photometry exists only at DiaObject positions (Fink `/fp`,
  Lasair `diaForcedSourcesList`); none of the three serves per-visit
  magLim or visit footprints for the 2026 era, so "window was
  observed but produced no alert" cannot be established from brokers
  alone. **Partial route found:** Babamul alert candidates carry
  `diffmaglim` — any alert from the *same visit* near a window
  position (a wider cone) yields that visit's difference-image depth
  and proves the visit happened; alert-poor visits/fields remain
  invisible. Full closure still needs Rubin visit/ConsDB metadata,
  scheduler visit logs, or helpdesk asks. This is the blocking item
  for a coverage ledger — same class as the §5.13 magLim-referenced
  threshold caveat.
- **No archive injections** in the alert stream (DP2's
  `pixelFlags_injected` has a Fink per-night counter, so injected
  sources are flagged if present) → anything built here stays
  threshold/screening class per the C1 injection-calibration rule.
- Alert history in packets is ~12 months of prior sources + ~2 weeks
  of precovery forced photometry — fine for vetting, not a substitute
  for the Year-1 image-level survey.

## 5. Recommended roles

- **Fink** — the workhorse for historical screening now: anonymous,
  5° cones, full archive by position, forced photometry, FITS stamps.
  Adapter cost lowest of the three.
- **Lasair** — the standing-watch instrument (watchmap MOCs +
  filters → Kafka/email) and SQL/Sherlock annotation source; needs
  one free account + token (user action).
- **Babamul** — best per-epoch filtered cone-search (JD windows +
  drb) and a second opinion/cross-check on Fink hits, plus unified
  ZTF+LSST access; needs email signup (user action).

Non-git credentials (`.env`, all present and verified 2026-08-27):
`LASAIR_API_TOKEN`, `BABAMUL_API_TOKEN`, `BABAMUL_KAFKA_USER`,
`BABAMUL_KAFKA_PW`. Fink needs none. All three brokers returned
identical diaObject IDs at the test position — cross-broker identity
is safe to rely on.
