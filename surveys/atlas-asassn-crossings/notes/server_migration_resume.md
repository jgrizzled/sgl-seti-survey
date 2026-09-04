---
title: "ATLAS/ASAS-SN crossings — server migration & resume checklist"
date: 2026-08-25
status: "RESUMED 2026-09-03 on the new server; survey COMPLETE 2026-09-04 (report/atlas_asassn_crossings.md)"
---

# Resume checklist (new machine)

**Resolved 2026-09-03.** All four "needs" verified on the server
(repo + `../sglseti` 19c8167 clean, `.env` token, `uv` env Python
3.12.2, ATLAS queue 200 / ASAS-SN port 9006 open); `runs/` is a
symlink to the 8 TB SSD (`/mnt/ssd2/astro/sgl-seti-survey-runs`).
First actions 1–4 done the same day (era end MJD 61286; ASAS-SN
ceiling unchanged at 2025-06-16; no radec-list batching exists —
schema snapshot under `runs/atlas-asassn-crossings/docs_snapshot/`;
windowed task = ~2.2 min; `scripts/atlas_drain.py` running the
308-task fleet, complete 2026-09-04, 0 failures). Item 5 done (`scripts/asassn_ledger.py`); item 6 continues from the threshold freeze.

State at pause (2026-08-25): recon complete
(`atlas_asassn_recon_2026-08-25.md`), era scoping done
(`../results/era_scope_v0.json`), **hypotheses v1.0 frozen** (D1–D8
adopted as recommended). Nothing after the freeze has run: no
coverage stage, no survey-position light curves in the repo. The
recon's two probe light curves lived only in a session scratchpad on
the macbook and are intentionally not carried over (D1: everything is
re-pulled fresh under snapshot discipline).

## What the new machine needs

1. **Repo checkouts**: this repo + `../sglseti` adjacent (editable
   path dep in `pyproject.toml`; pin/record its commit as usual).
2. **Secrets — the one thing git does not carry**: create `.env` in
   the repo root with `ATLAS_API_TOKEN=<token>` (gitignored; copy
   from the macbook or regenerate at
   https://fallingstar-data.com/forcedphot — the token page allows
   viewing/replacing). No other credentials exist in this survey
   (ASAS-SN Sky Patrol v2 is anonymous).
3. **Environment**: `uv sync` (Python ≥ 3.12.2). The ASAS-SN client
   is pinned in pyproject (`skypatrol==0.6.21`, import name
   `pyasassn`); `astro-datalab==2.22.1` and `pandas` are already
   pinned from DECam.
4. **Network reachability smoke test** (some of these are
   institution-network sensitive; ASAS-SN v2 is plain HTTP on port
   9006 — confirm the server's egress policy allows it):

   ```
   .venv/bin/python -c "
   import sys; sys.path.insert(0, 'surveys/atlas-asassn-crossings/scripts')
   import atlas_api, requests
   h = atlas_api.headers()
   r = requests.get('https://fallingstar-data.com/forcedphot/queue/', headers=h)
   print('ATLAS queue:', r.status_code)   # expect 200
   from pyasassn.client import SkyPatrolClient
   c = SkyPatrolClient(verbose=False)
   print('SkyPatrol catalogs OK:', len(c.catalogs.catalog_names()) if hasattr(c.catalogs,'catalog_names') else 'client init OK')
   "
   ```

   Plus the multiprocessing gotcha check: any script using
   `download=True` must run as a file under
   `if __name__ == '__main__':` (spawn-loops from stdin).

## First actions on resume (coverage stage)

1. Record the **ATLAS era end date** = resume date (hypotheses §2)
   and **re-measure the ASAS-SN v2 DB ceiling** (was JD 2460841 ≈
   2025-06-15 at recon; a reprocessing banner was up — the ceiling
   may have moved).
2. Re-run `scripts/era_scope.py` with the recorded era end (extend
   `ERAS` accordingly) — event counts grow with the era.
3. Probe **`radec`-list batching** (open item from the recon: batch
   size limit; whether list tasks parallelize within one account)
   before committing to the ~1,170-task per-event pull pattern —
   execution measured serial at ~25–70 min per *full-history* task;
   windowed per-event tasks (±110 d) should be much faster, but
   measure one before launching the fleet.
4. Build the coverage-stage scripts on `scripts/atlas_api.py`
   (submit → sidecar-snapshotted fetch under `runs/atlas-asassn-
   crossings/`), per-event windowed tasks per hypotheses §6,
   **wolf-359 `forced_dev`** (D1) wired into the split logic from the
   start.
5. ASAS-SN coverage-fraction ledger pull (v2 epoch lists at
   catalogued sources adjacent to each survey position) — coverage
   records only, per hypotheses §1.
6. Then the usual chain: coverage → threshold freeze (seed + dev/
   confirmatory split per D8) → dev → confirmatory → completeness →
   report.

## Ops notes for a long-running server pull

- The ATLAS queue is serial per account: the production pull is a
  days-long background drain. `callback_url` (must be public HTTPS,
  fired once, no retry) is available if the server has a public
  endpoint; otherwise slow polling (≥ 30 s) is fine.
- Result URLs are short-lived: fetch-then-DELETE immediately
  (already the `atlas_api.wait_and_fetch` behavior).
- Tasks can also be inspected at
  `fallingstar-data.com/forcedphot/queue/` (web) and
  `queuepositions.json` (API) if a pull needs triage mid-drain.
- ATLAS had a multi-day outage recovered 2026-08-21 — build the
  drain to be resumable (sidecar files double as the done-list).
