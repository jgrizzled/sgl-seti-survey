# Rubin DP2 crossings dev stage (2026-08-26; pseudo-units only, D8)

`scripts/dev_stage.py` → `results/dev_v1.json` (v1.0-rule run
preserved as `dev_v1_rule10as_superseded.json`); snapshots
`runs/rubin-crossings/dev_v1/`. A `BlindGuardClient` asserted the
confirmatory visit id never appeared in any dev query.

## Pseudo-unit family and null

63 visits within 2° of the ross-128 axis point in the DP2 era; 1
in-window (excluded, blind), 62 off-window pseudo-units run through
the identical frozen chain (z-family at each epoch, footprint gate,
static exclusion, S_det, 8 controls):

- 53 searched; 9 off-detector (dither/edge epochs).
- **Association rate S_det > 0: 0/53. Exceedances: 0.** The 60″-cone
  DiaSource counts are 0–1 per visit — the sparse-substrate regime;
  the 1/9 exchangeable control-crossing budget is very conservative
  here.
- Static objects within 90″ of the axis point: 21; per-epoch live
  z-positions 3–5 of 5.

## Amendment v1.1 (found here)

The v1.0 control rule (≥ 10″ locus-avoidance) is geometrically
impossible for inner-z offsets (~2–5″ from the axis point): several
pseudo-units — and by the same geometry the confirmatory unit — had
zero valid control patterns. Amended to ≥ 2.5″ (= 2.5 × r_assoc;
non-overlap of the 1″ association apertures with margin);
hypotheses.md v1.1; freeze re-bound
(`sha256:9384e97a…0b44da4`). Under v1.1: 161/424 control patterns
needed extra rotation, all resolved; null unchanged (0/53).

## Positive control — PASS

SkyBoT cone (obs code X05, TAI→UTC converted epoch) at the position
of an `ssObjectId`-linked DiaSource on off-window visit
2025071700513 (r): predicted position of **asteroid 2006 SE393**
0.308″ from the catalog DiaSource; the frozen association chain
recovered it at **S_det = 31.3** and the associated `diaSourceId`
equals the catalog's own solar-system link. Two prior attempts where
the 3′ SkyBoT cone contained only an unrelated object 144″ away are
recorded as non-matches (`positive_control.attempts`), not chain
failures.

## Archive-injection census — completeness route (ii)

0 `pixelFlags_injected*` DiaSources in 11 off-window r visits at the
field: DP2 (Early) provides no usable injection population here. The
confirmatory constraint is therefore a **magLim-referenced threshold
statement, labeled not injection-calibrated** (freeze §8, C1 rule —
no exclusion claim is made on it).

All dev gates pass → blind confirmatory authorized under
v1.0 + v1.1.
