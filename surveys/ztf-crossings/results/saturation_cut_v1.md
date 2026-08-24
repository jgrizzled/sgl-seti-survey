# Channel-A saturation cut v1 — results

Run 2026-08-23 by `scripts/saturation_cut.py` over the 86 channel-A era
targets. Frozen rule (per band, appended to hypotheses freeze v1.0):
**excluded** m_est < 13.0 (ZTF ~12.5 mag saturation + 0.5 mag margin
absorbing the ±0.3 mag color-transform accuracy), **marginal**
13.0–13.5 (searchable, flagged), **ok** ≥ 13.5. Estimates from Gaia DR3
G and G−RP via a piecewise dwarf-sequence table (CNS5 census snapshot;
CNS5-id fallback); seven Gaia-saturated classics plus EZ Aqr
(V≈12.2, M5V colors) and GJ 11068 (2MASS J=10.63 at 6.8 pc, ~M6.5
colors) from a sourced hand table. Late-T/Y objects (Luhman 16 A/B,
WISE 0855) are `optically_dark`: they cannot saturate and their
channel-A sensitivity statement is the empty-field limit, not a
contrast limit.

| band | excluded | marginal | ok | searchable targets | surviving covered event rows |
|---|---|---|---|---|---|
| g | 59 | 5 | 22 | 27 | 136 |
| r | 66 | 2 | 18 | 20 | 90 |
| i | 71 | 0 | 14 | 14 | 5 |

Searchable = ok + marginal. Full per-target table:
`saturation_cut_v1.ecsv`; machine summary: `saturation_cut_v1.json`.

Notes:
- The blended channel-A search population is entirely M dwarfs, white
  dwarfs (eps-ind-b is the T dwarf companion — bright-limit safe), and
  the optically-dark objects; every F/G/K endpoint and the bright
  classics are excluded in all bands, as expected at ≤ 10 pc.
- Teegarden's star is g-ok / r-marginal — relevant because it is also
  one of the two channel-B grazing-radius targets; channel B is an
  empty-field search at the antipode, so this cut does not apply there.
- The chromatic discriminator needs both bands: 20 targets are
  searchable in g **and** r; 7 more are g-only (their r-band excess
  test degrades to g-band-only phase-locking).
- i-band survives on only 5 covered event rows; it stays opportunistic
  as frozen.
