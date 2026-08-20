# targets — universal (survey-agnostic) target portfolio

The target list all surveys draw from. Membership and identity are
decided here from the network prior + host flags + universal
searchability; each survey applies a thin overlay
(`surveys/<archive>/targets/`) that grades, orders, and possibly
excludes corridors for that survey — overlays never add or redefine
targets.

- `config.yaml` — all knobs (horizon, basket sizes, vote thresholds,
  desirability thresholds, engineering weights, orbit-availability
  map, track budget)
- `scripts/pull_census.py` — snapshots CNS5 within the horizon and the
  engineering joins (TIC, Gaia DR3 main + astrophysical params,
  Kervella+2022 PMa, HGCA, Boro Saikia+2018 log R'HK, CARMENES Hα,
  ROSAT 2RXS) into `census/`
- `scripts/build_universal_list.py` — the builder (Stages 0–4 of the
  simplified ranking method in `notes/sgl_seti_star_ranking_methods.md`
  plus the Stage 2b engineering/desirability layer from
  `notes/picky_network_hypothesis_sgl_seti.md`)
- `science_interest.yaml` — curated basket D (tier 1 auto-included,
  tier 2 budget-permitting), one line of rationale + reference each
- `census/` — dated catalog snapshots (never edited in place)
- `universal_v2.json` / `universal_v2.md` — the current versioned
  output; `universal_v1.*` is the 8 pc list the WISE survey ran on

## Baskets (v2)

| Basket | Question it answers | Source |
| --- | --- | --- |
| Distance Core | nearest N systems, no exclusions | distance rank |
| Robust Geometric Neighbor | natural Sun-neighbor among *all* systems | ≥3 of 7 sparse-graph families |
| Selective Network Neighbor | Sun-neighbor once undesirable hosts are pruned | same vote on desirability level 1 and level 2 subgraphs |
| Engineering Backbone | best relay hosts regardless of geometry | rank-sum of M/R², power, quietness, cleanliness, lifetime, distance among level-2 dwarfs |
| Science Interest | systems a network might join for the science | curated YAML |
| Historical Neighbor | past/future close approach | linear motion ±1 Myr |
| Compact Lens Wildcard | compact-remnant architectures | isolated white dwarfs |

*Desirability* (the network's presumed taste: multiplicity, detected
acceleration, activity, evolutionary state) and the *solution gate*
(our ability to predict the track) are separate columns and must stay
that way — v1's "selective" basket conflated them.

Design rules: the distance core is protected (no exclusions, ever);
network-prior, engineering, and archive-searchability scores stay
separate columns; every member carries the basket labels that justify
its inclusion; unknown ≠ favorable (activity/acceleration classes carry
their evidence strings); targets failing a solution gate are *deferred
with a named unblocking condition*, not dropped. Feeds `registries/`
curation in priority order.
