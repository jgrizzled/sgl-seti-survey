# targets — universal (survey-agnostic) target portfolio

The target list all surveys draw from. Membership and identity are
decided here from the network prior + host flags + universal
searchability; each survey applies a thin overlay
(`surveys/<archive>/targets/`) that grades, orders, and possibly
excludes corridors for that survey — overlays never add or redefine
targets.

- `config.yaml` — all knobs (horizon, basket sizes, vote thresholds,
  orbit-availability map, track budget)
- `scripts/build_universal_list.py` — the builder (Stages 0–4 of the
  simplified ranking method in `notes/sgl_seti_star_ranking_methods.md`)
- `census/` — snapshotted catalog pulls (CNS5 via VizieR)
- `universal_v1.json` / `universal_v1.md` — the versioned output:
  systems with basket labels, adjacency votes, host flags, solution
  gates, and universal searchability columns

Design rules: the distance core is protected (no exclusions, ever);
network-prior scores and archive searchability stay separate columns;
every member carries the basket labels that justify its inclusion;
targets failing a solution gate are *deferred with a named unblocking
condition*, not dropped. Feeds `registries/` curation in priority
order.
