# Dev stage v1 — 2026-09-07

T = 3.451 (145 trials, FWER 0.05).

- dev unit **A gj-1276 1.2Rsun**: S_stack_L 0.04, S_event_L 0.80, S_stack_H 0.57, S_event_H 1.68, S_burst 2.08 (15 searchable windows; exceed T: [])
- dev unit **A gj-908 0.1AU**: S_stack_L 0.11, S_event_L 3.16, S_stack_H 2.57, S_event_H 2.07, S_burst 2.10 (18 searchable windows; exceed T: [])
- dev unit **B gj-908 0.1AU**: S_stack_L 0.36, S_event_L 1.58, S_stack_H 0.42, S_event_H 0.86, S_burst 1.58 (18 searchable windows; exceed T: [])
- D1 window B van-maanen 2020-04-02 1.2Rsun: S_stack_L 0.23, S_event_L 0.23, S_stack_H 0.00, S_event_H 0.00, S_burst 0.40; window n_L 4 vs λ 4.11, n_H 0 vs λ 0.032
- D1 window B van-maanen 2020-04-02 2.5Rsun: S_stack_L 0.23, S_event_L 0.23, S_stack_H 0.00, S_event_H 0.00, S_burst 0.24; window n_L 12 vs λ 12.41, n_H 0 vs λ 0.082
- D1 window B van-maanen 2020-04-02 0.1AU: S_stack_L 0.15, S_event_L 0.15, S_stack_H 0.00, S_event_H 0.00, S_burst 0.11; window n_L 105 vs λ 110.20, n_H 0 vs λ 1.442
- control **ctrl-grb130427a**: S_stack_L 300.00, S_event_L 300.00, S_stack_H 206.38, S_event_H 206.38, S_burst 300.00; expected ['S_event_L', 'S_event_H', 'S_burst'] > T → **PASS**; window n_L 240 vs λ 4.74
- control **ctrl-3c454.3**: S_stack_L 300.00, S_event_L 300.00, S_stack_H 300.00, S_event_H 300.00, S_burst 300.00; expected ['S_event_L'] > T → **PASS**; window n_L 7353 vs λ 816.26

## Machinery checks

- Gate census (raw → SOURCE+zenith → interval-gated): A-gj-1276 142165/70180/49056; A-gj-908 137607/70279/52728; A-ross-128 192508/71717/52363; A-ross-154 655527/338314/269148; A-teegarden 386523/127131/96356; A-van-maanen 166119/73562/53177; A-wolf-359 175278/78122/54833; B-gj-1276 176175/79188/55749; B-gj-908 193746/71481/53596; B-ross-128 133869/69353/50774; B-ross-154 563858/162702/122761; B-teegarden 773037/129330/98841; B-van-maanen 470216/210919/166963; B-wolf-359 142542/71180/49137
- Recon-file identity: {'ours_gated_in_span': 70, 'found_in_recon_file': 70, 'recon_file_photons': 150}
- Moon/Sun exclusion: 7981.1 ks excluded of 33600 ks gated over 137/580 windows
- IRF: L A_eff on-axis 0.517 m², containment 0.61 (cosθ 0.8) / 0.49 (0.5), ⟨E⟩ 256 MeV; H A_eff on-axis 0.894 m², containment 0.92 (cosθ 0.8) / 0.84 (0.5), ⟨E⟩ 5723 MeV; U A_eff on-axis 0.555 m², containment 0.40 (cosθ 0.8) / 0.30 (0.5), ⟨E⟩ 801 MeV