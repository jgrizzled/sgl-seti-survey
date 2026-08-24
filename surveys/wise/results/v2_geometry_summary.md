---
title: "WISE v2 — geometry: covariance envelopes and independent check"
date: 2026-08-22
---

# Covariance propagation (hypotheses v2.0 §1.4)

sglseti seeded Monte Carlo, N = 2000, seed 20260822, at z = 550 and 10,000 AU and epochs (2010.0, 2017.0, 2024.0); observer `wise-l1b-spacecraft`. Cross-track decision threshold 0.5 x FWHM(W1) = 3.05".

**3 of 176 endpoint-role cells exceed the threshold** (max sigma_xt_99 over all cells: 16.7831").

| endpoint | role | provider | max sigma_xt_99 ["] | max sigma_at_99 ["] | max r_99 ["] | cross-track dim |
| --- | --- | --- | --- | --- | --- | --- |
| 61-cyg-a | rx | linear_astrometry_v1 | 0.0012 | 0.0011 | 0.0013 | no |
| 61-cyg-a | tx | linear_astrometry_v1 | 0.0205 | 0.0653 | 0.0689 | no |
| 61-cyg-b | rx | linear_astrometry_v1 | 0.0006 | 0.0005 | 0.0006 | no |
| 61-cyg-b | tx | linear_astrometry_v1 | 0.0099 | 0.0309 | 0.0325 | no |
| 61-vir | rx | linear_astrometry_v1 | 0.0027 | 0.0031 | 0.0035 | no |
| 61-vir | tx | linear_astrometry_v1 | 0.2518 | 0.1212 | 0.2791 | no |
| 82-eri | rx | linear_astrometry_v1 | 0.0059 | 0.0053 | 0.0070 | no |
| 82-eri | tx | linear_astrometry_v1 | 0.2294 | 0.3917 | 0.4536 | no |
| alpha-cen-a | rx | two_body_orbit_v1 | 0.1688 | 0.1762 | 0.1980 | no |
| alpha-cen-a | tx | two_body_orbit_v1 | 0.2363 | 0.2635 | 0.2981 | no |
| alpha-cen-b | rx | two_body_orbit_v1 | 0.1657 | 0.1727 | 0.2025 | no |
| alpha-cen-b | tx | two_body_orbit_v1 | 0.2460 | 0.2635 | 0.3026 | no |
| barnard-star | rx | linear_astrometry_v1 | 0.0006 | 0.0007 | 0.0008 | no |
| barnard-star | tx | linear_astrometry_v1 | 0.0076 | 0.0237 | 0.0248 | no |
| eps-eri | rx | linear_astrometry_v1 | 0.0083 | 0.0096 | 0.0108 | no |
| eps-eri | tx | linear_astrometry_v1 | 0.0526 | 0.0665 | 0.0794 | no |
| eps-ind-a | rx | linear_astrometry_v1 | 0.0018 | 0.0018 | 0.0021 | no |
| eps-ind-a | tx | linear_astrometry_v1 | 0.0205 | 0.1044 | 0.1054 | no |
| eps-ind-b | rx | linear_astrometry_v1 | 4.2297 | 4.0117 | 4.9970 | YES |
| eps-ind-b | tx | linear_astrometry_v1 | 16.7831 | 16.0079 | 19.9009 | YES |
| ez-aqr | rx | linear_astrometry_v1 | 2.1761 | 2.2482 | 2.6399 | no |
| ez-aqr | tx | linear_astrometry_v1 | 7.7581 | 7.6390 | 9.2323 | YES |
| fomalhaut | rx | linear_astrometry_v1 | 0.0319 | 0.0423 | 0.0453 | no |
| fomalhaut | tx | linear_astrometry_v1 | 0.1285 | 0.1799 | 0.2057 | no |
| gj-1002 | rx | linear_astrometry_v1 | 0.0008 | 0.0013 | 0.0014 | no |
| gj-1002 | tx | linear_astrometry_v1 | 0.0275 | 0.0271 | 0.0388 | no |
| gj-1061 | rx | linear_astrometry_v1 | 0.0008 | 0.0007 | 0.0009 | no |
| gj-1061 | tx | linear_astrometry_v1 | 0.0067 | 0.0034 | 0.0071 | no |
| gj-1087 | rx | linear_astrometry_v1 | 0.0004 | 0.0004 | 0.0005 | no |
| gj-1087 | tx | linear_astrometry_v1 | 0.0038 | 0.0218 | 0.0221 | no |
| gj-11068 | rx | linear_astrometry_v1 | 0.3364 | 0.3214 | 0.3973 | no |
| gj-11068 | tx | linear_astrometry_v1 | 1.3549 | 1.3502 | 1.6340 | no |
| gj-1111 | rx | linear_astrometry_v1 | 0.0010 | 0.0014 | 0.0015 | no |
| gj-1111 | tx | linear_astrometry_v1 | 0.0084 | 0.0158 | 0.0175 | no |
| gj-11547 | rx | linear_astrometry_v1 | 0.0011 | 0.0011 | 0.0014 | no |
| gj-11547 | tx | linear_astrometry_v1 | 0.0223 | 0.0073 | 0.0227 | no |
| gj-1221 | rx | linear_astrometry_v1 | 0.0022 | 0.0030 | 0.0038 | no |
| gj-1221 | tx | linear_astrometry_v1 | 0.0687 | 0.0957 | 0.1178 | no |
| gj-12724 | rx | linear_astrometry_v1 | 0.0019 | 0.0021 | 0.0024 | no |
| gj-12724 | tx | linear_astrometry_v1 | 0.0991 | 0.0238 | 0.1015 | no |
| gj-1276 | rx | linear_astrometry_v1 | 0.0018 | 0.0022 | 0.0028 | no |
| gj-1276 | tx | linear_astrometry_v1 | 0.1147 | 0.1483 | 0.1885 | no |
| gj-13157 | rx | linear_astrometry_v1 | 0.0004 | 0.0004 | 0.0004 | no |
| gj-13157 | tx | linear_astrometry_v1 | 0.0042 | 0.0053 | 0.0060 | no |
| gj-2012 | rx | linear_astrometry_v1 | 0.0004 | 0.0007 | 0.0007 | no |
| gj-2012 | tx | linear_astrometry_v1 | 0.0128 | 0.0370 | 0.0388 | no |
| gj-2066 | rx | linear_astrometry_v1 | 0.0005 | 0.0004 | 0.0006 | no |
| gj-2066 | tx | linear_astrometry_v1 | 0.0075 | 0.0098 | 0.0119 | no |
| gj-229-a | rx | linear_astrometry_v1 | 0.0002 | 0.0004 | 0.0004 | no |
| gj-229-a | tx | linear_astrometry_v1 | 0.0013 | 0.0069 | 0.0069 | no |
| gj-251 | rx | linear_astrometry_v1 | 0.0006 | 0.0005 | 0.0007 | no |
| gj-251 | tx | linear_astrometry_v1 | 0.0093 | 0.0083 | 0.0121 | no |
| gj-293 | rx | linear_astrometry_v1 | 0.0011 | 0.0019 | 0.0022 | no |
| gj-293 | tx | linear_astrometry_v1 | 0.0580 | 0.1056 | 0.1206 | no |
| gj-3112 | rx | linear_astrometry_v1 | 0.0008 | 0.0007 | 0.0010 | no |
| gj-3112 | tx | linear_astrometry_v1 | 0.0548 | 0.0459 | 0.0712 | no |
| gj-318 | rx | linear_astrometry_v1 | 0.0004 | 0.0017 | 0.0017 | no |
| gj-318 | tx | linear_astrometry_v1 | 0.0142 | 0.0995 | 0.1006 | no |
| gj-3306 | rx | linear_astrometry_v1 | 0.0012 | 0.0009 | 0.0015 | no |
| gj-3306 | tx | linear_astrometry_v1 | 0.0820 | 0.0632 | 0.1035 | no |
| gj-338-a | rx | linear_astrometry_v1 | 0.0004 | 0.0004 | 0.0004 | no |
| gj-338-a | tx | linear_astrometry_v1 | 0.0051 | 0.0222 | 0.0227 | no |
| gj-338-b | rx | linear_astrometry_v1 | 0.0004 | 0.0004 | 0.0004 | no |
| gj-338-b | tx | linear_astrometry_v1 | 0.0046 | 0.0246 | 0.0250 | no |
| gj-3512 | rx | linear_astrometry_v1 | 0.0005 | 0.0005 | 0.0006 | no |
| gj-3512 | tx | linear_astrometry_v1 | 0.0335 | 0.0540 | 0.0632 | no |
| gj-367 | rx | linear_astrometry_v1 | 0.0003 | 0.0003 | 0.0004 | no |
| gj-367 | tx | linear_astrometry_v1 | 0.0152 | 0.0037 | 0.0157 | no |
| gj-514 | rx | linear_astrometry_v1 | 0.0004 | 0.0006 | 0.0006 | no |
| gj-514 | tx | linear_astrometry_v1 | 0.0132 | 0.0386 | 0.0407 | no |
| gj-518 | rx | linear_astrometry_v1 | 0.0026 | 0.0031 | 0.0041 | no |
| gj-518 | tx | linear_astrometry_v1 | 0.1712 | 0.1963 | 0.2605 | no |
| gj-526 | rx | linear_astrometry_v1 | 0.0004 | 0.0007 | 0.0008 | no |
| gj-526 | tx | linear_astrometry_v1 | 0.0065 | 0.0281 | 0.0288 | no |
| gj-54 | rx | linear_astrometry_v1 | 0.0011 | 0.0008 | 0.0012 | no |
| gj-54 | tx | linear_astrometry_v1 | 0.0045 | 0.0115 | 0.0116 | no |
| gj-581 | rx | linear_astrometry_v1 | 0.0006 | 0.0008 | 0.0009 | no |
| gj-581 | tx | linear_astrometry_v1 | 0.0116 | 0.0218 | 0.0242 | no |
| gj-588 | rx | linear_astrometry_v1 | 0.0004 | 0.0006 | 0.0007 | no |
| gj-588 | tx | linear_astrometry_v1 | 0.0134 | 0.0208 | 0.0249 | no |
| gj-625 | rx | linear_astrometry_v1 | 0.0004 | 0.0005 | 0.0005 | no |
| gj-625 | tx | linear_astrometry_v1 | 0.0042 | 0.0052 | 0.0062 | no |
| gj-66-a | rx | linear_astrometry_v1 | 0.0007 | 0.0007 | 0.0008 | no |
| gj-66-a | tx | linear_astrometry_v1 | 0.0066 | 0.0120 | 0.0126 | no |
| gj-66-b | rx | linear_astrometry_v1 | 0.0007 | 0.0006 | 0.0008 | no |
| gj-66-b | tx | linear_astrometry_v1 | 0.0063 | 0.0114 | 0.0123 | no |
| gj-667-c | rx | linear_astrometry_v1 | 0.0004 | 0.0005 | 0.0006 | no |
| gj-667-c | tx | linear_astrometry_v1 | 0.0188 | 0.0241 | 0.0302 | no |
| gj-674 | rx | linear_astrometry_v1 | 0.0008 | 0.0006 | 0.0009 | no |
| gj-674 | tx | linear_astrometry_v1 | 0.0096 | 0.0043 | 0.0099 | no |
| gj-682 | rx | linear_astrometry_v1 | 0.0006 | 0.0006 | 0.0008 | no |
| gj-682 | tx | linear_astrometry_v1 | 0.0033 | 0.0154 | 0.0154 | no |
| gj-687 | rx | linear_astrometry_v1 | 0.0005 | 0.0006 | 0.0006 | no |
| gj-687 | tx | linear_astrometry_v1 | 0.0051 | 0.0090 | 0.0105 | no |
| gj-783 | rx | linear_astrometry_v1 | 0.0018 | 0.0019 | 0.0021 | no |
| gj-783 | tx | linear_astrometry_v1 | 0.0776 | 0.0799 | 0.1128 | no |
| gj-784 | rx | linear_astrometry_v1 | 0.0004 | 0.0004 | 0.0005 | no |
| gj-784 | tx | linear_astrometry_v1 | 0.0067 | 0.0105 | 0.0123 | no |
| gj-832 | rx | linear_astrometry_v1 | 0.0004 | 0.0004 | 0.0005 | no |
| gj-832 | tx | linear_astrometry_v1 | 0.0080 | 0.0036 | 0.0086 | no |
| gj-876 | rx | linear_astrometry_v1 | 0.0007 | 0.0009 | 0.0009 | no |
| gj-876 | tx | linear_astrometry_v1 | 0.0126 | 0.0100 | 0.0160 | no |
| gj-908 | rx | linear_astrometry_v1 | 0.0005 | 0.0007 | 0.0008 | no |
| gj-908 | tx | linear_astrometry_v1 | 0.0215 | 0.0094 | 0.0232 | no |
| gj-915 | rx | linear_astrometry_v1 | 0.0009 | 0.0006 | 0.0011 | no |
| gj-915 | tx | linear_astrometry_v1 | 0.0476 | 0.0284 | 0.0550 | no |
| gj-9193 | rx | linear_astrometry_v1 | 0.0030 | 0.0041 | 0.0051 | no |
| gj-9193 | tx | linear_astrometry_v1 | 0.1031 | 0.1405 | 0.1743 | no |
| gj65-a | rx | two_body_orbit_v1 | 0.4030 | 0.4124 | 0.4739 | no |
| gj65-a | tx | two_body_orbit_v1 | 1.2327 | 1.3187 | 1.5341 | no |
| gj65-b | rx | two_body_orbit_v1 | 0.3887 | 0.4145 | 0.4734 | no |
| gj65-b | tx | two_body_orbit_v1 | 1.2462 | 1.3177 | 1.5372 | no |
| groombridge-34-a | rx | linear_astrometry_v1 | 0.0003 | 0.0003 | 0.0003 | no |
| groombridge-34-a | tx | linear_astrometry_v1 | 0.0031 | 0.0131 | 0.0134 | no |
| groombridge-34-b | rx | linear_astrometry_v1 | 0.0004 | 0.0005 | 0.0005 | no |
| groombridge-34-b | tx | linear_astrometry_v1 | 0.0047 | 0.0176 | 0.0181 | no |
| hd-219134 | rx | linear_astrometry_v1 | 0.0009 | 0.0010 | 0.0012 | no |
| hd-219134 | tx | linear_astrometry_v1 | 0.0297 | 0.0681 | 0.0735 | no |
| kapteyn-star | rx | linear_astrometry_v1 | 0.0004 | 0.0004 | 0.0005 | no |
| kapteyn-star | tx | linear_astrometry_v1 | 0.0361 | 0.0127 | 0.0382 | no |
| lacaille-8760 | rx | linear_astrometry_v1 | 0.0005 | 0.0006 | 0.0007 | no |
| lacaille-8760 | tx | linear_astrometry_v1 | 0.0190 | 0.0275 | 0.0332 | no |
| lacaille-9352 | rx | linear_astrometry_v1 | 0.0003 | 0.0004 | 0.0005 | no |
| lacaille-9352 | tx | linear_astrometry_v1 | 0.0024 | 0.0246 | 0.0247 | no |
| lalande-21185 | rx | linear_astrometry_v1 | 0.0006 | 0.0005 | 0.0007 | no |
| lalande-21185 | tx | linear_astrometry_v1 | 0.0174 | 0.0017 | 0.0174 | no |
| lhs-1723 | rx | linear_astrometry_v1 | 0.0005 | 0.0005 | 0.0006 | no |
| lhs-1723 | tx | linear_astrometry_v1 | 0.0029 | 0.0110 | 0.0111 | no |
| lp-145-141 | rx | linear_astrometry_v1 | 0.0021 | 0.0077 | 0.0080 | no |
| lp-145-141 | tx | linear_astrometry_v1 | 0.0413 | 0.1482 | 0.1541 | no |
| ltt-1445-a | rx | linear_astrometry_v1 | 0.0006 | 0.0006 | 0.0007 | no |
| ltt-1445-a | tx | linear_astrometry_v1 | 0.0039 | 0.0098 | 0.0099 | no |
| luhman16-a | rx | two_body_orbit_v1 | 0.3516 | 0.6439 | 0.6565 | no |
| luhman16-a | tx | two_body_orbit_v1 | 0.3951 | 0.7570 | 0.7701 | no |
| luhman16-b | rx | two_body_orbit_v1 | 0.4494 | 0.7924 | 0.8024 | no |
| luhman16-b | tx | two_body_orbit_v1 | 0.4996 | 0.8683 | 0.8797 | no |
| luyten-star | rx | linear_astrometry_v1 | 0.0009 | 0.0008 | 0.0011 | no |
| luyten-star | tx | linear_astrometry_v1 | 0.0129 | 0.0350 | 0.0371 | no |
| procyon-a | rx | two_body_orbit_v1 | 1.6299 | 1.6227 | 1.8761 | no |
| procyon-a | tx | two_body_orbit_v1 | 2.7759 | 2.7860 | 3.2089 | no |
| procyon-b | rx | two_body_orbit_v1 | 1.6255 | 1.6196 | 1.8734 | no |
| procyon-b | tx | two_body_orbit_v1 | 2.8652 | 2.7657 | 3.2634 | no |
| proxima-cen | rx | linear_astrometry_v1 | 0.0011 | 0.0007 | 0.0012 | no |
| proxima-cen | tx | linear_astrometry_v1 | 0.0031 | 0.0050 | 0.0057 | no |
| ross-128 | rx | linear_astrometry_v1 | 0.0005 | 0.0007 | 0.0008 | no |
| ross-128 | tx | linear_astrometry_v1 | 0.0053 | 0.0068 | 0.0083 | no |
| ross-154 | rx | linear_astrometry_v1 | 0.0007 | 0.0008 | 0.0009 | no |
| ross-154 | tx | linear_astrometry_v1 | 0.0024 | 0.0040 | 0.0041 | no |
| ross-248 | rx | linear_astrometry_v1 | 0.0006 | 0.0007 | 0.0008 | no |
| ross-248 | tx | linear_astrometry_v1 | 0.0119 | 0.0053 | 0.0128 | no |
| sigma-dra | rx | linear_astrometry_v1 | 0.0074 | 0.0061 | 0.0083 | no |
| sigma-dra | tx | linear_astrometry_v1 | 0.1252 | 0.2118 | 0.2434 | no |
| sirius-a | rx | two_body_orbit_v1 | 1.7275 | 1.7013 | 1.9936 | no |
| sirius-a | tx | two_body_orbit_v1 | 2.6366 | 2.6236 | 3.0965 | no |
| sirius-b | rx | two_body_orbit_v1 | 1.7377 | 1.6920 | 1.9990 | no |
| sirius-b | tx | two_body_orbit_v1 | 2.6930 | 2.6663 | 3.1473 | no |
| struve-2398-a | rx | linear_astrometry_v1 | 0.0006 | 0.0005 | 0.0007 | no |
| struve-2398-a | tx | linear_astrometry_v1 | 0.0065 | 0.0087 | 0.0106 | no |
| struve-2398-b | rx | linear_astrometry_v1 | 0.0008 | 0.0007 | 0.0009 | no |
| struve-2398-b | tx | linear_astrometry_v1 | 0.0092 | 0.0116 | 0.0146 | no |
| tau-cet | rx | linear_astrometry_v1 | 0.0069 | 0.0106 | 0.0118 | no |
| tau-cet | tx | linear_astrometry_v1 | 0.1874 | 0.1429 | 0.2308 | no |
| teegarden | rx | linear_astrometry_v1 | 0.0020 | 0.0018 | 0.0023 | no |
| teegarden | tx | linear_astrometry_v1 | 0.1107 | 0.0548 | 0.1235 | no |
| van-maanen | rx | linear_astrometry_v1 | 0.0060 | 0.0006 | 0.0060 | no |
| van-maanen | tx | linear_astrometry_v1 | 0.1032 | 0.0034 | 0.1032 | no |
| wise-0855 | rx | linear_astrometry_v1 | 0.7097 | 0.7260 | 0.8453 | no |
| wise-0855 | tx | linear_astrometry_v1 | 2.0804 | 2.2788 | 2.6567 | no |
| wolf-1061 | rx | linear_astrometry_v1 | 0.0005 | 0.0007 | 0.0007 | no |
| wolf-1061 | tx | linear_astrometry_v1 | 0.0094 | 0.0043 | 0.0100 | no |
| wolf-1069 | rx | linear_astrometry_v1 | 0.0021 | 0.0018 | 0.0023 | no |
| wolf-1069 | tx | linear_astrometry_v1 | 0.0186 | 0.0762 | 0.0763 | no |
| wolf-359 | rx | linear_astrometry_v1 | 0.0014 | 0.0017 | 0.0019 | no |
| wolf-359 | tx | linear_astrometry_v1 | 0.0272 | 0.0180 | 0.0329 | no |
| wolf-437 | rx | linear_astrometry_v1 | 0.0007 | 0.0009 | 0.0009 | no |
| wolf-437 | tx | linear_astrometry_v1 | 0.0302 | 0.0284 | 0.0415 | no |

# Independent ephemeris check (plan §2.4)

astropy-only anti-star locus (SkyCoord.apply_space_motion for the star direction — at t for Rx, at t + 2d/c for Tx — plus Earth-centre parallax) for the 156 linear-astrometry endpoint-roles, vs the sglseti `tusay2022_eq5_7_v1` locus with the same observer: median residual 0.0127", max 3.2904". The Tx loci agree to ≲ 0.05"; the Rx residual grows linearly with z and with the star's proper motion (≈ 2 µ z / c, the relay light-time term of the model, which the naive construction omits) — 0.18" at 550 AU and 3.3" at 10,000 AU for Barnard's Star (µ = 10.4"/yr), ≲ 0.4" for µ ≲ 1.3"/yr. Both are below the W1 PSF and are properties of the model, not of its implementation.

| endpoint | role | max residual ["] |
| --- | --- | --- |
| 61-cyg-a | rx | 1.6710 |
| 61-cyg-a | tx | 0.0185 |
| 61-cyg-b | rx | 1.6383 |
| 61-cyg-b | tx | 0.0177 |
| 61-vir | rx | 0.4772 |
| 61-vir | tx | 0.0006 |
| 82-eri | rx | 0.9871 |
| 82-eri | tx | 0.0147 |
| barnard-star | rx | 3.2904 |
| barnard-star | tx | 0.0614 |
| eps-eri | rx | 0.3083 |
| eps-eri | tx | 0.0009 |
| eps-ind-a | rx | 1.4893 |
| eps-ind-a | tx | 0.0100 |
| eps-ind-b | rx | 1.4817 |
| eps-ind-b | tx | 0.0100 |
| ez-aqr | rx | 1.0310 |
| ez-aqr | tx | 0.0001 |
| fomalhaut | rx | 0.1163 |
| fomalhaut | tx | 0.0005 |
| gj-1002 | rx | 0.6515 |
| gj-1002 | tx | 0.0044 |
| gj-1061 | rx | 0.2637 |
| gj-1061 | tx | 0.0001 |
| gj-1087 | rx | 0.3249 |
| gj-1087 | tx | 0.0228 |
| gj-11068 | rx | 0.0396 |
| gj-11068 | tx | 0.0006 |
| gj-1111 | rx | 0.4019 |
| gj-1111 | tx | 0.0008 |
| gj-11547 | rx | 0.4866 |
| gj-11547 | tx | 0.0009 |
| gj-1221 | rx | 0.5323 |
| gj-1221 | tx | 0.0001 |
| gj-12724 | rx | 0.8380 |
| gj-12724 | tx | 0.0025 |
| gj-1276 | rx | 0.8152 |
| gj-1276 | tx | 0.0002 |
| gj-13157 | rx | 0.0281 |
| gj-13157 | tx | 0.0000 |
| gj-2012 | rx | 0.1898 |
| gj-2012 | tx | 0.0000 |
| gj-2066 | rx | 0.1205 |
| gj-2066 | tx | 0.0013 |
| gj-229-a | rx | 0.2314 |
| gj-229-a | tx | 0.0002 |
| gj-251 | rx | 0.2620 |
| gj-251 | tx | 0.0010 |
| gj-293 | rx | 0.6612 |
| gj-293 | tx | 0.0001 |
| gj-3112 | rx | 0.3385 |
| gj-3112 | tx | 0.0001 |
| gj-318 | rx | 0.5420 |
| gj-318 | tx | 0.0001 |
| gj-3306 | rx | 0.4967 |
| gj-3306 | tx | 0.0001 |
| gj-338-a | rx | 0.5209 |
| gj-338-a | tx | 0.0010 |
| gj-338-b | rx | 0.5394 |
| gj-338-b | tx | 0.0012 |
| gj-3512 | rx | 0.4129 |
| gj-3512 | tx | 0.0006 |
| gj-367 | rx | 0.2353 |
| gj-367 | tx | 0.0019 |
| gj-514 | rx | 0.4924 |
| gj-514 | tx | 0.0013 |
| gj-518 | rx | 1.2328 |
| gj-518 | tx | 0.0003 |
| gj-526 | rx | 0.7261 |
| gj-526 | tx | 0.0020 |
| gj-54 | rx | 0.4312 |
| gj-54 | tx | 0.0009 |
| gj-581 | rx | 0.3875 |
| gj-581 | tx | 0.0006 |
| gj-588 | rx | 0.4948 |
| gj-588 | tx | 0.0018 |
| gj-625 | rx | 0.1471 |
| gj-625 | tx | 0.0003 |
| gj-66-a | rx | 0.0832 |
| gj-66-a | tx | 0.0003 |
| gj-66-b | rx | 0.0979 |
| gj-66-b | tx | 0.0004 |
| gj-667-c | rx | 0.3643 |
| gj-667-c | tx | 0.0004 |
| gj-674 | rx | 0.3322 |
| gj-674 | tx | 0.0002 |
| gj-682 | rx | 0.3714 |
| gj-682 | tx | 0.0021 |
| gj-687 | rx | 0.4143 |
| gj-687 | tx | 0.0020 |
| gj-783 | rx | 0.5196 |
| gj-783 | tx | 0.0113 |
| gj-784 | rx | 0.2513 |
| gj-784 | tx | 0.0014 |
| gj-832 | rx | 0.2588 |
| gj-832 | tx | 0.0006 |
| gj-876 | rx | 0.3703 |
| gj-876 | tx | 0.0001 |
| gj-908 | rx | 0.4388 |
| gj-908 | tx | 0.0052 |
| gj-915 | rx | 0.2914 |
| gj-915 | tx | 0.0001 |
| gj-9193 | rx | 0.7520 |
| gj-9193 | tx | 0.0001 |
| groombridge-34-a | rx | 0.9237 |
| groombridge-34-a | tx | 0.0019 |
| groombridge-34-b | rx | 0.9116 |
| groombridge-34-b | tx | 0.0017 |
| hd-219134 | rx | 0.6626 |
| hd-219134 | tx | 0.0020 |
| kapteyn-star | rx | 2.7358 |
| kapteyn-star | tx | 0.1128 |
| lacaille-8760 | rx | 1.0927 |
| lacaille-8760 | tx | 0.0039 |
| lacaille-9352 | rx | 2.1809 |
| lacaille-9352 | tx | 0.0032 |
| lalande-21185 | rx | 1.5224 |
| lalande-21185 | tx | 0.0218 |
| lhs-1723 | rx | 0.2428 |
| lhs-1723 | tx | 0.0018 |
| lp-145-141 | rx | 0.8488 |
| lp-145-141 | tx | 0.0001 |
| ltt-1445-a | rx | 0.1445 |
| ltt-1445-a | tx | 0.0001 |
| luyten-star | rx | 1.1813 |
| luyten-star | tx | 0.0036 |
| proxima-cen | rx | 1.2209 |
| proxima-cen | tx | 0.0045 |
| ross-128 | rx | 0.4319 |
| ross-128 | tx | 0.0022 |
| ross-154 | rx | 0.2113 |
| ross-154 | tx | 0.0004 |
| ross-248 | rx | 0.5048 |
| ross-248 | tx | 0.0066 |
| sigma-dra | rx | 0.5813 |
| sigma-dra | tx | 0.0027 |
| struve-2398-a | rx | 0.7024 |
| struve-2398-a | tx | 0.0002 |
| struve-2398-b | rx | 0.7369 |
| struve-2398-b | tx | 0.0002 |
| tau-cet | rx | 0.6080 |
| tau-cet | tx | 0.0016 |
| teegarden | rx | 1.6203 |
| teegarden | tx | 0.0174 |
| van-maanen | rx | 0.9426 |
| van-maanen | tx | 0.0417 |
| wise-0855 | rx | 2.5777 |
| wise-0855 | tx | 0.0002 |
| wolf-1061 | rx | 0.3757 |
| wolf-1061 | tx | 0.0013 |
| wolf-1069 | rx | 0.1905 |
| wolf-1069 | tx | 0.0019 |
| wolf-359 | rx | 1.4913 |
| wolf-359 | tx | 0.0050 |
| wolf-437 | rx | 0.3505 |
| wolf-437 | tx | 0.0012 |
