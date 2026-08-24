"""v2 engine: re-score a v1 asteroid positive control with the v2 rule.

The v1 control scripts fetched cutouts centred on the Horizons position
of a numbered asteroid at each frame and stacked a 5 x 5 residual grid
against 8 offset controls. v2 re-samples the SAME retained cutouts on
the 48-offset ring, applies the per-frame cap and the single-epoch
clip of the profile, and reports R, R~ = R/q95(ring), the rank
statement, the peak offset and the recovered stack magnitude against
the prediction — the blind recovery test of hypotheses §1.5."""

from __future__ import annotations

import json

import numpy as np

from sglsurvey import nulls
from sglsurvey.v2.build import sample_pix

RESID_GRID = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])


def rescore(P, frames: list[dict], maps: list, positions: list[tuple], out_dir, label: str,
            predicted_mag: dict | None = None, extra: dict | None = None, clip_sigma="profile") -> dict:
    """``frames``: ok-frame dicts (band, mjd, v_pred); ``maps``: FluxMaps
    aligned with them; ``positions``: (ra, dec) of the asteroid per
    frame; ``predicted_mag``: band -> predicted magnitude at ZP_REF
    scale (from Horizons V plus colour terms) for the throughput line."""
    ring = np.array(P.ring)
    offsets = np.vstack([[0.0, 0.0], ring])
    nr, nt = len(RESID_GRID), len(offsets)
    F, V, G, band_of, wf = [], [], [], [], []
    for fr, fm, (ra, dec) in zip(frames, maps, positions):
        if fm is None or fm.magzp is None:
            continue
        scale = 10.0 ** ((P.zp_ref - fm.magzp) / 2.5)
        cosd = np.cos(np.deg2rad(dec))
        pp = fm.world2pix(np.array([ra, ra + 10 / 3600 / cosd, ra]), np.array([dec, dec, dec + 10 / 3600]))
        J = np.stack([(pp[1] - pp[0]) / 10.0, (pp[2] - pp[0]) / 10.0], axis=1); p0 = pp[0]
        off = np.empty((nt, nr, nr, 2))
        off[..., 0] = RESID_GRID[None, :, None] + offsets[:, 0][:, None, None]
        off[..., 1] = RESID_GRID[None, None, :] + offsets[:, 1][:, None, None]
        pix = p0 + np.einsum("ij,...j->...i", J, off)
        f, v, g = sample_pix(fm, pix[..., 0].ravel(), pix[..., 1].ravel())
        F.append(f.reshape(nt, nr, nr) * scale); V.append(np.minimum(v.reshape(nt, nr, nr) * scale ** 2, 1e30))
        G.append(g.reshape(nt, nr, nr)); band_of.append(fr["band"])
        sel = np.isfinite(fm.var) & (fm.good_frac >= P.min_good_frac) & (fm.var > 0)
        wf.append(1.0 / (np.median(fm.var[sel]) * scale * scale) if sel.sum() > 100 else np.nan)
    F, V, G = np.array(F), np.array(V), np.array(G)
    band_of = np.array(band_of); wf = np.array(wf)
    designated = np.array([1 + i for i in P.designated])
    clip = P.clip_sigma if clip_sigma == "profile" else clip_sigma
    summary = {"asteroid": label, "n_frames_ok": int(len(F)), "rule": "v2 ring-48 / per-frame cap", "clip_sigma": clip}
    if extra:
        summary.update(extra)
    for band in P.bands:
        sel = band_of == band
        if sel.sum() < P.min_epochs:
            continue
        cap = P.weight_cap * float(np.nanmedian(wf[sel]))
        fc = np.full(int(sel.sum()), cap, np.float32)
        f, v, g = F[sel].transpose(1, 0, 2, 3), V[sel].transpose(1, 0, 2, 3), G[sel].transpose(1, 0, 2, 3)
        mx = np.array([nulls.grid_max(nulls.stack_S(f[t], v[t], g[t], frame_cap=fc, clip_sigma=clip))[0]
                       for t in range(nt)])
        R, T = nulls.exceedance_ratios(mx, designated)
        S0, A0, B0, n0, _ = nulls.stack_S(f[0], v[0], g[0], return_parts=True, frame_cap=fc, clip_sigma=clip)
        s_max, node = nulls.grid_max(S0)
        if node is None:
            continue
        flux = A0[node] / B0[node]
        mag = P.zp_ref - 2.5 * np.log10(flux) if flux > 0 else None
        ring_R = R[1:]
        q95 = float(np.nanquantile(ring_R, P.norm_quantile))
        rs = nulls.rank_statement(R[0], ring_R)
        pm = (predicted_mag or {}).get(band)
        summary[band] = {
            "n_epochs": int(sel.sum()), "S_max_real": float(s_max), "S_centre": float(S0[2, 2]),
            "T_designated": float(T), "R_real": float(R[0]), "ring_R_median": float(np.nanmedian(ring_R)),
            "ring_R_q95": q95, "R_norm_real": float(R[0] / q95) if q95 > 0 else None, "rank_statement": rs,
            "peak_offset_arcsec": [float(RESID_GRID[node[0]]), float(RESID_GRID[node[1]])],
            "recovered_stack_mag": mag, "predicted_mag": pm,
            "throughput": {"median_dmag": (mag - pm) if (mag is not None and pm is not None) else None},
            "throughput_note": "(recovered stack mag minus Horizons-predicted mag with solar colours)",
        }
        print(f"[{label} {band}] N={sel.sum()} S_max={s_max:.1f} T={T:.2f} R={R[0]:.2f} q95={q95:.2f} "
              f"R~={summary[band]['R_norm_real']:.2f} rank p={rs['p_cell']:.3f} mag {mag if mag is None else round(mag, 2)} "
              f"vs pred {pm}", flush=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    return summary
