"""Depth-vs-relay-distance figure for the Pan-STARRS1 pilot report (pure SVG,
no plotting dependency — same approach as the WISE report figure).

Usage: uv run python surveys/panstarrs/scripts/make_figures.py
"""

from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
CAL = REPO / "runs" / "panstarrs" / "calib_v1"
OUT = REPO / "report" / "figures" / "ps1_survey_v1_m90.svg"
PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b",
           "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


class _Colors(dict):
    def __missing__(self, k):
        self[k] = PALETTE[len(self) % len(PALETTE)]
        return self[k]


COLORS = _Colors()
SHOW_ENDPOINTS = {"ross-128", "proxima-cen", "eps-ind-a", "eps-ind-b",
                  "barnard-star", "sirius-a", "alpha-cen-a", "82-eri"}


def main() -> None:
    d = np.load(CAL / "m90_curves.npz")
    z = np.load(CAL / "tensors" / "ross-128__rx.npz")["z_grid"]
    W, H, L, R, T, B = 760, 460, 70, 20, 40, 60
    zmin, zmax, mmin, mmax = 550.0, 10000.0, 18.0, 24.0

    def X(zz):
        return L + (np.log10(zz) - np.log10(zmin)) / (np.log10(zmax) - np.log10(zmin)) * (W - L - R)

    def Y(m):
        return T + (mmax - m) / (mmax - mmin) * (H - T - B)

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{W}" height="{H}" fill="white"/>']
    for zz in (550, 1000, 2000, 5000, 10000):
        out.append(f'<line x1="{X(zz):.1f}" y1="{T}" x2="{X(zz):.1f}" y2="{H-B}" stroke="#ddd"/>')
        out.append(f'<text x="{X(zz):.1f}" y="{H-B+16}" text-anchor="middle">{zz}</text>')
    for m in range(int(mmin), int(mmax) + 1):
        out.append(f'<line x1="{L}" y1="{Y(m):.1f}" x2="{W-R}" y2="{Y(m):.1f}" stroke="#ddd"/>')
        out.append(f'<text x="{L-6}" y="{Y(m)+4:.1f}" text-anchor="end">{m}</text>')
    out.append(f'<text x="{(L+W-R)/2:.0f}" y="{H-B+36}" text-anchor="middle">relay distance z [AU]</text>')
    out.append(f'<text transform="translate(16,{(T+H-B)/2:.0f}) rotate(-90)" text-anchor="middle">m90 [AB mag], 90% recovery, duty ≥ 0.5 (fainter ↑)</text>')
    out.append(f'<text x="{L}" y="{T-18}" font-size="14" font-weight="bold">PS1 survey v1 — injection-recovery depths, selected endpoints (star-calibrated ZP)</text>')
    legend_y = T + 10
    legend_x = L + 10
    for k in sorted(d.files):
        if not k.endswith("__0.5"):
            continue
        e, r, b, _ = k.split("__")
        m = d[k]
        if not np.isfinite(m).any() or e not in SHOW_ENDPOINTS:
            continue
        pts = [(X(zz), Y(np.clip(mm, mmin, mmax))) for zz, mm in zip(z, m)
               if np.isfinite(mm)]
        # break the polyline at gaps
        segs, cur, last = [], [], None
        for zz, mm in zip(z, m):
            if np.isfinite(mm):
                cur.append(f"{X(zz):.1f},{Y(np.clip(mm, mmin, mmax)):.1f}")
            elif cur:
                segs.append(cur); cur = []
        if cur:
            segs.append(cur)
        dash = {"g": "", "r": " stroke-dasharray=\"6,4\"", "i": " stroke-dasharray=\"2,3\"", "z": " stroke-dasharray=\"8,3,2,3\"", "y": " stroke-dasharray=\"1,2\""}[b]
        op = "1.0" if r == "rx" else "0.45"
        for sg in segs:
            out.append(f'<polyline points="{" ".join(sg)}" fill="none" stroke="{COLORS[e]}" stroke-width="2" opacity="{op}"{dash}/>')
        out.append(f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x+30}" y2="{legend_y}" stroke="{COLORS[e]}" stroke-width="2" opacity="{op}"{dash}/>')
        out.append(f'<text x="{legend_x+36}" y="{legend_y+4}" font-size="10">{e} {r} {b}</text>')
        legend_y += 13
        if legend_y > H - B - 10:
            legend_y, legend_x = T + 10, legend_x + 150
    out.append("</svg>")
    OUT.write_text("\n".join(out))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
