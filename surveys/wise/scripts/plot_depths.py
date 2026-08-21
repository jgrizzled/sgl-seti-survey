"""Regenerate report/figures/depth_by_endpoint.svg from the calibration
m90 curves (duty 0.5, 90% recovery). Four panels (W1..W4); one row per
endpoint sorted by W1 Rx median depth; whisker = range over 550–10,000
AU; filled marker = Rx median, open marker = Tx median. Single hue;
role is carried by marker fill, not colour.

Usage: uv run python surveys/wise/scripts/plot_depths.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
M90 = REPO / "runs" / "wise" / "calib_v1" / "m90_curves.npz"
OUT = REPO / "report" / "figures" / "depth_by_endpoint.svg"

ROW = 16
LEFT, TOP = 158, 46
PANEL_W, GAP_X, GAP_Y = 410, 54, 70
BANDS = [("W1", "3.4 μm", 11, 18), ("W2", "4.6 μm", 10, 18),
         ("W3", "12 μm", 7, 15), ("W4", "22 μm", 4, 13)]


def main():
    m = np.load(M90, allow_pickle=True)
    eps = sorted({k.split("__")[0] for k in m.files})

    def stats(e, role, band):
        k = f"{e}__{role}__{band}__0.5"
        if k not in m.files:
            return None
        c = m[k][np.isfinite(m[k])]
        if len(c) == 0:
            return None
        return float(np.median(c)), float(c.min()), float(c.max())

    order = sorted(eps, key=lambda e: -(stats(e, "rx", "W1") or (0,))[0])
    n = len(order)
    panel_h = n * ROW + 12
    width = LEFT + 2 * PANEL_W + GAP_X + 40
    height = TOP + 2 * panel_h + GAP_Y + 70
    L = [f'<svg viewBox="0 0 {width} {height}" width="{width}" '
         'xmlns="http://www.w3.org/2000/svg" role="img" '
         'aria-label="90 percent recovery depth by endpoint and band">',
         '<style>text{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;fill:#333}'
         '.pt{font-size:15px;font-weight:600}.ps{font-size:13px;fill:#777;font-weight:400}'
         '.tk{font-size:10.5px;fill:#999}.ax{font-size:11px;fill:#666}.rl{font-size:10px}'
         '.gr{stroke:#e8e8e4;stroke-width:1}.wh{stroke:#c4c4be;stroke-width:1.8;stroke-linecap:round}'
         '.rx{fill:#b4551f}.tx{fill:#fdfdfc;stroke:#b4551f;stroke-width:1.6}</style>',
         f'<rect x="0" y="0" width="{width}" height="{height}" fill="#fdfdfc" rx="8"/>',
         f'<text x="{LEFT}" y="22" class="pt">90% recovery depth by endpoint and band'
         '<tspan class="ps"> · duty ≥ 0.5 · filled = Rx, open = Tx, bar = range over '
         f'550–10,000 AU · calibration v0.2.0 · {n} endpoints</tspan></text>']
    for bi, (band, lam, lo, hi) in enumerate(BANDS):
        px = LEFT + (bi % 2) * (PANEL_W + GAP_X)
        py = TOP + (bi // 2) * (panel_h + GAP_Y)
        sx = lambda v, lo=lo, hi=hi, px=px: px + (v - lo) / (hi - lo) * PANEL_W
        L.append(f'<text x="{px}" y="{py - 8}" class="pt">{band}'
                 f'<tspan class="ps"> · {lam}</tspan></text>')
        for t in range(lo, hi + 1):
            L.append(f'<line x1="{sx(t):.1f}" y1="{py}" x2="{sx(t):.1f}" '
                     f'y2="{py + panel_h - 12}" class="gr"/>')
            L.append(f'<text x="{sx(t):.1f}" y="{py + panel_h + 2}" class="tk" '
                     f'text-anchor="middle">{t}</text>')
        L.append(f'<text x="{px + PANEL_W / 2:.0f}" y="{py + panel_h + 20}" '
                 'class="ax" text-anchor="middle">Vega mag → deeper</text>')
        for i, e in enumerate(order):
            y = py + 8 + i * ROW
            if bi % 2 == 0:
                L.append(f'<text x="{px - 8}" y="{y + 3.5}" class="rl" '
                         f'text-anchor="end">{e}</text>')
            rx, tx = stats(e, "rx", band), stats(e, "tx", band)
            if rx:
                a, b = max(rx[1], lo), min(rx[2], hi)
                L.append(f'<line x1="{sx(a):.1f}" y1="{y}" x2="{sx(b):.1f}" '
                         f'y2="{y}" class="wh"/>')
                L.append(f'<circle cx="{sx(min(max(rx[0], lo), hi)):.1f}" cy="{y}" r="3.1" class="rx"/>')
            if tx:
                L.append(f'<circle cx="{sx(min(max(tx[0], lo), hi)):.1f}" cy="{y}" r="3.1" class="tx"/>')
    L.append("</svg>")
    OUT.write_text("\n".join(L))
    print(f"wrote {OUT.relative_to(REPO)}: {n} endpoints")


if __name__ == "__main__":
    main()
