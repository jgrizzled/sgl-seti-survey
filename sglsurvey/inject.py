"""Image-level source injection for WISE L1b cutouts (v2 plan §4).

Pieces:

* :class:`PRFGrid` — the IRSA second-pass single-exposure PSF templates
  (``wise-w?-psf-wpro-09x09-XXxYY.fits``: a 9 x 9 focal-plane grid of
  641 x 641 unit-sum templates sampled at 1/8 detector pixel, centre at
  1-based CRPIX 321). ``render`` places a template at a sub-pixel
  position of a cutout, returning a stamp in detector pixels whose sum
  is the fraction of the unit source flux falling inside the stamp.
* spectrum → DN: WISE Vega magnitude of a source with a declared
  spectrum and flux density at the band's isophotal wavelength, via the
  Explanatory Supplement §IV.4.h colour corrections (Table 2, Wright et
  al. 2010) and zero-magnitude flux densities (Table 1, F*_nu0), then
  DN through the frame MAGZP.
* :func:`stamp_response` — the exact change of the matched-filter flux
  map caused by adding a stamp to the image, computed locally from the
  map's stored kernel / usable-pixel mask / denominator (the estimator
  is linear in the image; only the background median is not, and a
  point source moves a sigma-clipped median by a negligible amount —
  measured by the invariant tests).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PRF_OVERSAMPLE = 8
PRF_GRID = 9
#: Native frame size per band (W4 is 508 x 508 at 5.52"/pixel).
FRAME_PIX = {"W1": 1016, "W2": 1016, "W3": 1016, "W4": 508}
#: Template size and 0-based centre (W1-W3: 641 px, CRPIX 321; W4: 127 px,
#: CRPIX 64) — both 8x oversampled in their band's native pixel.
PRF_SHAPE = {"W1": 641, "W2": 641, "W3": 641, "W4": 127}
PRF_CENTRE = {"W1": 320, "W2": 320, "W3": 320, "W4": 63}

#: Explanatory Supplement IV.4.h Table 1, F*_nu0 [Jy] (for use with f_c).
F_NU0_STAR_JY = {"W1": 306.682, "W2": 170.663, "W3": 29.045, "W4": 8.284}
#: Isophotal wavelengths [um] (Wright et al. 2010 / Jarrett et al. 2011).
LAMBDA_ISO_UM = {"W1": 3.3526, "W2": 4.6028, "W3": 11.5608, "W4": 22.0883}
#: Vega -> AB offsets (Table 3).
AB_OFFSET = {"W1": 2.699, "W2": 3.339, "W3": 5.174, "W4": 6.620}
#: Colour corrections f_c (Table 2): spectrum -> per-band factor.
COLOR_CORRECTIONS = {
    "powerlaw_nu-2":   {"W1": 1.0000, "W2": 1.0000, "W3": 1.0000, "W4": 1.0000},
    "flat_fnu":        {"W1": 0.9907, "W2": 0.9935, "W3": 0.9169, "W4": 0.9905},
    "powerlaw_nu-1":   {"W1": 0.9921, "W2": 0.9943, "W3": 0.9373, "W4": 0.9926},
    "g2v":             {"W1": 1.0049, "W2": 1.0193, "W3": 1.0024, "W4": 1.0012},
    "k2v":             {"W1": 1.0038, "W2": 1.0512, "W3": 1.0030, "W4": 1.0013},
    "blackbody_283K":  {"W1": 1.3917, "W2": 1.1124, "W3": 0.8791, "W4": 0.9865},
    "blackbody_400K":  {"W1": 1.1316, "W2": 1.0229, "W3": 0.8622, "W4": 0.9903},
}


def _bb_300k():
    """300 K is not tabulated: interpolate 283 K and 400 K in log T."""
    t = math.log(300 / 283) / math.log(400 / 283)
    a, b = COLOR_CORRECTIONS["blackbody_283K"], COLOR_CORRECTIONS["blackbody_400K"]
    return {k: round(a[k] + t * (b[k] - a[k]), 4) for k in a}


COLOR_CORRECTIONS["blackbody_300K"] = _bb_300k()

#: Spectral families used per band set (plan §4.3).
DEFAULT_SPECTRUM = {"W1": "flat_fnu", "W2": "flat_fnu",
                    "W3": "blackbody_300K", "W4": "blackbody_300K"}


# -- spectrum <-> magnitude <-> DN ---------------------------------------
def vega_mag_from_fnu(band: str, f_nu_jy: float, spectrum: str) -> float:
    """WISE Vega magnitude of a source whose spectrum is ``spectrum`` and
    whose flux density at the isophotal wavelength is ``f_nu_jy``
    (Supplement Eq. 2 inverted: F_nu = F*_nu0 10^(-m/2.5) / f_c)."""
    fc = COLOR_CORRECTIONS[spectrum][band]
    return -2.5 * math.log10(fc * f_nu_jy / F_NU0_STAR_JY[band])


def fnu_from_vega_mag(band: str, mag: float, spectrum: str) -> float:
    fc = COLOR_CORRECTIONS[spectrum][band]
    return F_NU0_STAR_JY[band] * 10 ** (-0.4 * mag) / fc


def dn_from_vega_mag(mag: float, magzp: float) -> float:
    """Total DN of a source of Vega magnitude ``mag`` in a frame of zero
    point ``magzp`` (single-exposure relative photometric ZP)."""
    return 10 ** (0.4 * (magzp - mag))


# -- PRF ------------------------------------------------------------------
@dataclass
class PRFGrid:
    band: str
    templates: dict      # (ex, ey) -> (641, 641) float32 unit-sum
    sha256: str | None = None

    @classmethod
    def load(cls, band: str, prf_dir: Path, elements=None) -> "PRFGrid":
        """Load the 9 x 9 templates of ``band`` from ``prf_dir`` (the
        unpacked ``pass2_w?_psf.tar.gz``). ``elements`` restricts to a
        subset of (ex, ey) 0-based indices."""
        import hashlib

        from astropy.io import fits

        b = band[-1]
        templates = {}
        h = hashlib.sha256()
        for ex in range(PRF_GRID):
            for ey in range(PRF_GRID):
                if elements is not None and (ex, ey) not in elements:
                    continue
                name = f"wise-w{b}-psf-wpro-09x09-{ex + 1:02d}x{ey + 1:02d}.fits"
                path = Path(prf_dir) / name
                data = np.asarray(fits.getdata(path), dtype=np.float32)
                n = PRF_SHAPE[band]
                if data.shape != (n, n):
                    raise ValueError(f"{name}: unexpected shape {data.shape}")
                h.update(path.read_bytes())
                templates[(ex, ey)] = data / data.sum()
        return cls(band=band, templates=templates, sha256=h.hexdigest())

    def element_of(self, x_frame: float, y_frame: float) -> tuple[int, int]:
        """0-based (ex, ey) of the 9 x 9 grid for a full-frame pixel
        position (element 01x01 is the lower-left corner)."""
        step = FRAME_PIX[self.band] / PRF_GRID
        ex = int(min(max(math.floor(x_frame / step), 0), PRF_GRID - 1))
        ey = int(min(max(math.floor(y_frame / step), 0), PRF_GRID - 1))
        return ex, ey

    def mean_template(self) -> np.ndarray:
        m = np.mean(np.stack(list(self.templates.values())), axis=0)
        return m / m.sum()

    def render(self, x: float, y: float, element=None, half: int = 16,
               template: np.ndarray | None = None):
        """Stamp (2*half+1)^2 of a unit-flux source at cutout pixel
        position (x, y) (0-based, pixel centres at integers).

        Returns (stamp, ox, oy): stamp[b, a] is the flux in cutout pixel
        (ox + a, oy + b). The template is bilinearly interpolated at the
        sub-pixel offsets and scaled by PRF_OVERSAMPLE^2 so that a
        template summing to 1 over its 1/8-pixel grid gives a stamp
        summing to the enclosed flux fraction in detector pixels.
        """
        t = template if template is not None else self.templates[element]
        ox, oy = int(round(x)) - half, int(round(y)) - half
        n = 2 * half + 1
        ii = ox + np.arange(n)            # cutout x of stamp columns
        jj = oy + np.arange(n)            # cutout y of stamp rows
        c = PRF_CENTRE[self.band]
        u = c + PRF_OVERSAMPLE * (ii - x)   # template column (float)
        w = c + PRF_OVERSAMPLE * (jj - y)   # template row
        return (_bilinear_grid(t, u, w) * PRF_OVERSAMPLE ** 2).astype(np.float32), ox, oy


def _bilinear_grid(t: np.ndarray, u: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Bilinear samples of ``t`` on the outer grid of columns ``u`` and
    rows ``w`` (float indices); zero outside the template."""
    nrow, ncol = t.shape
    u0 = np.floor(u).astype(int)
    w0 = np.floor(w).astype(int)
    fu = (u - u0)[None, :]
    fw = (w - w0)[:, None]
    oku = (u0 >= 0) & (u0 < ncol - 1)
    okw = (w0 >= 0) & (w0 < nrow - 1)
    u0c = np.clip(u0, 0, ncol - 2)
    w0c = np.clip(w0, 0, nrow - 2)
    a = t[w0c[:, None], u0c[None, :]]
    b = t[w0c[:, None], u0c[None, :] + 1]
    c = t[w0c[:, None] + 1, u0c[None, :]]
    d = t[w0c[:, None] + 1, u0c[None, :] + 1]
    out = (a * (1 - fu) * (1 - fw) + b * fu * (1 - fw)
           + c * (1 - fu) * fw + d * fu * fw)
    out *= (okw[:, None] & oku[None, :])
    return out


def add_stamp(img: np.ndarray, stamp: np.ndarray, ox: int, oy: int,
              amplitude: float) -> np.ndarray:
    """``img`` + amplitude * stamp placed at (ox, oy), clipped to bounds."""
    out = np.array(img, dtype=float, copy=True)
    ny, nx = img.shape
    h, w = stamp.shape
    y0, y1 = max(0, oy), min(ny, oy + h)
    x0, x1 = max(0, ox), min(nx, ox + w)
    if y1 > y0 and x1 > x0:
        out[y0:y1, x0:x1] += amplitude * stamp[y0 - oy:y1 - oy, x0 - ox:x1 - ox]
    return out


# -- exact local response of the matched filter --------------------------
@dataclass
class ResponseWindow:
    """Change of the flux map per unit injected flux, on a window."""

    delta: np.ndarray   # (h, w) d(flux) per unit source flux
    ox: int             # cutout x of delta[:, 0]
    oy: int             # cutout y of delta[0, :]

    def sample(self, x, y):
        """Bilinear sample at cutout pixel positions (0 outside)."""
        x = np.asarray(x, dtype=float) - self.ox
        y = np.asarray(y, dtype=float) - self.oy
        h, w = self.delta.shape
        out = np.zeros(np.broadcast(x, y).shape)
        ok = (x >= 0) & (x <= w - 1.001) & (y >= 0) & (y <= h - 1.001) \
            & np.isfinite(x) & np.isfinite(y)
        if ok.any():
            xs, ys = x[ok], y[ok]
            x0 = np.floor(xs).astype(int)
            y0 = np.floor(ys).astype(int)
            fx, fy = xs - x0, ys - y0
            d = self.delta
            out[ok] = (d[y0, x0] * (1 - fx) * (1 - fy) + d[y0, x0 + 1] * fx * (1 - fy)
                       + d[y0 + 1, x0] * (1 - fx) * fy + d[y0 + 1, x0 + 1] * fx * fy)
        return out


def stamp_response(fm, stamp: np.ndarray, ox: int, oy: int) -> ResponseWindow:
    """Exact flux-map change for a unit-amplitude ``stamp`` added to the
    image of ``fm`` (a FluxMap built with ``keep_inputs=True``).

    The estimator's numerator is the cross-correlation of the masked,
    background-subtracted image with the kernel; adding s*stamp adds
    s * xcorr(stamp * good, kernel) to it, and the denominator map is
    unchanged, so delta(flux) = xcorr(stamp * good, kernel) / denom.
    """
    from scipy.signal import fftconvolve

    if fm.denom is None or fm.kernel is None or fm.good is None:
        raise ValueError("FluxMap lacks stored inputs (build with keep_inputs)")
    ny, nx = fm.denom.shape
    h, w = stamp.shape
    # clip the stamp to the image
    y0, y1 = max(0, oy), min(ny, oy + h)
    x0, x1 = max(0, ox), min(nx, ox + w)
    if y1 <= y0 or x1 <= x0:
        return ResponseWindow(np.zeros((1, 1)), 0, 0)
    s = stamp[y0 - oy:y1 - oy, x0 - ox:x1 - ox] * fm.good[y0:y1, x0:x1]
    k = fm.kernel                     # already flipped for cross-correlation
    ky, kx = k.shape
    full = fftconvolve(s, k, mode="full")   # == _xcorr's FFT product, uncropped
    # _xcorr crops the 'full' product at ((ky-1)//2, (kx-1)//2)
    wy0 = y0 - (ky - 1) // 2
    wx0 = x0 - (kx - 1) // 2
    # intersect with the image
    fy0, fy1 = max(0, wy0), min(ny, wy0 + full.shape[0])
    fx0, fx1 = max(0, wx0), min(nx, wx0 + full.shape[1])
    num = full[fy0 - wy0:fy1 - wy0, fx0 - wx0:fx1 - wx0]
    den = fm.denom[fy0:fy1, fx0:fx1]
    # Positions whose PSF weight is mostly on masked pixels have tiny
    # denominators; the stack never uses them (good_frac >= 0.7) but a
    # bilinear lookup next to them would be polluted, so zero them.
    floor = 0.25 * float(np.nanmax(fm.denom))
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.where(den > floor, num / den, 0.0)
    return ResponseWindow(delta.astype(np.float64), fx0, fy0)


# -- analytic PSF renderer (ZTF / PS1, which ship no per-frame PRF) ----------
@dataclass
class MoffatPSF:
    """Circular Moffat profile I(r) ∝ (1 + (r/α)²)^(−β) with the frame's
    seeing FWHM (pixels) and β = 3 (a ground-based default), integrated
    over pixels by 4 x 4 sub-sampling. Used where the archive ships no
    empirical PRF; the chain's throughput on it is what the asteroid
    positive control measures against the archive's own photometry."""

    fwhm_pix: float
    beta: float = 3.0
    band: str = ""
    sub: int = 4

    @property
    def alpha(self) -> float:
        return self.fwhm_pix / (2.0 * math.sqrt(2.0 ** (1.0 / self.beta) - 1.0))

    def render(self, x: float, y: float, half: int = 16, element=None):
        ox, oy = int(round(x)) - half, int(round(y)) - half
        n = 2 * half + 1
        s = self.sub
        off = (np.arange(s) + 0.5) / s - 0.5
        xx = (ox + np.arange(n))[:, None] + off[None, :]          # (n, s) sub-pixel x
        yy = (oy + np.arange(n))[:, None] + off[None, :]
        dx2 = ((xx - x) ** 2).reshape(n * s)
        dy2 = ((yy - y) ** 2).reshape(n * s)
        r2 = dy2[:, None] + dx2[None, :]
        prof = (1.0 + r2 / self.alpha ** 2) ** (-self.beta)
        stamp = prof.reshape(n, s, n, s).sum(axis=(1, 3))
        # normalise to the total Moffat flux (analytic), so the stamp sum
        # is the enclosed fraction
        total = math.pi * self.alpha ** 2 / (self.beta - 1.0) * s * s
        return (stamp / total).astype(np.float32), ox, oy


def stamp_response_multiphase(fm, stamp: np.ndarray, ox: int, oy: int,
                              flux_scale: float = 1.0) -> ResponseWindow:
    """Exact flux-map change for a stamp added to the image of a FluxMap
    whose maps are interleaved over ``upsample`` x ``upsample`` sub-pixel
    kernel phases (SPHEREx): one local response per phase, interleaved
    into a window on the upsampled grid (the coordinate system of
    ``FluxMap.world2pix``). ``flux_scale`` converts the map's native
    amplitude units to the stored units (SPHEREx: MJy/sr-sum -> uJy)."""
    from scipy.signal import fftconvolve

    UP = fm.upsample
    ny, nx = fm.good.shape
    h, w = stamp.shape
    y0, y1 = max(0, oy), min(ny, oy + h)
    x0, x1 = max(0, ox), min(nx, ox + w)
    if y1 <= y0 or x1 <= x0:
        return ResponseWindow(np.zeros((1, 1)), 0, 0)
    s = stamp[y0 - oy:y1 - oy, x0 - ox:x1 - ox] * fm.good[y0:y1, x0:x1]
    wins = {}
    for dy, dx, k, den in fm.phase_parts:
        ky, kx = k.shape
        full = fftconvolve(s, k, mode="full")
        wy0 = y0 - (ky - 1) // 2; wx0 = x0 - (kx - 1) // 2
        fy0, fy1 = max(0, wy0), min(ny, wy0 + full.shape[0])
        fx0, fx1 = max(0, wx0), min(nx, wx0 + full.shape[1])
        num = full[fy0 - wy0:fy1 - wy0, fx0 - wx0:fx1 - wx0]
        d = den[fy0:fy1, fx0:fx1]
        floor = 0.25 * float(np.nanmax(den))
        with np.errstate(divide="ignore", invalid="ignore"):
            wins[(dy, dx)] = (np.where(d > floor, num / d, 0.0) * flux_scale, fx0, fy0)
    (_, fx0, fy0) = next(iter(wins.values()))
    hh, ww = next(iter(wins.values()))[0].shape
    delta = np.zeros((hh * UP, ww * UP))
    for (dy, dx), (win, _, _) in wins.items():
        delta[dy::UP, dx::UP] = win[:hh, :ww]
    return ResponseWindow(delta, fx0 * UP, fy0 * UP)


@dataclass
class OversampledPSF:
    """Renderer for a PSF plane oversampled by ``oversample`` with the
    centre at ``centre`` (SPHEREx: 101 x 101 at 0.615", centre 50)."""

    plane: np.ndarray
    oversample: int = 10
    centre: int = 50
    band: str = ""

    def render(self, x: float, y: float, half: int = 8, element=None):
        t = self.plane / self.plane.sum()
        ox, oy = int(round(x)) - half, int(round(y)) - half
        n = 2 * half + 1
        ii = ox + np.arange(n); jj = oy + np.arange(n)
        u = self.centre + self.oversample * (ii - x)
        w = self.centre + self.oversample * (jj - y)
        return (_bilinear_grid(t, u, w) * self.oversample ** 2).astype(np.float32), ox, oy

    def mean_template(self):
        return self.plane / self.plane.sum()
