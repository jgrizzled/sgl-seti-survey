"""Spectral-archive family v1 machinery (hypotheses D4/D5, thresholds §2).

load_spectrum(inst, path)  -> dict(wave_nm [vacuum, barycentric], flux, err,
                                    mask, mjd_mid, meta)   (per-order arrays
                                    stacked; merge happens on the grid)
Grid(inst)                 -> common log-lambda grid for the instrument
to_grid(spec, grid)        -> (flux, err, mask) on the grid (order merge by
                              inverse-variance weights)
normalise(flux, err, mask) -> running-median (201 px) pseudo-continuum
                              normalisation
Ensemble                   -> leave-one-out templates, per-pixel sigma,
                              matched filter S(lambda), cell maxima
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage import convolve1d
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parents[1]
CFG = json.load(open(HERE / "configs" / "threshold_freeze_v1_1.json"))
C_KMS = 299792.458
RUNMED = CFG["normalisation"]["running_median_px"]
XS_BANDS = {"XSHOOTER-UVB": [300, 560], "XSHOOTER-VIS": [530, 1020], "XSHOOTER-NIR": [990, 2480]}


def inst_cfg(inst):
    if inst in CFG["instruments"]:
        return CFG["instruments"][inst]
    if inst.startswith("XSHOOTER"):
        c = dict(CFG["instruments"]["XSHOOTER"]); c["band_nm"] = XS_BANDS[inst]
        return c
    raise KeyError(inst)


def air_to_vacuum_nm(w_nm):
    """Morton (2000) IAU standard, valid 200-2000 nm (input air nm)."""
    s2 = (1e3 / w_nm) ** 2          # (1/um)^2
    n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
    return w_nm * n


# ------------------------------------------------------------------ loaders
def load_spectrum(inst, path):
    """Return dict with per-order lists (1D products = one order)."""
    c = inst_cfg(inst)
    h = fits.open(path, memmap=False)
    if c["archive"] == "ESO-phase3":
        d = h[1].data
        hd = h[0].header
        wunit = str(d.columns["WAVE"].unit or "angstrom").lower()
        w = np.asarray(d["WAVE"][0], float) / (10.0 if wunit.startswith("a") else 1.0)   # Angstrom -> nm (X-shooter IDPs are in nm)
        fcol = c.get("flux_column", "FLUX")
        if fcol not in d.columns.names:
            fcol = "FLUX"
        ecol = fcol.replace("FLUX", "ERR")
        if ecol not in d.columns.names:
            ecol = "ERR"
        qcol = fcol.replace("FLUX", "QUAL")
        if qcol not in d.columns.names:
            qcol = "QUAL" if "QUAL" in d.columns.names else None
        f = np.asarray(d[fcol][0], float); e = np.asarray(d[ecol][0], float)
        q = np.asarray(d[qcol][0]) if qcol else np.zeros_like(f, int)
        err_synth = False
        good_e = np.isfinite(e) & (e > 0)
        if good_e.mean() < 0.5:                 # HARPS s1d: no error vector (FLUXERR=-1)
            fin = np.isfinite(f) & (f > 0)
            e = np.sqrt(np.clip(f, 0, None) + 1.0)
            snr = float(hd.get("SNR", 0) or 0)
            if snr > 0 and fin.any():
                e *= np.nanmedian(f[fin] / e[fin]) / snr     # median(f/e) = header SNR
            err_synth = True
        mask = ~np.isfinite(f) | ~np.isfinite(e) | (e <= 0) | (q != 0)
        if c.get("air_to_vacuum"):
            w = air_to_vacuum_nm(w)
        # SPECSYS=BARYCENT for all ESO 1D products used here (checked at recon)
        if hd.get("SPECSYS", "BARYCENT") != "BARYCENT":
            berv = float(hd.get("HIERARCH ESO DRS BERV", hd.get("HIERARCH ESO QC BERV", 0.0)))
            w = w * (1 + berv / C_KMS)
        mjd_mid = float(hd["MJD-OBS"]) + 0.5 * float(hd.get("EXPTIME", 0)) / 86400
        meta = dict(snr=hd.get("SNR"), exptime=hd.get("EXPTIME"), berv=hd.get("HIERARCH ESO DRS BERV", hd.get("HIERARCH ESO QC BERV")),
                    object=hd.get("OBJECT"), specsys=hd.get("SPECSYS"), fluxcal=hd.get("FLUXCAL"), flux_column=fcol, err_synth=err_synth)
        orders = [(w, f, e, mask)]
    elif inst == "SPIRou":
        hd = h[1].header
        berv = float(hd["BERV"])
        F = h["FluxAB"].data; W = h["WaveAB"].data; B = h["BlazeAB"].data
        names = [x.name for x in h]
        R = h["Recon"].data if "Recon" in names else np.ones_like(F)
        OH = h["OHLine"].data if "OHLine" in names else np.zeros_like(F)
        if "MJDMID" in hd:
            mjd_mid = float(hd["MJDMID"])
        else:                                   # older APERO headers: MJD-OBS + half exposure (primary header)
            hp = h[0].header
            mjd_mid = float(hp.get("MJD-OBS", hd.get("MJD-OBS"))) + 0.5 * float(hp.get("EXPTIME", hd.get("EXPTIME", 0))) / 86400
        orders = []
        for o in range(F.shape[0]):
            w = np.asarray(W[o], float) * (1 + berv / C_KMS)      # observer -> barycentric
            f = np.asarray(F[o], float); b = np.asarray(B[o], float)
            fn = f / b                                            # deblazed
            # photon-ish error from the blaze-scaled flux (e-/px): sigma = sqrt(max(f,0)+ron^2)/b
            e = np.sqrt(np.clip(f, 0, None) + 30.0 ** 2) / b
            oh = np.asarray(OH[o], float)
            ohp = np.nanpercentile(oh, CFG["mask"]["spirou_ohline_pct"]) if np.isfinite(oh).any() else np.inf
            mask = (~np.isfinite(fn) | ~np.isfinite(e) | (b <= 0) | (np.asarray(R[o], float) < CFG["mask"]["spirou_recon_min"]) | (oh > ohp))
            orders.append((w, fn, e, mask))
        meta = dict(snr=hd.get("SPEMSNR"), exptime=hd.get("EXPTIME"), berv=berv, object=hd.get("OBJECT"), specsys="observer+BERV")
    elif inst == "CARMENES-VIS":
        hd = h[0].header
        berv = float(hd["HIERARCH CARACAL BERV"])
        S = h["SPEC"].data; Cn = h["CONT"].data; Sg = h["SIG"].data; W = h["WAVE"].data
        mjd_mid = float(hd["HIERARCH CARACAL MJD-OBS"]) + 0.5 * float(hd["EXPTIME"]) / 86400
        orders = []
        for o in range(S.shape[0]):
            w = np.asarray(W[o], float) / 10.0 * (1 + berv / C_KMS)
            f = np.asarray(S[o], float); e = np.asarray(Sg[o], float)
            mask = ~np.isfinite(f) | ~np.isfinite(e) | (e <= 0)
            orders.append((w, f, e, mask))
        meta = dict(snr=hd.get("HIERARCH CARACAL FOX SNR 50"), exptime=hd.get("EXPTIME"), berv=berv, object=hd.get("OBJECT"), specsys="observer+BERV")
    else:
        raise KeyError(inst)
    h.close()
    return dict(inst=inst, orders=orders, mjd_mid=mjd_mid, meta=meta)


# --------------------------------------------------------------------- grid
class Grid:
    def __init__(self, inst):
        c = inst_cfg(inst)
        self.inst = inst
        self.R = c["R"]; self.fwhm_px = c["fwhm_px"]
        lo, hi = c["band_nm"]
        self.step = 1.0 / (3.0 * self.R)
        n = int(np.ceil(np.log(hi / lo) / self.step))
        self.lnw = np.log(lo) + self.step * np.arange(n)
        self.w = np.exp(self.lnw)
        self.n = n

    def idx(self, w_nm):
        return np.clip(np.round((np.log(w_nm) - self.lnw[0]) / self.step).astype(int), 0, self.n - 1)

    def cell_slice(self, cell):
        iv = CFG["cells_nm"][cell]
        if iv is None:
            return slice(3, self.n - 3)
        return slice(int(self.idx(iv[0])), int(self.idx(iv[1])) + 1)


def to_grid(spec, grid):
    """Inverse-variance merge of orders onto the grid by linear interpolation."""
    num = np.zeros(grid.n); den = np.zeros(grid.n); cover = np.zeros(grid.n, bool)
    for (w, f, e, mask) in spec["orders"]:
        ok = ~mask & np.isfinite(w)
        if ok.sum() < 10:
            continue
        # order may be descending
        order = np.argsort(w)
        w_s, f_s, e_s, ok_s = w[order], f[order], e[order], ok[order]
        lo, hi = grid.idx(w_s[0]), grid.idx(w_s[-1])
        gw = grid.w[lo:hi + 1]
        fi = np.interp(gw, w_s[ok_s], f_s[ok_s])
        ei = np.interp(gw, w_s[ok_s], e_s[ok_s])
        # mark grid pixels whose nearest source pixel is masked (interp over gaps > 2 native px)
        src_ok_pos = w_s[ok_s]
        j = np.searchsorted(src_ok_pos, gw)
        j = np.clip(j, 1, len(src_ok_pos) - 1)
        gap = src_ok_pos[j] - src_ok_pos[j - 1]
        native = np.median(np.diff(w_s))
        good = gap < 2.5 * abs(native)
        wgt = np.where(good & (ei > 0), 1.0 / ei ** 2, 0.0)
        num[lo:hi + 1] += wgt * fi; den[lo:hi + 1] += wgt
        cover[lo:hi + 1] |= good
    flux = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    err = np.where(den > 0, 1 / np.sqrt(np.where(den > 0, den, 1)), np.nan)
    mask = ~cover | ~np.isfinite(flux)
    # +-3 px around masked runs (order gaps, defects)
    m = mask.copy()
    for s in (1, 2, 3):
        m[s:] |= mask[:-s]; m[:-s] |= mask[s:]
    return flux, err, m


def normalise(flux, err, mask):
    """Running-median pseudo-continuum normalisation. v1.1 A2 (per pixel, per
    spectrum): pixels whose local continuum has SNR = cont / running-median(err)
    below CFG pixel_snr_min are masked (zero-flux blue orders of M dwarfs,
    dead regions) -- the frozen SNR >= 5 gate at pixel granularity."""
    f = np.where(mask, np.nan, flux)
    e = np.where(mask, np.nan, err)
    cont = pd.Series(f).rolling(RUNMED, center=True, min_periods=RUNMED // 4).median().to_numpy()
    emed = pd.Series(e).rolling(RUNMED, center=True, min_periods=RUNMED // 4).median().to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        nf = f / cont; ne = err / cont
        lowsnr = ~(cont / emed >= CFG.get("pixel_snr_min", 5.0))
    bad = mask | ~np.isfinite(nf) | ~np.isfinite(ne) | (cont <= 0) | lowsnr
    return np.where(bad, np.nan, nf), np.where(bad, np.nan, ne), bad


def _gauss(x, a, mu, sig, c):
    return a * np.exp(-0.5 * ((x - mu) / sig) ** 2) + c


def gauss_fit(z, peak, fwhm_px, half=8):
    lo, hi = max(0, peak - half), min(len(z), peak + half + 1)
    x = np.arange(lo, hi); y = z[lo:hi]
    ok = np.isfinite(y)
    if ok.sum() < 6 or not np.isfinite(z[peak]):
        return None
    try:
        p, _ = curve_fit(_gauss, x[ok], y[ok], p0=[z[peak], peak, fwhm_px / 2.3548, 0.0], maxfev=2000)
    except Exception:
        return None
    if not np.all(np.isfinite(p)):
        return None
    return {"amp_z": float(p[0]), "mu_px": float(p[1]), "fwhm_px": float(abs(p[2]) * 2.3548), "fwhm_ratio": float(abs(p[2]) * 2.3548 / fwhm_px)}


def xcorr_shift(x, y, maxlag=60):
    x = x - np.nanmean(x); y = y - np.nanmean(y)
    x = np.where(np.isfinite(x), x, 0); y = np.where(np.isfinite(y), y, 0)
    lags = list(range(-maxlag, maxlag + 1))
    cc = [np.dot(x[maxlag:-maxlag], np.roll(y, l)[maxlag:-maxlag]) for l in lags]
    return int(lags[int(np.argmax(cc))])


STELLAR_REGIONS_NM = {"CARMENES-VIS": [(600, 620), (700, 720)], "SPIRou": [(1000, 1030), (1580, 1620)],
                      "NIRPS": [(1000, 1030), (1580, 1620)], "HARPS": [(600, 620), (640, 660)],
                      "ESPRESSO": [(600, 620), (700, 720)], "XSHOOTER-NIR": [(1000, 1030), (1580, 1620)],
                      "XSHOOTER-VIS": [(600, 620), (700, 720)], "XSHOOTER-UVB": [(430, 450), (500, 520)]}


def gauss_kernel(fwhm_px):
    sig = fwhm_px / 2.3548
    half = int(np.ceil(3 * sig))
    x = np.arange(-half, half + 1)
    k = np.exp(-0.5 * (x / sig) ** 2)
    return k / k.sum()


def stellar_line_mask(grid, rv_km_s, tol_km_s=None):
    """v1.1: grid pixels within +-tol of any frozen stellar emission line (star frame)."""
    tol = CFG["stellar_line_tolerance_km_s"] if tol_km_s is None else tol_km_s
    m = np.zeros(grid.n, bool)
    for l0 in CFG["stellar_lines_nm_vacuum"].values():
        lc = l0 * (1 + rv_km_s / C_KMS)
        lo, hi = lc * (1 - tol / C_KMS), lc * (1 + tol / C_KMS)
        if hi < grid.w[0] or lo > grid.w[-1]:
            continue
        m[int(grid.idx(lo)):int(grid.idx(hi)) + 1] = True
    return m


def despike(z, thresh=4.0, ratio=0.4):
    """v1.1: NaN isolated single-pixel spikes (|z_i| > thresh with both
    neighbours below ratio*|z_i|) -- a PSF-consistency precondition; a real
    line of FWHM >= 2.5 px keeps neighbours >= 0.5 of the peak."""
    z = z.copy()
    a = np.abs(z)
    left = np.roll(a, 1); right = np.roll(a, -1)
    left[0] = 0; right[-1] = 0
    spike = (a > thresh) & (np.nan_to_num(left) < ratio * a) & (np.nan_to_num(right) < ratio * a)
    z[spike] = np.nan
    return z, int(spike.sum())


# ----------------------------------------------------------------- ensemble
class Ensemble:
    """Null ensemble of normalised spectra on the grid (N x G)."""

    def __init__(self, grid, norm_flux, norm_err, masks, extra_mask=None):
        self.grid = grid
        self.F = np.asarray(norm_flux)            # NaN where masked
        if extra_mask is not None:
            self.F[:, extra_mask] = np.nan
        self.extra_mask = extra_mask if extra_mask is not None else np.zeros(grid.n, bool)
        self.E = np.asarray(norm_err)
        self.M = np.asarray(masks)
        self.N = self.F.shape[0]
        self.k = gauss_kernel(grid.fwhm_px)
        self.template = np.nanmedian(self.F, axis=0)
        # leave-one-out templates and residual scatter
        loo = np.empty_like(self.F)
        for i in range(self.N):
            loo[i] = np.nanmedian(np.delete(self.F, i, axis=0), axis=0)
        self.loo = loo
        res = self.F - loo
        mad = np.nanmedian(np.abs(res - np.nanmedian(res, axis=0)), axis=0) * 1.4826
        self.mad = np.where(np.isfinite(mad), mad, np.nan)
        phot = np.nanmedian(self.E, axis=0)
        sig = np.fmax(mad, phot)
        sig = np.where(np.isfinite(sig) & (sig > 0), sig, np.nan)
        self.sigma = sig
        # pixel usable if defined in >= 80% of nulls, sigma finite, and (v1.1 A2)
        # ensemble-median photon SNR per pixel >= 5
        with np.errstate(invalid="ignore", divide="ignore"):
            self.snr_pix = np.nanmedian(self.F / self.E, axis=0)
        self.pix_ok = (np.isfinite(self.F).sum(axis=0) >= 0.8 * self.N) & np.isfinite(sig) & (self.snr_pix >= CFG.get("pixel_snr_min", 5.0))
        self.mask_frac = 1 - np.isfinite(self.F).mean(axis=0)

    def matched(self, z):
        """Matched-filter S over the grid from a z array (NaN = masked)."""
        zz = np.where(np.isfinite(z), z, 0.0)
        wgt = np.isfinite(z).astype(float)
        num = convolve1d(zz, self.k, mode="constant")
        den = np.sqrt(convolve1d(wgt, self.k ** 2, mode="constant"))
        with np.errstate(invalid="ignore", divide="ignore"):
            S = num / den
        # require >= 60% of the kernel weight present
        frac = convolve1d(wgt, self.k, mode="constant")
        S[frac < 0.6] = np.nan
        return S

    def sigma_for(self, e):
        """Per-pixel scale for one spectrum: max(ensemble MAD scatter, the
        spectrum's own propagated photon error) -- thresholds section 2."""
        if e is None:
            return self.sigma
        return np.fmax(self.mad, np.where(np.isfinite(e), e, np.nan))

    def z_single(self, f, template=None, e=None):
        t = self.template if template is None else template
        z = (f - t) / self.sigma_for(e)
        z[~self.pix_ok | self.extra_mask] = np.nan
        z, self.last_spikes = despike(z)
        return z

    def z_mean(self, F_list, E_list=None):
        Fm = np.asarray(F_list)
        n_eff = np.isfinite(Fm).sum(axis=0)
        mean = np.nanmean(Fm, axis=0)
        if E_list is not None:
            erms = np.sqrt(np.nanmean(np.asarray(E_list) ** 2, axis=0))
            sig = np.fmax(self.mad, erms)
        else:
            sig = self.sigma
        z = (mean - self.template) / (sig / np.sqrt(np.fmax(n_eff, 1)))
        z[(n_eff < max(1, int(0.5 * len(F_list)))) | ~self.pix_ok | self.extra_mask] = np.nan
        z, self.last_spikes = despike(z)
        return z

    def cell_max(self, S, cell, z=None):
        """v1.1 A3: PSF-consistent cell maximum. If z is given, the maximum is
        taken over local maxima of S whose Gaussian fit to z has FWHM within
        CFG shape bounds of the instrumental FWHM (>= 6 finite px); peaks that
        fail are skipped. Falls back to the plain maximum below the peak floor."""
        sl = self.grid.cell_slice(cell)
        seg = S[sl]
        if not np.isfinite(seg).any():
            return np.nan, None
        if z is None:
            i = int(np.nanargmax(seg))
            return float(seg[i]), int(sl.start + i)
        floor = CFG.get("peak_floor_S", 4.0)
        segz = np.where(np.isfinite(seg), seg, -np.inf)
        pk, _ = find_peaks(segz, height=floor, distance=3)
        pk = pk[np.argsort(segz[pk])[::-1]][:CFG.get("max_peaks_fit", 200)]
        lo_r, hi_r = CFG["shape_fwhm_bounds"]
        for p in pk:
            g = int(sl.start + p)
            fit = gauss_fit(z, g, self.grid.fwhm_px)
            if fit is not None and lo_r <= fit["fwhm_ratio"] <= hi_r and abs(fit["mu_px"] - g) <= 2.0:
                return float(seg[p]), g
        # no compliant peak: report the largest sub-floor value (or floor if all peaks failed)
        below = np.where(segz < floor, segz, -np.inf)
        if np.isfinite(below).any() and below.max() > -np.inf:
            i = int(np.argmax(below))
            return float(min(below[i], floor)), int(sl.start + i)
        return float(floor), None

    def null_stats(self, cells):
        """Per null spectrum (LOO template) the PSF-consistent cell maxima."""
        out = {c: np.full(self.N, np.nan) for c in cells}
        peaks = {c: np.full(self.N, -1, int) for c in cells}
        for i in range(self.N):
            z = self.z_single(self.F[i], self.loo[i], self.E[i])
            S = self.matched(z)
            for c in cells:
                out[c][i], p = self.cell_max(S, c, z)
                peaks[c][i] = -1 if p is None else p
        return out, peaks

    def wavelength_shift_px(self, f, regions_nm):
        """v1.1 A1: stellar-region cross-correlation shift of a normalised spectrum
        against the ensemble template (max |shift| over regions, px)."""
        shifts = []
        for lo, hi in regions_nm:
            if hi < self.grid.w[0] or lo > self.grid.w[-1]:
                continue
            sl = slice(int(self.grid.idx(lo)), int(self.grid.idx(hi)))
            x = self.template[sl]; y = f[sl]
            if np.isfinite(x).mean() < 0.5 or np.isfinite(y).mean() < 0.5:
                continue
            shifts.append(xcorr_shift(x, y, 60))
        return (max(shifts, key=abs) if shifts else 0), shifts
