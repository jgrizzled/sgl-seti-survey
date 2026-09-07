"""Shared machinery for the high-energy (Fermi-LAT) crossings survey
(hypotheses.md freeze v1.0). Everything data-facing lives here:

  spacecraft table (weekly files -> cached npz), per-position interval
  gate (theta <= 60, DATA_QUAL > 0, LAT_CONFIG = 1, Moon/Sun > 8 deg),
  CALDB effective area per lane and PSF containment / sampling,
  photon loading with the shared gate, windows and pseudo-windows,
  the five frozen statistics (S_stack_L/H, S_event_L/H, S_burst).

Nothing here decides anything: units, thresholds and splits are read
from the freeze (hypotheses.md D1-D11) via the FREEZE dict below.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body, get_sun
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from scipy.stats import poisson

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import recon_scan as R  # noqa: E402  (events, windows, positions, sep_deg)

RUN = REPO / "runs" / "highenergy-crossings" / "v1"
SC_DIR = REPO / "runs" / "highenergy-crossings" / "recon" / "fermi_spacecraft"
CALDB = RUN / "caldb"
PHOT = RUN / "photons"
RES = Path(__file__).resolve().parents[1] / "results"

MET_MJD0 = 51910.0 + 7.428703703703703e-4     # 2001-01-01T00:00:00 UTC; MET is in TT-like SI s
RSUN_AU = 0.00465047
KM_PER_AU = 1.495978707e8
CM_PER_AU = 1.495978707e13
MEV_ERG = 1.602176634e-6

FREEZE = {
    "era_utc": ("2008-08-04T15:43:36", "2026-09-07T00:13:13"),
    "lanes": {"L": {"emin": 100.0, "emax": 1000.0, "r_ap": 2.0},
              "H": {"emin": 1000.0, "emax": 300000.0, "r_ap": 0.7},
              "U": {"emin": 100.0, "emax": 300000.0, "r_ap": 1.0}},
    "gamma": 2.0,
    "theta_max": 60.0, "theta_strict": 50.0,
    "zenith_max": 100.0, "zenith_strict": 90.0,
    "moon_sun_excl": 8.0, "moon_sun_strict": 16.0,
    "window_livetime_min_s": 1000.0,
    "pw_span_d": {"1.2Rsun": 60.0, "2.5Rsun": 60.0, "0.1AU": 120.0},
    "pw_excl_factor": 1.5, "pw_min_expo_frac": 0.5,
    "boxcar_s": 1000.0, "boxcar_step_s": 100.0, "boxcar_min_live_s": 100.0,
    "dispersion_gate": (0.05, 0.10),
    "n_trials": 145, "fwer": 0.05,   # amendment v1.1: 29 confirmatory units x 5
    "seed": 20260907,
    "dev_targets": ["gj-908"],
    "forced_dev": [("B", "van-maanen", "2020-04-02")],
    # amendment v1.1: unit demoted after its statistics were printed during
    # machinery testing before the threshold freeze
    "forced_dev_units": [("A", "gj-1276", "1.2Rsun")],
    "rungs": ["1.2Rsun", "2.5Rsun", "0.1AU"],
}
FREEZE["alpha_trial"] = 1.0 - (1.0 - FREEZE["fwer"]) ** (1.0 / FREEZE["n_trials"])
FREEZE["T"] = float(-np.log10(FREEZE["alpha_trial"]))
RUNG_R_AU = {"1.2Rsun": 1.2 * RSUN_AU, "2.5Rsun": 2.5 * RSUN_AU, "0.1AU": 0.1}


def met_to_mjd(met):
    return MET_MJD0 + np.asarray(met) / 86400.0


def mjd_to_met(mjd):
    return (np.asarray(mjd) - MET_MJD0) * 86400.0


def sep_deg_vec(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = (np.radians(np.asarray(x, dtype=float)) for x in (ra1, de1, ra2, de2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


# ------------------------------------------------------------ spacecraft
_SC = None


def spacecraft():
    """All weekly spacecraft files concatenated (cached npz)."""
    global _SC
    if _SC is not None:
        return _SC
    cache = RUN / "spacecraft_all.npz"
    if cache.exists():
        z = np.load(cache)
        _SC = {k: z[k] for k in z.files}
        return _SC
    cols = {"START": [], "STOP": [], "RA_SCZ": [], "DEC_SCZ": [], "RA_ZENITH": [], "DEC_ZENITH": [],
            "LIVETIME": [], "DATA_QUAL": [], "LAT_CONFIG": [], "LAT_MODE": [], "ROCK_ANGLE": []}
    files = sorted(SC_DIR.glob("lat_spacecraft_weekly_w*_p310_v001.fits"))
    for f in files:
        with fits.open(f) as h:
            d = h["SC_DATA"].data
            for k in cols:
                cols[k].append(np.asarray(d[k]))
    sc = {k: np.concatenate(v) for k, v in cols.items()}
    o = np.argsort(sc["START"])
    sc = {k: v[o] for k, v in sc.items()}
    # drop exact duplicates at week boundaries
    keep = np.r_[True, np.diff(sc["START"]) > 0]
    sc = {k: v[keep] for k, v in sc.items()}
    for k in ("RA_SCZ", "DEC_SCZ", "RA_ZENITH", "DEC_ZENITH", "LIVETIME", "ROCK_ANGLE"):
        sc[k] = sc[k].astype(np.float32)
    for k in ("DATA_QUAL", "LAT_CONFIG", "LAT_MODE"):
        sc[k] = sc[k].astype(np.int16)
    sc["n_files"] = np.array([len(files)])
    np.savez(cache, **sc)
    _SC = sc
    return sc


# ------------------------------------------------------------ Moon / Sun
_EPH = None


def ephem():
    """Geocentric Moon and Sun unit vectors on a 1-h grid over the era (cached)."""
    global _EPH
    if _EPH is not None:
        return _EPH
    cache = RUN / "moon_sun_1h.npz"
    if cache.exists():
        z = np.load(cache)
        _EPH = {k: z[k] for k in z.files}
        return _EPH
    t0 = Time(FREEZE["era_utc"][0]).mjd - 200.0
    t1 = Time(FREEZE["era_utc"][1]).mjd + 200.0
    mjd = np.arange(t0, t1, 1.0 / 24.0)
    t = Time(mjd, format="mjd", scale="utc")
    moon = get_body("moon", t)
    sun = get_sun(t)

    def uv(c):
        ra = c.ra.rad; de = c.dec.rad
        return np.stack([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)], axis=1).astype(np.float32)
    _EPH = {"mjd": mjd, "moon": uv(moon), "sun": uv(sun)}
    np.savez(cache, **_EPH)
    return _EPH


def body_sep(mjd, ra, dec, body):
    e = ephem()
    i = np.clip(np.searchsorted(e["mjd"], mjd) - 1, 0, len(e["mjd"]) - 2)
    f = ((mjd - e["mjd"][i]) / (e["mjd"][i + 1] - e["mjd"][i]))[:, None]
    v = e[body][i] * (1 - f) + e[body][i + 1] * f
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    ra, de = np.radians(ra), np.radians(dec)
    p = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])
    return np.degrees(np.arccos(np.clip(v @ p, -1, 1)))


# ------------------------------------------------------------ IRF
_AEFF = None
_PSF = None


def aeff_tables():
    global _AEFF
    if _AEFF is None:
        h = fits.open(CALDB / "aeff_P8R3_SOURCE_V2_FB.fits")
        out = {}
        for side in ("FRONT", "BACK"):
            d = h[f"EFFECTIVE AREA_{side}"].data
            E_lo = d["ENERG_LO"][0]; E_hi = d["ENERG_HI"][0]
            ct_lo = d["CTHETA_LO"][0]; ct_hi = d["CTHETA_HI"][0]
            A = d["EFFAREA"][0].reshape(len(ct_lo), len(E_lo))
            out[side] = (E_lo, E_hi, ct_lo, ct_hi, A)
        _AEFF = out
    return _AEFF


def aeff_lane(lane):
    """Gamma-weighted front+back effective area (m^2) vs cos(theta) for a lane;
    returns (ct_mid, A_lane)."""
    tabs = aeff_tables()
    E_lo, E_hi, ct_lo, ct_hi, _ = tabs["FRONT"]
    emin, emax = FREEZE["lanes"][lane]["emin"], FREEZE["lanes"][lane]["emax"]
    g = FREEZE["gamma"]
    lo = np.maximum(E_lo, emin); hi = np.minimum(E_hi, emax)
    m = hi > lo
    # photon weight per bin for a power law: int E^-g dE
    w = np.zeros_like(E_lo)
    w[m] = (lo[m] ** (1 - g) - hi[m] ** (1 - g)) / (g - 1)
    A = tabs["FRONT"][4] + tabs["BACK"][4]
    A_l = (A * w[None, :]).sum(axis=1) / w.sum()
    return 0.5 * (ct_lo + ct_hi), A_l


def aeff_interp(lane, costheta):
    ct, A = aeff_lane(lane)
    return np.interp(costheta, ct, A, left=0.0, right=A[-1])


def lane_mean_energy_mev(lane):
    emin, emax = FREEZE["lanes"][lane]["emin"], FREEZE["lanes"][lane]["emax"]
    return np.log(emax / emin) / (1.0 / emin - 1.0 / emax)   # Gamma = 2


def psf_tables():
    global _PSF
    if _PSF is None:
        h = fits.open(CALDB / "psf_P8R3_SOURCE_V2_FB.fits")
        out = {}
        for side in ("FRONT", "BACK"):
            d = h[f"RPSF_{side}"].data
            E_lo = d["ENERG_LO"][0]; E_hi = d["ENERG_HI"][0]
            ct_lo = d["CTHETA_LO"][0]; ct_hi = d["CTHETA_HI"][0]
            shape = (len(ct_lo), len(E_lo))
            pars = {k: d[k][0].reshape(shape) for k in ("NCORE", "NTAIL", "SCORE", "STAIL", "GCORE", "GTAIL")}
            sc = h[f"PSF_SCALING_PARAMS_{side}"].data["PSFSCALE"][0]
            out[side] = (E_lo, E_hi, ct_lo, ct_hi, pars, sc)
        _PSF = out
    return _PSF


def _king(x, sigma, gamma):
    return (1.0 / (2 * np.pi * sigma ** 2)) * (1 - 1 / gamma) * (1 + x ** 2 / (2 * gamma * sigma ** 2)) ** (-gamma)


def psf_radial_cdf(E_mev, costheta, side, r_deg_grid):
    """Cumulative containment vs radius (deg) for one energy / inclination."""
    E_lo, E_hi, ct_lo, ct_hi, p, sc = psf_tables()[side]
    ie = int(np.clip(np.searchsorted(E_hi, E_mev), 0, len(E_lo) - 1))
    ic = int(np.clip(np.searchsorted(ct_hi, costheta), 0, len(ct_lo) - 1))
    c0, c1, beta = sc[0], sc[1], sc[2]          # PSFSCALE = (c0, c1, beta) in radians, beta = -0.8 as stored
    Sp = np.degrees(np.sqrt((c0 * (E_mev / 100.0) ** beta) ** 2 + c1 ** 2))   # deg
    xs = np.linspace(0, max(r_deg_grid.max(), 90.0) / Sp, 20001)
    pdf = p["NCORE"][ic, ie] * _king(xs, p["SCORE"][ic, ie], p["GCORE"][ic, ie]) \
        + p["NTAIL"][ic, ie] * p["NCORE"][ic, ie] * _king(xs, p["STAIL"][ic, ie], p["GTAIL"][ic, ie])
    cdf = np.cumsum(2 * np.pi * xs * pdf) * (xs[1] - xs[0])
    cdf /= cdf[-1]
    return np.interp(r_deg_grid / Sp, xs, cdf)


def lane_containment(lane, costheta=0.8, n_e=12):
    """Gamma-weighted, front+back-averaged containment fraction inside the
    lane aperture, at a representative inclination."""
    emin, emax, r = (FREEZE["lanes"][lane][k] for k in ("emin", "emax", "r_ap"))
    g = FREEZE["gamma"]
    edges = np.geomspace(emin, emax, n_e + 1)
    w = (edges[:-1] ** (1 - g) - edges[1:] ** (1 - g)) / (g - 1)
    Es = np.sqrt(edges[:-1] * edges[1:])
    tot = 0.0
    for E, wi in zip(Es, w):
        # weight front/back by their effective areas at this energy
        tabs = aeff_tables()
        E_lo, E_hi, ct_lo, ct_hi, Af = tabs["FRONT"]
        Ab = tabs["BACK"][4]
        ie = int(np.clip(np.searchsorted(E_hi, E), 0, len(E_lo) - 1))
        ic = int(np.clip(np.searchsorted(ct_hi, costheta), 0, len(ct_lo) - 1))
        af, ab = Af[ic, ie], Ab[ic, ie]
        cf = psf_radial_cdf(E, costheta, "FRONT", np.array([r]))[0]
        cb = psf_radial_cdf(E, costheta, "BACK", np.array([r]))[0]
        tot += wi * ((af * cf + ab * cb) / max(af + ab, 1e-9))
    return float(tot / w.sum())


def psf_sample(rng, E_mev, costheta, side):
    """Draw one PSF offset (deg) by inverse CDF."""
    r_grid = np.geomspace(1e-3, 60.0, 400)
    cdf = psf_radial_cdf(E_mev, costheta, side, r_grid)
    u = rng.uniform(0, cdf[-1])
    return float(np.interp(u, cdf, r_grid))


# ------------------------------------------------------------ interval gate
def intervals(ra, dec, mjd0, mjd1, strict=False):
    """Spacecraft intervals overlapping [mjd0, mjd1] with per-interval gate
    and per-lane exposure (m^2 s), pro-rated to the overlap."""
    sc = spacecraft()
    met0, met1 = mjd_to_met(mjd0), mjd_to_met(mjd1)
    i0 = np.searchsorted(sc["STOP"], met0)
    i1 = np.searchsorted(sc["START"], met1)
    sl = slice(i0, i1)
    st = sc["START"][sl].astype(float); sp = sc["STOP"][sl].astype(float)
    if len(st) == 0:
        return None
    frac = np.clip((np.minimum(sp, met1) - np.maximum(st, met0)) / np.maximum(sp - st, 1e-9), 0, 1)
    theta = sep_deg_vec(sc["RA_SCZ"][sl], sc["DEC_SCZ"][sl], ra, dec)
    zen = sep_deg_vec(sc["RA_ZENITH"][sl], sc["DEC_ZENITH"][sl], ra, dec)
    mid = met_to_mjd(0.5 * (st + sp))
    moon = body_sep(mid, ra, dec, "moon")
    sun = body_sep(mid, ra, dec, "sun")
    th_max = FREEZE["theta_strict"] if strict else FREEZE["theta_max"]
    ex = FREEZE["moon_sun_strict"] if strict else FREEZE["moon_sun_excl"]
    gate = ((theta <= th_max) & (sc["DATA_QUAL"][sl] > 0) & (sc["LAT_CONFIG"][sl] == 1)
            & (moon > ex) & (sun > ex))
    live = sc["LIVETIME"][sl].astype(float) * frac
    ct = np.cos(np.radians(theta))
    out = {"start": st, "stop": sp, "frac": frac, "theta": theta, "zenith": zen, "moon": moon, "sun": sun,
           "gate": gate, "live": np.where(gate, live, 0.0), "live_raw": live, "costheta": ct,
           "mode": sc["LAT_MODE"][sl], "rock": sc["ROCK_ANGLE"][sl]}
    for lane in FREEZE["lanes"]:
        out["expo_" + lane] = np.where(gate, live * aeff_interp(lane, ct), 0.0)
    out["moon_excluded_s"] = float((live * ((moon <= ex) | (sun <= ex))).sum())
    return out


# ------------------------------------------------------------ photons
_PH = {}


def photons(key, ra, dec):
    """Gated photons for a position key (A-<target>, B-<target>, ctrl-*):
    SOURCE class, front+back, zenith <= 100, inside a gated interval.
    Returns dict of arrays sorted by time (MET), with sep (deg) to the
    channel position and the interval index."""
    if key in _PH:
        return _PH[key]
    cache = PHOT / key / "gated.npz"
    if cache.exists():
        z = np.load(cache)
        _PH[key] = {k: z[k] for k in z.files}
        return _PH[key]
    man = json.loads((PHOT / key / "manifest.json").read_text())
    cols = {"TIME": [], "RA": [], "DEC": [], "ENERGY": [], "ZENITH_ANGLE": [], "EVENT_CLASS": [], "EVENT_TYPE": []}
    for f in man["files"]:
        if "_PH" not in f["file"]:
            continue
        with fits.open(PHOT / key / f["file"]) as h:
            d = h["EVENTS"].data
            for k in cols:
                v = d[k]
                if k in ("EVENT_CLASS", "EVENT_TYPE"):
                    # bit arrays (32X) -> integer bitmask
                    v = np.asarray(v)
                    if v.ndim == 2:
                        v = (v[:, ::-1] * (1 << np.arange(v.shape[1]))).sum(axis=1)
                cols[k].append(np.asarray(v))
    ph = {k: np.concatenate(v) for k, v in cols.items()}
    o = np.argsort(ph["TIME"])
    ph = {k: v[o] for k, v in ph.items()}
    ok = ((ph["EVENT_CLASS"].astype(np.int64) & 128) > 0) & ((ph["EVENT_TYPE"].astype(np.int64) & 3) > 0) \
        & (ph["ZENITH_ANGLE"] <= FREEZE["zenith_max"])
    ph = {k: v[ok] for k, v in ph.items()}
    sc = spacecraft()
    idx = np.searchsorted(sc["START"], ph["TIME"], side="right") - 1
    idx = np.clip(idx, 0, len(sc["START"]) - 1)
    inside = ph["TIME"] < sc["STOP"][idx]
    theta = sep_deg_vec(sc["RA_SCZ"][idx], sc["DEC_SCZ"][idx], ra, dec)
    mid = met_to_mjd(0.5 * (sc["START"][idx] + sc["STOP"][idx]))
    moon = body_sep(mid, ra, dec, "moon"); sun = body_sep(mid, ra, dec, "sun")
    gate = inside & (theta <= FREEZE["theta_max"]) & (sc["DATA_QUAL"][idx] > 0) & (sc["LAT_CONFIG"][idx] == 1) \
        & (moon > FREEZE["moon_sun_excl"]) & (sun > FREEZE["moon_sun_excl"])
    out = {k: v[gate] for k, v in ph.items()}
    out["sep"] = sep_deg_vec(out["RA"], out["DEC"], ra, dec)
    out["theta"] = theta[gate]
    out["moon"] = moon[gate]
    out["mjd"] = met_to_mjd(out["TIME"])
    out["n_raw"] = np.array([len(o)]); out["n_class"] = np.array([int(ok.sum())]); out["n_gated"] = np.array([len(out["TIME"])])
    np.savez(cache, **out)
    _PH[key] = out
    return out


def count(ph, mjd0, mjd1, lane):
    lo = np.searchsorted(ph["mjd"], mjd0); hi = np.searchsorted(ph["mjd"], mjd1)
    if hi <= lo:
        return 0, np.zeros(0)
    e = ph["ENERGY"][lo:hi]; s = ph["sep"][lo:hi]
    L = FREEZE["lanes"][lane]
    m = (e >= L["emin"]) & (e < L["emax"]) & (s <= L["r_ap"])
    return int(m.sum()), ph["mjd"][lo:hi][m]


# ------------------------------------------------------------ windows
def unit_windows(ch, tid, rung, events=None):
    """Searchable-candidate windows of a unit in the era: list of dicts with
    event_id, t_ca_mjd, mjd0, mjd1, b_rsun."""
    events = events or R.load_events()
    tab = events[ch]
    sub = tab[tab["target_id"] == tid]
    era = [Time(x, format="isot", scale="utc").mjd for x in FREEZE["era_utc"]]
    out = []
    for ev in sub:
        w = R.windows_for(ev)[rung]
        if w is None:
            continue
        if w[0] < era[0] or w[1] > era[1]:
            continue
        out.append({"event_id": str(ev["event_id"]), "t_ca_utc": str(ev["t_ca_utc"]), "t_ca_mjd": float(ev["mjd"]),
                    "mjd0": float(w[0]), "mjd1": float(w[1]), "dur_d": float(2 * w[2]),
                    "b_rsun": float(ev["b_min_au"]) / RSUN_AU})
    return out


def pseudo_offsets(rung, dur_d):
    D = max(dur_d, 1.0)
    span = FREEZE["pw_span_d"][rung]
    js = []
    j = 1
    while j * D <= span:
        if j * D >= FREEZE["pw_excl_factor"] * D:
            js += [j, -j]
        j += 1
    return sorted(js, key=abs)


def is_forced_dev(ch, tid, t_ca_utc):
    return any(ch == c and tid == t and t_ca_utc.startswith(d) for c, t, d in FREEZE["forced_dev"])


# ------------------------------------------------------------ statistics
def nlog10_sf(n, lam):
    """-log10 P(N >= n | lam), Poisson upper tail; 0 when n == 0."""
    n = int(n)
    if n <= 0:
        return 0.0
    p = poisson.sf(n - 1, lam)
    return float(-np.log10(max(p, 1e-300)))


def burst_stat(ph, ra, dec, mjd0, mjd1, rate_U, iv=None):
    """S_burst: max over sliding 1-ks boxcars in the union band of the
    look-elsewhere-corrected Poisson tail. rate_U in counts per (m^2 s)."""
    if iv is None:
        iv = intervals(ra, dec, mjd0, mjd1)
    if iv is None or iv["live"].sum() < FREEZE["boxcar_min_live_s"]:
        return 0.0, None
    box = FREEZE["boxcar_s"] / 86400.0; step = FREEZE["boxcar_step_s"] / 86400.0
    n_ph, tph = count(ph, mjd0, mjd1, "U")
    st = met_to_mjd(iv["start"]); sp = met_to_mjd(iv["stop"])
    # cumulative exposure/livetime vs time for fast boxcar sums
    cum_e = np.r_[0.0, np.cumsum(iv["expo_U"])]; cum_l = np.r_[0.0, np.cumsum(iv["live"])]
    n_box_total = max(iv["live"].sum() / FREEZE["boxcar_s"], 1.0)
    best = 0.0; best_rec = None
    t = mjd0
    while t + box <= mjd1 + 1e-9:
        i0 = np.searchsorted(sp, t); i1 = np.searchsorted(st, t + box)
        live = cum_l[i1] - cum_l[i0]
        if live >= FREEZE["boxcar_min_live_s"]:
            expo = cum_e[i1] - cum_e[i0]
            n = int(((tph >= t) & (tph < t + box)).sum())
            lam = rate_U * expo
            if n > 0:
                p = poisson.sf(n - 1, lam)
                pc = 1.0 - (1.0 - p) ** n_box_total
                s = float(-np.log10(max(pc, 1e-300)))
                if s > best:
                    best = s; best_rec = {"t0": Time(t, format="mjd").isot, "n": n, "lam": float(lam), "live_s": float(live)}
        t += step
    return best, best_rec


def window_counts(ph, ra, dec, mjd0, mjd1, strict=False):
    iv = intervals(ra, dec, mjd0, mjd1, strict=strict)
    rec = {"live_s": 0.0, "moon_excluded_s": 0.0}
    for lane in FREEZE["lanes"]:
        rec["expo_" + lane] = 0.0
        rec["n_" + lane] = 0
    if iv is None:
        rec["iv"] = None
        return rec
    rec["live_s"] = float(iv["live"].sum()); rec["moon_excluded_s"] = iv["moon_excluded_s"]
    rec["live_raw_s"] = float(iv["live_raw"].sum())
    for lane in FREEZE["lanes"]:
        rec["expo_" + lane] = float(iv["expo_" + lane].sum())
        rec["n_" + lane] = count(ph, mjd0, mjd1, lane)[0]
    rec["iv"] = iv
    return rec


def unit_analysis(ch, tid, rung, ph, ra, dec, windows, strict=False, exclude_event_ids=(), inject=None):
    """Full per-unit computation: per-window counts/exposure, pseudo-window
    rates, the five statistics, the dispersion gate and the pseudo-stack
    ensemble. `inject` (optional) is a callable(window, lane) -> extra
    counts to add to the real counts (completeness)."""
    rows = []
    for w in windows:
        if w["event_id"] in exclude_event_ids:
            continue
        rec = window_counts(ph, ra, dec, w["mjd0"], w["mjd1"], strict=strict)
        searchable = rec["live_s"] >= FREEZE["window_livetime_min_s"]
        pw = []
        if searchable:
            for j in pseudo_offsets(rung, w["dur_d"]):
                off = j * max(w["dur_d"], 1.0)
                r = window_counts(ph, ra, dec, w["mjd0"] + off, w["mjd1"] + off, strict=strict)
                r.pop("iv", None)
                r["offset_d"] = off
                r["kept"] = r["expo_L"] >= FREEZE["pw_min_expo_frac"] * rec["expo_L"] and r["live_s"] > 0
                pw.append(r)
        row = {"event_id": w["event_id"], "t_ca_utc": w["t_ca_utc"], "b_rsun": w["b_rsun"], "mjd0": w["mjd0"], "mjd1": w["mjd1"],
               "dur_d": w["dur_d"], "searchable": searchable, "live_s": rec["live_s"], "live_raw_s": rec.get("live_raw_s", 0.0),
               "moon_excluded_s": rec["moon_excluded_s"], "n_pw": int(sum(r["kept"] for r in pw)),
               "pw": pw, "iv": rec.get("iv")}
        for lane in FREEZE["lanes"]:
            row["expo_" + lane] = rec["expo_" + lane]
            row["n_" + lane] = rec["n_" + lane] + (inject(w, lane) if inject else 0)
        rows.append(row)
    stats = {}
    diag = {}
    searchable = [r for r in rows if r["searchable"] and r["n_pw"] >= 4]
    for lane in ("L", "H", "U"):
        # per-window rates from the kept pseudo-windows
        for r in rows:
            kept = [p for p in r["pw"] if p["kept"]]
            se = sum(p["expo_" + lane] for p in kept); sn = sum(p["n_" + lane] for p in kept)
            # amendment v1.1: Jeffreys floor (sn + 1/2)/se so that an empty
            # pseudo-window pool never yields lambda = 0 (lane H has ~0.01
            # counts per grazing window)
            r["rate_" + lane] = ((sn + 0.5) / se) if se > 0 else np.nan
            r["lam_" + lane] = r["rate_" + lane] * r["expo_" + lane] if se > 0 else np.nan
        if lane == "U":
            continue
        n_tot = sum(r["n_" + lane] for r in searchable); lam_tot = sum(r["lam_" + lane] for r in searchable)
        stats["S_stack_" + lane] = nlog10_sf(n_tot, lam_tot) if searchable else np.nan
        ev = [(nlog10_sf(r["n_" + lane], r["lam_" + lane]), r["t_ca_utc"]) for r in searchable]
        stats["S_event_" + lane] = max(ev)[0] if ev else np.nan
        diag["S_event_" + lane + "_window"] = max(ev)[1] if ev else None
        diag["n_tot_" + lane] = n_tot; diag["lam_tot_" + lane] = float(lam_tot)
        # dispersion gate + pseudo-stack ensemble (members = offset index)
        pvals = []
        for r in searchable:
            kept = [p for p in r["pw"] if p["kept"]]
            for p in kept:
                # leave-one-out rate
                se = sum(q["expo_" + lane] for q in kept if q is not p); sn = sum(q["n_" + lane] for q in kept if q is not p)
                if se > 0:
                    lam = (sn + 0.5) / se * p["expo_" + lane]      # amendment v1.1 floor
                    pvals.append(poisson.sf(p["n_" + lane] - 1, lam) if p["n_" + lane] > 0 else 1.0)
        pvals = np.array(pvals)
        frac = float((pvals < FREEZE["dispersion_gate"][0]).mean()) if len(pvals) else np.nan
        diag["disp_frac_" + lane] = frac; diag["n_pw_" + lane] = int(len(pvals))
        diag["pw_exceed_T_" + lane] = int((-np.log10(np.maximum(pvals, 1e-300)) > FREEZE["T"]).sum()) if len(pvals) else 0
        # pseudo-stack ensemble: for offset index j, sum over windows
        offsets = sorted(set(p["offset_d"] / max(r["dur_d"], 1.0) for r in searchable for p in r["pw"] if p["kept"]))
        ens = []
        for j in offsets:
            n = 0; lam = 0.0; ok = False
            for r in searchable:
                for p in r["pw"]:
                    if p["kept"] and abs(p["offset_d"] / max(r["dur_d"], 1.0) - j) < 1e-6:
                        n += p["n_" + lane]; lam += r["rate_" + lane] * p["expo_" + lane]; ok = True
            if ok:
                ens.append(nlog10_sf(n, lam))
        diag["stack_ens_" + lane] = ens
    # burst
    sb = 0.0; sb_rec = None
    for r in searchable:
        s, rec = burst_stat(ph, ra, dec, r["mjd0"], r["mjd1"], r["rate_U"], iv=r["iv"])
        if inject is not None and hasattr(inject, "burst"):
            s2, rec2 = inject.burst(r, ph, ra, dec)
            if s2 > s:
                s, rec = s2, rec2
        if s > sb:
            sb, sb_rec = s, dict(rec or {}, window=r["t_ca_utc"])
    stats["S_burst"] = sb if searchable else np.nan
    diag["S_burst_rec"] = sb_rec
    for r in rows:
        r.pop("iv", None)
    return {"channel": ch, "target": tid, "rung": rung, "n_windows": len(rows), "n_searchable": len(searchable),
            "stats": stats, "diag": diag, "windows": rows}


def all_units(events=None):
    events = events or R.load_events()
    pos = R.positions(events)
    units = []
    for (ch, tid), (ra, de, n) in sorted(pos.items()):
        for rung in FREEZE["rungs"]:
            ws = unit_windows(ch, tid, rung, events)
            if not ws:
                continue
            units.append({"channel": ch, "target": tid, "rung": rung, "ra": ra, "dec": de, "key": f"{ch}-{tid}",
                          "dev": tid in FREEZE["dev_targets"] or (ch, tid, rung) in FREEZE["forced_dev_units"],
                          "n_windows": len(ws)})
    return units


def power_w(flux_ph_cm2_s, lane, rung):
    """P = F_E * pi b^2 (flat-top footprint of radius b_rung), watts."""
    fe = flux_ph_cm2_s * lane_mean_energy_mev(lane) * MEV_ERG          # erg cm^-2 s^-1
    b_cm = RUNG_R_AU[rung] * CM_PER_AU
    return fe * np.pi * b_cm ** 2 * 1e-7


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if (isinstance(o, float) or isinstance(o, np.floating)) and np.isnan(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    return o
