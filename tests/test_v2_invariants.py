"""Invariant tests for the v2 shared code (wise v2_plan §4.8, §8.3).

Fast, data-free tests run everywhere; tests that need the v1 products
under runs/ are skipped when they are absent.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from sglsurvey import inject, manifest, nulls
from sglsurvey.photometry import (_gaussian_kernel, flux_map_from_arrays,
                                  matched_filter)

REPO = Path(__file__).resolve().parents[1]
PRF_DIR = REPO / "runs" / "wise" / "v2" / "prf"
V1_TENSORS = REPO / "runs" / "wise" / "calib_v1" / "tensors"


# -- kernels / PRF --------------------------------------------------------
def test_gaussian_kernel_unit_sum():
    k = _gaussian_kernel(6.1 / 2.75, 6)
    assert abs(k.sum() - 1.0) < 1e-9
    assert k.shape == (13, 13) and k[6, 6] == k.max()


@pytest.mark.skipif(not (PRF_DIR / "w1").exists(), reason="PRF files not fetched")
def test_prf_render_normalisation_and_shift():
    prf = inject.PRFGrid.load("W1", PRF_DIR / "w1", elements={(4, 4)})
    s0, ox, oy = prf.render(50.0, 50.0, element=(4, 4))
    assert 0.97 < s0.sum() <= 1.0001            # enclosed flux in 33 px
    assert np.unravel_index(s0.argmax(), s0.shape) == (50 - oy, 50 - ox)
    s1, ox1, oy1 = prf.render(50.5, 50.0, element=(4, 4))
    # sub-pixel placement: the stamp centroid moves by exactly the shift
    # (the empirical template is itself slightly asymmetric, so compare
    # centroids rather than pixel symmetry)
    def cx(s, ox):
        xx = ox + np.arange(s.shape[1])
        return (s.sum(axis=0) * xx).sum() / s.sum()
    assert abs((cx(s1, ox1) - cx(s0, ox)) - 0.5) < 0.02
    assert abs(s1.sum() - s0.sum()) < 0.01


def test_element_of_corners():
    g = inject.PRFGrid(band="W1", templates={})
    assert g.element_of(0, 0) == (0, 0)
    assert g.element_of(1015, 1015) == (8, 8)
    assert g.element_of(508, 508) == (4, 4)
    g4 = inject.PRFGrid(band="W4", templates={})
    assert g4.element_of(507, 0) == (8, 0)


# -- spectrum / magnitudes ------------------------------------------------
def test_dn_mag_round_trip():
    for mag in (10.0, 14.3, 17.9):
        dn = inject.dn_from_vega_mag(mag, 20.752)
        assert abs(20.752 - 2.5 * np.log10(dn) - mag) < 1e-9


def test_fnu_mag_round_trip_and_colour_corrections():
    for band in ("W1", "W2", "W3", "W4"):
        for spec in ("flat_fnu", "blackbody_300K", "g2v"):
            f = inject.fnu_from_vega_mag(band, 13.0, spec)
            assert abs(inject.vega_mag_from_fnu(band, f, spec) - 13.0) < 1e-9
    bb = inject.COLOR_CORRECTIONS["blackbody_300K"]
    a, b = inject.COLOR_CORRECTIONS["blackbody_283K"], inject.COLOR_CORRECTIONS["blackbody_400K"]
    for k in bb:
        assert min(a[k], b[k]) - 1e-9 <= bb[k] <= max(a[k], b[k]) + 1e-9
    assert inject.COLOR_CORRECTIONS["powerlaw_nu-2"]["W3"] == 1.0


# -- injection linearity --------------------------------------------------
def _synthetic_map(seed=0, ny=120, nx=110):
    rng = np.random.default_rng(seed)
    img = 30.0 + rng.normal(0, 3.0, (ny, nx))
    img[40:43, 60:63] += 400.0                     # a star
    var = np.full(img.shape, 9.0)
    good = np.ones(img.shape, bool)
    good[70:80, 20:30] = False                     # a masked block
    good[55, 56] = False
    return img, var, good


def test_zero_flux_injection_changes_nothing():
    img, var, good = _synthetic_map()
    fm0 = flux_map_from_arrays(img, var, good, None, "W1", 57000.0, keep_inputs=True)
    fm1 = flux_map_from_arrays(img + 0.0, var, good, None, "W1", 57000.0)
    assert np.array_equal(np.nan_to_num(fm0.flux), np.nan_to_num(fm1.flux))


@pytest.mark.skipif(not (PRF_DIR / "w1").exists(), reason="PRF files not fetched")
def test_stamp_response_equals_full_reconvolution():
    prf = inject.PRFGrid.load("W1", PRF_DIR / "w1", elements={(4, 4)})
    img, var, good = _synthetic_map()
    fm0 = flux_map_from_arrays(img, var, good, None, "W1", 57000.0, keep_inputs=True)
    rng = np.random.default_rng(3)
    for _ in range(5):
        x, y = rng.uniform(10, 100), rng.uniform(10, 110)
        amp = 300.0
        stamp, ox, oy = prf.render(x, y, element=(4, 4))
        fm1 = flux_map_from_arrays(inject.add_stamp(img, stamp, ox, oy, amp),
                                   var, good, None, "W1", 57000.0)
        rw = inject.stamp_response(fm0, stamp, ox, oy)
        full = np.nan_to_num(fm1.flux - fm0.flux)
        usable = fm0.good_frac >= nulls.MIN_GOOD_FRAC   # as sampled by the stack
        sl = (slice(rw.oy, rw.oy + rw.delta.shape[0]), slice(rw.ox, rw.ox + rw.delta.shape[1]))
        win_err = np.abs(full[sl] - amp * rw.delta)[usable[sl]]
        assert np.max(win_err) < 2e-3 * amp
        outside = usable.copy()
        outside[sl] = False
        assert np.max(np.abs(full[outside])) < 2e-3 * amp
        # response at the source position is the Gaussian-vs-PRF throughput
        r0 = rw.sample(x, y)
        assert 0.6 < r0 < 1.0
        # the injected source landing on masked pixels is attenuated
    stamp, ox, oy = prf.render(25.0, 75.0, element=(4, 4))   # inside the mask block
    rw = inject.stamp_response(fm0, stamp, ox, oy)
    assert rw.sample(25.0, 75.0) < 0.3


# -- interval tiling ------------------------------------------------------
def reciprocal_midpoint_edges(z_nodes: np.ndarray, n_intervals: int,
                              z_min: float = 550.0, z_max: float = 10000.0):
    """Interval edges at reciprocal-distance midpoints tiling [z_min, z_max]."""
    q = np.sort(1.0 / np.asarray(z_nodes))
    blocks = np.array_split(np.arange(len(q)), n_intervals)
    edges = [1.0 / q[0]]
    for b in blocks[:-1]:
        edges.append(2.0 / (q[b[-1]] + q[b[-1] + 1]))
    edges.append(1.0 / q[-1])
    edges = np.array(edges)
    edges[0], edges[-1] = max(edges[0], z_max), min(edges[-1], z_min)
    return edges[::-1] if edges[0] > edges[-1] else edges


def test_intervals_tile_without_gaps():
    z = 1.0 / np.linspace(1 / 550.0, 1 / 10000.0, 64)
    e = reciprocal_midpoint_edges(z, 8)
    assert e[0] == 550.0 and e[-1] == 10000.0
    assert np.all(np.diff(e) > 0)
    # log-uniform prior mass covered == 1 exactly
    mass = sum(np.log(e[i + 1] / e[i]) for i in range(len(e) - 1)) / np.log(10000 / 550)
    assert abs(mass - 1.0) < 1e-12


# -- null statistics ------------------------------------------------------
def test_exceedance_ratios_loo():
    mx = np.array([10.0, 8.0, 6.0, 5.0])
    R, T = nulls.exceedance_ratios(mx, np.array([1, 2, 3]))
    assert T == 8.0 and R[0] == 10.0 / 8.0
    assert R[1] == 8.0 / 6.0 and R[2] == 6.0 / 8.0 and R[3] == 5.0 / 8.0


def test_phase_scramble_single_phase_invariant():
    rng = np.random.default_rng(0)
    E, nz, nm = 40, 16, 5
    f = rng.normal(0, 1, (E, nz, nm, nm)).astype(np.float32)
    f[:, 3, 2, 2] += 5.0
    v = np.ones_like(f); g = np.ones_like(f)
    phase = np.zeros(E, int)
    parts = nulls.phase_parts(f, v, g, phase)
    S = nulls.stack_S(f, v, g)
    sm = nulls.phase_scrambled_maxima(parts, 10, rng)
    assert np.allclose(sm, np.nanmax(S), atol=1e-4)
    phase[E // 2:] = 1
    parts = nulls.phase_parts(f, v, g, phase)
    sm = nulls.phase_scrambled_maxima(parts, 50, rng)
    assert np.median(sm) < np.nanmax(S)   # coherence broken -> lower maxima


def test_fwer_and_rank_and_bh():
    rng = np.random.default_rng(1)
    cells = [nulls.CellNull(key=f"c{i}", R_real=1.0 + 0.1 * i,
                            pooled=rng.normal(1.0, 0.1, 300)) for i in range(10)]
    res = nulls.fwer_threshold(cells, alpha=0.05, n_draws=2000, seed=0)
    assert res["n_cells"] == 10 and res["R_fwer"] > 1.0
    assert res["global_p"]["c9"] < res["global_p"]["c0"]
    rs = nulls.rank_statement(1.3, cells[0].pooled)
    assert rs["n_controls"] == 300 and 0 <= rs["p_cell"] <= 1
    lo, hi = nulls.wilson_interval(5, 100, 1.96)
    assert lo < 0.05 < hi
    q = nulls.bh_qvalues({"a": 0.001, "b": 0.04, "c": 0.5})
    assert q["a"] <= q["b"] <= q["c"] and q["a"] == 0.003


# -- manifests ------------------------------------------------------------
def test_manifest_and_stale_detection(tmp_path):
    f1, f2 = tmp_path / "a.bin", tmp_path / "b.bin"
    f1.write_bytes(b"abc"); f2.write_bytes(b"def")
    cache = manifest.HashCache(tmp_path / "cache.json")
    man = manifest.build_manifest({"a": f1, "b": f2}, stage="t", schema_version="v1",
                                  cache=cache)
    assert man["n_files"] == 2 and man["manifest_hash"].startswith("sha256:")
    digests = [e["sha256"] for e in man["files"].values()]
    manifest.check_product_inputs(manifest.combined_hash(digests), digests, "p")
    f2.write_bytes(b"xyz")
    cache2 = manifest.HashCache(tmp_path / "cache.json")
    assert cache2.sha256(f2) != man["files"]["b"]["sha256"]
    with pytest.raises(manifest.StaleProductError):
        manifest.check_product_inputs(manifest.combined_hash(digests),
                                      [cache2.sha256(f1), cache2.sha256(f2)], "p")


# -- records / corridors ----------------------------------------------------
def test_records_v2_fields_are_additive():
    import dataclasses

    from sglsurvey.records import Candidate, Constraint
    names = {f.name for f in dataclasses.fields(Constraint)}
    assert {"completeness_kind", "ci_68", "n_injections", "prf_model",
            "spectrum_model", "motion_bound_norm", "supersedes"} <= names
    names = {f.name for f in dataclasses.fields(Candidate)}
    assert {"annotations", "global_p_value", "rank_statement", "rejection_test"} <= names


def test_corridor_table():
    from sglsurvey.corridors import ALL_CORRIDORS, ALL_ENDPOINTS, CORRIDOR_OF, MEMBERS
    assert len(ALL_ENDPOINTS) == 88 and len(ALL_CORRIDORS) == 77
    assert all(e in MEMBERS[c] for e, c in CORRIDOR_OF.items())


# -- v1 products -------------------------------------------------------------
@pytest.mark.skipif(not (V1_TENSORS / "gj-625__rx.npz").exists(), reason="v1 tensors absent")
def test_stack_reproduces_v1_threshold_report():
    rep = json.loads((REPO / "runs/wise/calib_v1/threshold_report.json").read_text())
    d = np.load(V1_TENSORS / "gj-625__rx.npz")
    eb = d["band_idx"] == 1
    mx = nulls.cell_maxima(d["f"][:, eb], d["v"][:, eb], d["g"][:, eb])
    R, T = nulls.exceedance_ratios(mx, np.arange(1, 9))
    assert abs(T - rep["gj-625/rx/W1"]["T"]) < 1e-3
    assert abs(mx[0] - rep["gj-625/rx/W1"]["real_max_S"]) < 1e-3
