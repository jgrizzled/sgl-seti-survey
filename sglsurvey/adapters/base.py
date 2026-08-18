"""Adapter protocols: the boundary between the survey pipeline and any
specific archive (plan §6). The pipeline and sglseti see footprints only
through ``FootprintTest.contains`` — the exact signature that
``sglseti.covered_z_intervals(contains=...)`` consumes."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Protocol, Sequence, runtime_checkable

from sglsurvey.records import Observation
from sglsurvey.snapshots import SnapshotStore


@dataclass(frozen=True)
class MjdRange:
    start_mjd_utc: float
    stop_mjd_utc: float


@dataclass(frozen=True)
class ConeRegion:
    """Pilot-scale discovery region; MOC/polygon forms come later."""

    ra_deg: float
    dec_deg: float
    radius_deg: float


DiscoveryRegion = ConeRegion


@dataclass(frozen=True)
class CutoutSpec:
    ra_deg: float
    dec_deg: float
    size_pix: int


@dataclass(frozen=True)
class LocalProduct:
    kind: str  # e.g. "int", "msk", "unc"
    path: Path
    checksum: str  # verified against the archive-published checksum


@dataclass(frozen=True)
class ProductSet:
    observation_id: str
    products: tuple[LocalProduct, ...]
    cutout: CutoutSpec | None = None


@runtime_checkable
class FootprintTest(Protocol):
    def contains(self, ra_deg: float, dec_deg: float) -> bool: ...


class ArchiveAdapter(Protocol):
    """One implementation per archive family (e.g. IRSA WISE L1b)."""

    archive_id: str

    def discover(self, region: DiscoveryRegion, time_range: MjdRange,
                 store: SnapshotStore) -> Iterator[Observation]:
        """Query archive metadata, snapshot raw responses verbatim, and
        yield normalized Observation records."""
        ...

    def nominal_footprint(self, obs: Observation) -> FootprintTest:
        """Cheap containment test from inventory metadata (corner
        polygon / nominal WCS). Coarse pass only — not usability."""
        ...

    def exact_footprint(self, obs: Observation,
                        products: ProductSet) -> FootprintTest:
        """Exact spherical WCS + mask + quality cuts over downloaded
        products; defines usable pixels (plan §3.3)."""
        ...

    def fetch(self, obs: Observation, kinds: Sequence[str], dest: Path,
              *, cutout: CutoutSpec | None = None) -> ProductSet:
        """Download products (full frames or cutouts), verify published
        checksums, and return local paths with hashes."""
        ...
