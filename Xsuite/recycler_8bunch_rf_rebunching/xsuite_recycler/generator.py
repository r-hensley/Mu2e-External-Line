"""Compact statistical generator for the historical longitudinal input.

The model intentionally contains only the structure supported by the checked
historical distribution: two trains on a common microbunch grid and common,
Gaussian within-microbunch time and momentum residuals.  It does not add
microbunch-specific jitter or non-Gaussian tails.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math

import numpy as np

from .provenance import array_sha256


MODEL_NAME = "independent-gaussian-microbunch-grid"
MODEL_VERSION = 1
N_MICROBUNCHES = 168
MICROBUNCHES_PER_BUNCH = 21
N_BUNCHES = N_MICROBUNCHES // MICROBUNCHES_PER_BUNCH
# The fitted input occupies two trains of 84 RF-grid slots. Grid indices
# 84--125 form the 42-slot gap between them.
OCCUPIED_GRID_INDICES = np.concatenate(
    (
        np.arange(0, 84, dtype=np.int16),
        np.arange(126, 210, dtype=np.int16),
    )
)


def _require_finite(name: str, value: float) -> None:
    """Raise a field-specific error when a floating-point value is not finite."""
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")


@dataclass(frozen=True)
class MicrobunchGeneratorConfig:
    """Parameters of the fitted 168-microbunch Gaussian input model.

    ``time_origin_ns`` locates the first occupied grid slot and ``spacing_ns``
    separates successive RF-grid slots.  ``sigma_time_ns`` and ``sigma_dpop``
    set the within-microbunch Gaussian widths; ``correlation`` is the Pearson
    correlation of those two residuals rather than a shift of the grid centers.
    """

    particles_per_microbunch: int = 6_400
    time_origin_ns: float = -593.7974058027177
    spacing_ns: float = 19.154736447320314
    sigma_time_ns: float = 5.0249773002015825
    mean_dpop: float = 8.20465626423e-8
    sigma_dpop: float = 1.759398503806553e-4
    correlation: float = 0.0
    seed: int = 202_208

    def __post_init__(self) -> None:
        """Validate counts, RNG seed, finite values, widths, and correlation."""
        if isinstance(self.particles_per_microbunch, (bool, np.bool_)) or not isinstance(
            self.particles_per_microbunch, (int, np.integer)
        ):
            raise TypeError("particles_per_microbunch must be an integer")
        if self.particles_per_microbunch <= 0:
            raise ValueError("particles_per_microbunch must be positive")

        if isinstance(self.seed, (bool, np.bool_)) or not isinstance(
            self.seed, (int, np.integer)
        ):
            raise TypeError("seed must be an integer")
        if self.seed < 0 or self.seed >= 2**64:
            raise ValueError("seed must be in the range [0, 2**64)")

        for name in (
            "time_origin_ns",
            "spacing_ns",
            "sigma_time_ns",
            "mean_dpop",
            "sigma_dpop",
            "correlation",
        ):
            _require_finite(name, float(getattr(self, name)))

        if self.spacing_ns <= 0:
            raise ValueError("spacing_ns must be positive")
        if self.sigma_time_ns <= 0:
            raise ValueError("sigma_time_ns must be positive")
        if self.sigma_dpop <= 0:
            raise ValueError("sigma_dpop must be positive")
        if not -1.0 <= self.correlation <= 1.0:
            raise ValueError("correlation must lie in [-1, 1]")

    @property
    def n_particles(self) -> int:
        """Return the total population across all 168 microbunches."""
        return N_MICROBUNCHES * int(self.particles_per_microbunch)

    def canonical_metadata(self) -> dict[str, object]:
        """Return the complete, JSON-serializable model specification."""
        return {
            "model": MODEL_NAME,
            "model_version": MODEL_VERSION,
            "particles_per_microbunch": int(self.particles_per_microbunch),
            "time_origin_ns": float(self.time_origin_ns),
            "spacing_ns": float(self.spacing_ns),
            "sigma_time_ns": float(self.sigma_time_ns),
            "mean_dpop": float(self.mean_dpop),
            "sigma_dpop": float(self.sigma_dpop),
            "correlation": float(self.correlation),
            "seed": int(self.seed),
            "occupied_grid_indices": [int(value) for value in OCCUPIED_GRID_INDICES],
            "microbunches_per_bunch": MICROBUNCHES_PER_BUNCH,
        }

    @property
    def fingerprint(self) -> str:
        """SHA-256 of the canonical model configuration."""
        payload = json.dumps(
            self.canonical_metadata(),
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        return sha256(payload).hexdigest()


@dataclass(frozen=True)
class GeneratedLongitudinalDistribution:
    """Generated coordinates and exact labels needed for validation/tracking.

    ``dt_s`` and ``dpop`` contain one value per particle.  ``center_time_s``
    repeats the corresponding ideal center for every particle, whereas
    ``microbunch_center_times_s`` contains the 168 unique centers.  Microbunch
    and bunch indices retain the structural origin of each random sample. The
    configuration fingerprint identifies model inputs; the realized-coordinate
    digest additionally identifies the exact sampled ``dt_s`` and ``dpop``
    arrays produced by the recorded NumPy RNG implementation.
    """

    dt_s: np.ndarray
    dpop: np.ndarray
    microbunch_index: np.ndarray
    bunch_index: np.ndarray
    center_time_s: np.ndarray
    microbunch_center_times_s: np.ndarray
    occupied_grid_indices: np.ndarray
    config_fingerprint: str
    realized_coordinate_sha256: str
    metadata: dict[str, object]

    @property
    def n_particles(self) -> int:
        """Return the number of generated macroparticles."""
        return int(self.dt_s.size)


def generate_longitudinal_distribution(
    config: MicrobunchGeneratorConfig | None = None,
) -> GeneratedLongitudinalDistribution:
    """Sample the compact 168-microbunch longitudinal model.

    Particles are returned in contiguous microbunch blocks.  The eight bunch
    labels are structural labels only: every consecutive set of 21 occupied
    microbunches belongs to one bunch.
    """
    if config is None:
        config = MicrobunchGeneratorConfig()
    if not isinstance(config, MicrobunchGeneratorConfig):
        raise TypeError("config must be a MicrobunchGeneratorConfig")

    particles_per_microbunch = int(config.particles_per_microbunch)
    microbunch_index = np.repeat(
        np.arange(N_MICROBUNCHES, dtype=np.int16), particles_per_microbunch
    )
    # Every consecutive 21 occupied microbunches is one structural bunch, so
    # the 168 occupied microbunches produce eight equal bunch labels.
    bunch_index = (microbunch_index // MICROBUNCHES_PER_BUNCH).astype(
        np.int8, copy=False
    )

    microbunch_center_times_s = (
        float(config.time_origin_ns)
        + OCCUPIED_GRID_INDICES.astype(np.float64) * float(config.spacing_ns)
    ) * 1e-9
    center_time_s = np.repeat(
        microbunch_center_times_s, particles_per_microbunch
    )

    rng = np.random.default_rng(int(config.seed))
    standardized_time = rng.standard_normal(config.n_particles)
    independent_momentum = rng.standard_normal(config.n_particles)

    dt_s = center_time_s + float(config.sigma_time_ns) * 1e-9 * standardized_time
    rho = float(config.correlation)
    # Standard correlated-normal construction: the shared time deviate creates
    # the requested correlation, while max guards |rho|=1 against roundoff.
    standardized_momentum = (
        rho * standardized_time
        + math.sqrt(max(0.0, 1.0 - rho * rho)) * independent_momentum
    )
    dpop = float(config.mean_dpop) + float(config.sigma_dpop) * standardized_momentum

    # Hash the realized random coordinates separately from the model
    # fingerprint. The same configuration can otherwise produce different
    # arrays if a future NumPy release changes RNG implementation details.
    dt_s = np.ascontiguousarray(dt_s)
    dpop = np.ascontiguousarray(dpop)
    realized_coordinate_sha256 = array_sha256(dt_s, dpop)

    metadata = config.canonical_metadata()
    metadata.update(
        {
            "config_fingerprint": config.fingerprint,
            "n_microbunches": N_MICROBUNCHES,
            "n_bunches": N_BUNCHES,
            "n_particles": config.n_particles,
            "numpy_version": np.__version__,
            "rng_bit_generator": type(rng.bit_generator).__name__,
            "realized_coordinate_sha256": realized_coordinate_sha256,
        }
    )

    return GeneratedLongitudinalDistribution(
        dt_s=dt_s,
        dpop=dpop,
        microbunch_index=np.ascontiguousarray(microbunch_index),
        bunch_index=np.ascontiguousarray(bunch_index),
        center_time_s=np.ascontiguousarray(center_time_s),
        microbunch_center_times_s=np.ascontiguousarray(microbunch_center_times_s),
        occupied_grid_indices=OCCUPIED_GRID_INDICES.copy(),
        config_fingerprint=config.fingerprint,
        realized_coordinate_sha256=realized_coordinate_sha256,
        metadata=metadata,
    )
