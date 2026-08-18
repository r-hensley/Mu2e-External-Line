"""Load the historical BLonD distribution and map it into Xsuite coordinates."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .config import C_LIGHT_M_PER_S, RecyclerConfig
from .generator import GeneratedLongitudinalDistribution, MODEL_NAME, N_MICROBUNCHES
from .provenance import array_sha256


@dataclass(frozen=True)
class InitialDistribution:
    """Longitudinal particles and provenance prepared for Xsuite tracking.

    ``zeta_m`` and ``ptau`` are Xsuite longitudinal coordinates, while
    ``bunch_index`` preserves the eight structural bunch labels used for plots.
    ``macro_equivalent_weight`` scales a reduced sample back to the historical
    macroparticle count; ``physical_proton_weight`` instead scales it to the
    configured physical beam intensity. The SHA-256 fields distinguish the
    complete source/configuration, realized source coordinates, final tracking
    coordinates, and any selected row indices. Sampling fields state how that
    realized ensemble was chosen from the source.
    """

    zeta_m: np.ndarray
    ptau: np.ndarray
    bunch_index: np.ndarray
    source_rows: int
    macro_equivalent_weight: float
    physical_proton_weight: float
    source_path: str
    source_sha256: str
    source_coordinate_sha256: str
    tracking_coordinate_sha256: str
    selected_indices_sha256: str
    sampling_method: str
    sampling_seed: int | None
    sampling_strata: int
    source_kind: str = "file"
    source_metadata_json: str = "{}"

    @property
    def n_particles(self) -> int:
        """Return the number of tracked macroparticles in this sample."""
        return int(self.zeta_m.size)


def dt_to_zeta(dt_s, beta0: float):
    """Convert BLonD arrival-time offset in seconds to Xsuite ``zeta`` in metres.

    This project uses ``zeta=-beta0*c*dt``; a later-than-reference arrival
    therefore has negative ``zeta``.
    """
    return -beta0 * C_LIGHT_M_PER_S * np.asarray(dt_s)


def zeta_to_dt(zeta_m, beta0: float):
    """Convert Xsuite ``zeta`` in metres to BLonD arrival-time offset in seconds."""
    return -np.asarray(zeta_m) / (beta0 * C_LIGHT_M_PER_S)


def _stratified_indices(
    n_rows: int,
    requested: int,
    n_strata: int,
    seed: int,
) -> np.ndarray:
    """Select rows without replacement from equal contiguous file strata.

    The historical file is ordered by bunch and longitudinal structure.  Taking
    a controlled number from every contiguous stratum retains that structure
    more reliably than an unstratified random sample.
    """
    if n_strata <= 0:
        raise ValueError("n_strata must be positive")
    if n_strata > n_rows:
        raise ValueError("n_strata cannot exceed the number of input rows")
    if requested < n_strata:
        raise ValueError("Sample at least one particle from each sampling stratum")

    # array_split allows strata to differ by at most one row when the file size
    # is not exactly divisible by the requested stratum count.
    blocks = np.array_split(np.arange(n_rows, dtype=np.int64), n_strata)
    per_block = [requested // n_strata] * n_strata
    for ii in range(requested % n_strata):
        per_block[ii] += 1

    rng = np.random.default_rng(seed)
    chosen = []
    for block, count in zip(blocks, per_block):
        if count > block.size:
            raise ValueError("Requested sample is larger than an input stratum")
        take = np.sort(rng.choice(block, size=count, replace=False))
        chosen.append(take)
    return np.concatenate(chosen)


def _bunch_labels_for_rows(
    selected: np.ndarray,
    n_rows: int,
    n_bunches: int,
) -> np.ndarray:
    """Label selected ordered rows using the file's equal contiguous bunch blocks."""
    base, remainder = divmod(n_rows, n_bunches)
    sizes = np.full(n_bunches, base, dtype=np.int64)
    sizes[:remainder] += 1
    boundaries = np.cumsum(sizes)
    return np.searchsorted(boundaries, selected, side="right").astype(np.int16)


def load_initial_distribution(
    path: str | Path,
    machine: RecyclerConfig,
    max_particles: int | None = None,
    seed: int = 202208,
    sampling_strata: int | None = None,
) -> InitialDistribution:
    """Load, validate, optionally sample, and convert a historical input file.

    The text file must contain exactly two columns: the historical longitudinal
    coordinate in metres and fractional momentum offset ``dpop``.  Its complete
    contents are hashed for provenance even when only a stratified subset is
    tracked.  Returned coordinates use the conversions from the committed
    BLonD driver, and returned weights describe both the original macroensemble
    and the physical beam intensity.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Initial distribution not found: {path}")

    # Hash the full source before sampling so the result identifies the exact
    # historical input, not merely the selected subset.
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)

    raw = np.loadtxt(path, dtype=np.float64)
    if raw.ndim != 2 or raw.shape[1] != 2:
        raise ValueError(
            f"Expected exactly two columns (z [m], dpop); got shape {raw.shape}"
        )
    if not np.isfinite(raw).all():
        raise ValueError("Initial distribution contains non-finite values")

    n_rows = int(raw.shape[0])
    if n_rows < machine.n_bunches:
        raise ValueError("Initial distribution has fewer rows than bunches")

    # None, zero, and negative limits all mean "track every input row".
    if max_particles is None or max_particles <= 0 or max_particles >= n_rows:
        selected = np.arange(n_rows, dtype=np.int64)
        sampling_method = "all-rows"
        effective_sampling_seed = None
        effective_sampling_strata = 0
    else:
        n_strata = machine.n_bunches if sampling_strata is None else int(sampling_strata)
        # A whole number of strata per bunch avoids giving any of the eight
        # contiguous bunch blocks systematically different sampling coverage.
        if n_strata % machine.n_bunches != 0:
            raise ValueError("sampling_strata must be divisible by the eight bunches")
        selected = _stratified_indices(
            n_rows, max_particles, n_strata, seed
        )
        sampling_method = "contiguous-stratified-random-without-replacement"
        effective_sampling_seed = int(seed)
        effective_sampling_strata = n_strata
    bunch_index = _bunch_labels_for_rows(selected, n_rows, machine.n_bunches)

    z_input_m = raw[selected, 0]
    dpop_input = raw[selected, 1]
    del raw
    source_coordinate_sha256 = array_sha256(z_input_m, dpop_input)
    selected_indices_sha256 = array_sha256(selected)

    # Preserve the source driver's approximate 3e8 m/s conversion here before
    # mapping arrival time into Xsuite's zeta convention with the true beta0.
    dt_s = z_input_m / machine.input_time_denominator_m_per_s
    zeta_m = dt_to_zeta(dt_s, machine.beta0)

    # The committed BLonD code computes dE=dpop*p0c/0.9944. Xsuite defines
    # ptau=(E-E0)/p0c, hence the exact code-faithful mapping below.
    ptau = dpop_input / machine.blond_energy_beta
    n_selected = int(selected.size)
    zeta_m = np.ascontiguousarray(zeta_m)
    ptau = np.ascontiguousarray(ptau)
    bunch_index = np.ascontiguousarray(bunch_index)
    tracking_coordinate_sha256 = array_sha256(zeta_m, ptau, bunch_index)
    source_metadata = {
        "requested_max_particles": max_particles,
        "sampling_method": sampling_method,
        "sampling_seed": effective_sampling_seed,
        "sampling_strata": effective_sampling_strata,
        "selected_indices_sha256": selected_indices_sha256,
        "selected_source_coordinate_sha256": source_coordinate_sha256,
        "tracking_coordinate_sha256": tracking_coordinate_sha256,
    }

    return InitialDistribution(
        zeta_m=zeta_m,
        ptau=ptau,
        bunch_index=bunch_index,
        source_rows=n_rows,
        macro_equivalent_weight=n_rows / n_selected,
        physical_proton_weight=machine.intensity_protons / n_selected,
        source_path=str(path.resolve()),
        source_sha256=digest.hexdigest(),
        source_coordinate_sha256=source_coordinate_sha256,
        tracking_coordinate_sha256=tracking_coordinate_sha256,
        selected_indices_sha256=selected_indices_sha256,
        sampling_method=sampling_method,
        sampling_seed=effective_sampling_seed,
        sampling_strata=effective_sampling_strata,
        source_metadata_json=json.dumps(
            source_metadata,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ),
    )


def initial_distribution_from_generated(
    generated: GeneratedLongitudinalDistribution,
    machine: RecyclerConfig,
    *,
    reference_particles_per_microbunch: int = 6_400,
) -> InitialDistribution:
    """Map a generated sample into Xsuite coordinates with explicit scaling.

    A reduced generated ensemble represents the same 1,075,200-macroparticle
    reference population as the historical file, so its histogram weight is
    6,400 divided by its generated particles per microbunch.
    """
    if not isinstance(generated, GeneratedLongitudinalDistribution):
        raise TypeError("generated must be a GeneratedLongitudinalDistribution")
    if reference_particles_per_microbunch <= 0:
        raise ValueError("reference_particles_per_microbunch must be positive")
    if generated.n_particles <= 0 or generated.n_particles % N_MICROBUNCHES:
        raise ValueError("Generated sample must have equal nonzero microbunch populations")
    if not np.isfinite(generated.dt_s).all() or not np.isfinite(generated.dpop).all():
        raise ValueError("Generated distribution contains non-finite values")

    reference_rows = N_MICROBUNCHES * int(reference_particles_per_microbunch)
    zeta_m = np.ascontiguousarray(dt_to_zeta(generated.dt_s, machine.beta0))
    ptau = np.ascontiguousarray(generated.dpop / machine.blond_energy_beta)
    bunch_index = np.ascontiguousarray(generated.bunch_index, dtype=np.int16)
    tracking_coordinate_sha256 = array_sha256(zeta_m, ptau, bunch_index)
    return InitialDistribution(
        zeta_m=zeta_m,
        ptau=ptau,
        bunch_index=bunch_index,
        source_rows=reference_rows,
        macro_equivalent_weight=reference_rows / generated.n_particles,
        physical_proton_weight=machine.intensity_protons / generated.n_particles,
        source_path=f"generated://{MODEL_NAME}",
        source_sha256=generated.config_fingerprint,
        source_coordinate_sha256=generated.realized_coordinate_sha256,
        tracking_coordinate_sha256=tracking_coordinate_sha256,
        selected_indices_sha256="",
        sampling_method="generated-model",
        sampling_seed=int(generated.metadata["seed"]),
        sampling_strata=0,
        source_kind="generated",
        source_metadata_json=json.dumps(
            {
                **generated.metadata,
                "tracking_coordinate_sha256": tracking_coordinate_sha256,
            },
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ),
    )
