"""
Combine single-exposure SPHEREx measurements into a spectrum.

Each SPHEREx exposure samples one wavelength at the target position, so a
spectrum is built by stacking many exposures in wavelength bins.

SPHEREx revisits a position once per six-month survey pass, in a burst spanning
a few weeks, and one burst covers the whole wavelength range. So the epochs
group naturally into blocks, and each block is a complete spectrum. Everything
before the reference epoch is stacked into a single baseline, since the host is
constant; everything after is kept per block, because a transient evolves
between passes and averaging them together hides it.

Two error estimates are carried for every bin. 'flux_err' propagates the
VARIANCE extension and the background scatter; 'flux_scatter' is the standard
error of the individual measurements in the bin, so it measures what the
photometry actually repeats to. The formal errors underestimate the observed
repeatability, so by default they are rescaled by an inflation factor measured
from the data itself, via error_inflation_factor.
"""

import logging

import numpy as np
import pandas as pd

from obsutils.spherex.constants import BLOCK_GAP_DAYS, default_wavelength_edges

logger = logging.getLogger(__name__)

# One phase holds far fewer epochs than the whole set, so the per-bin threshold
# is relaxed when calibrating on it; the median over many bins carries the
# estimate. Below this many usable bins the phase-only estimate is unstable and
# every epoch is used instead.
CALIBRATION_MIN_PER_BIN = 3
MIN_CALIBRATION_BINS = 10

# An inflation factor above this usually means real variability is being
# absorbed into the noise estimate.
LARGE_INFLATION = 6.0


def error_inflation_factor(
    photometry: pd.DataFrame,
    edges: np.ndarray | None = None,
    min_per_bin: int = 6,
) -> float:
    """
    Measure how far the formal errors understate the epoch-to-epoch scatter.

    Within a bin the SED still has a gradient, and different epochs sample
    different wavelengths across the bin, so a straight line in log wavelength
    is removed before the residuals are compared with the formal errors. The
    returned factor is the square root of the median reduced chi-squared over
    bins, i.e. what the formal errors must be multiplied by for the scatter to
    look like noise.

    :param photometry: Table of measurements, from measure_images
    :param edges: Wavelength bin edges in micron
    :param min_per_bin: Bins with fewer measurements than this are ignored
    :return: Inflation factor, at least 1
    """
    if edges is None:
        edges = default_wavelength_edges()

    bin_index = np.digitize(photometry["wavelength"], edges) - 1

    reduced_chi2 = []
    for index in range(len(edges) - 1):
        selected = photometry[bin_index == index]
        if len(selected) < min_per_bin:
            continue

        flux = selected["flux"].to_numpy()
        error = selected["flux_err"].to_numpy()
        log_wavelength = np.log(selected["wavelength"].to_numpy())

        gradient = np.polyfit(log_wavelength, flux, 1)
        residual = flux - np.polyval(gradient, log_wavelength)

        reduced_chi2.append(np.sum((residual / error) ** 2) / (len(flux) - 2))

    if len(reduced_chi2) == 0:
        logger.warning("Too few measurements per bin to calibrate the errors")
        return 1.0

    factor = float(np.sqrt(np.median(reduced_chi2)))
    logger.info(
        f"Formal errors understate the epoch scatter by a factor {factor:.2f}, "
        f"from {len(reduced_chi2)} bins"
    )
    if factor > LARGE_INFLATION:
        logger.warning(
            f"An inflation factor of {factor:.1f} is large. If the source varies "
            f"within the epochs used, its variability is being counted as noise; "
            f"calibrate on epochs where the transient is absent instead."
        )

    return max(factor, 1.0)


def stack_in_wavelength(
    photometry: pd.DataFrame, edges: np.ndarray | None = None
) -> pd.DataFrame:
    """
    Inverse-variance weighted mean flux in each wavelength bin.

    'flux_err' is the formal error on the weighted mean; stack_by_phase
    rescales it against the epoch scatter, keeping the raw value as
    'flux_err_formal'. 'flux_scatter' is the standard error of the individual
    measurements in the bin, an independent estimate that needs no calibration
    but is itself noisy when a bin holds only a few epochs.

    :param photometry: Table of measurements, from measure_images
    :param edges: Wavelength bin edges in micron. Default spans the SPHEREx range.
    :return: One row per populated bin
    """
    if edges is None:
        edges = default_wavelength_edges()

    bin_index = np.digitize(photometry["wavelength"], edges) - 1

    rows = []
    for index in range(len(edges) - 1):
        selected = photometry[bin_index == index]
        if len(selected) == 0:
            continue

        weights = 1.0 / selected["flux_err"] ** 2
        formal = float(1.0 / np.sqrt(np.sum(weights)))
        scatter = (
            float(np.std(selected["flux"], ddof=1) / np.sqrt(len(selected)))
            if len(selected) > 1
            else np.nan
        )
        rows.append(
            {
                "wavelength": float(np.sqrt(edges[index] * edges[index + 1])),
                "wave_min": float(edges[index]),
                "wave_max": float(edges[index + 1]),
                "n": len(selected),
                "flux": float(np.sum(selected["flux"] * weights) / np.sum(weights)),
                "flux_err": formal,
                "flux_err_formal": formal,
                # stack_by_phase fills this in for bins too sparse to measure a
                # scatter of their own.
                "flux_scatter": scatter,
            }
        )

    return pd.DataFrame(rows)


def assign_blocks(photometry: pd.DataFrame) -> pd.Series:
    """
    Label each measurement with the survey pass it belongs to.

    Everything before the reference epoch becomes one 'pre' baseline. Passes
    after it are numbered in time order. Without a reference epoch there is no
    baseline to form, so every pass is labelled separately.

    :param photometry: Table of measurements, from measure_images
    :return: Block label per row, aligned to the input index
    """
    ordered = photometry.sort_values("mjd")
    epochs = ordered["mjd"].to_numpy()
    survey_pass = np.concatenate([[0], np.cumsum(np.diff(epochs) > BLOCK_GAP_DAYS)])

    labels = pd.Series(index=ordered.index, dtype=object)
    is_pre = (ordered["phase"] == "pre").to_numpy()

    if not is_pre.any():
        for number, value in enumerate(sorted(set(survey_pass)), start=1):
            labels[survey_pass == value] = f"visit {number}"
        return labels.reindex(photometry.index)

    labels[is_pre] = "pre"
    for number, value in enumerate(sorted(set(survey_pass[~is_pre])), start=1):
        labels[(survey_pass == value) & ~is_pre] = f"post {number}"

    return labels.reindex(photometry.index)


def stack_by_block(
    photometry: pd.DataFrame,
    edges: np.ndarray | None = None,
    reference_mjd: float | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Stack in wavelength per block, with calibrated errors.

    One inflation factor is measured and applied to every block, so they carry
    consistent errors. It is measured on the pre-reference epochs, where the
    transient is absent, falling back to every epoch when that phase covers too
    few wavelength bins to be stable.

    That matters for a source that really does vary: measured over all epochs,
    the transient's own variability is counted as noise and the errors inflate
    until the detection disappears. On the nova used to test this the factor is
    17 measured over all epochs against 4 measured before the eruption, which is
    the difference between a 3 sigma and a 52 sigma detection. Measuring it per
    block would fall into the same trap, so it stays global.

    :param photometry: Table of measurements, from measure_images
    :param edges: Wavelength bin edges in micron
    :param reference_mjd: Epoch the blocks are measured from, for 'phase_days'
    :return: Mapping of block label to its stack, baseline first
    """
    error_scale = _calibrate_errors(photometry, edges=edges)
    labels = assign_blocks(photometry)

    ordered_labels = sorted(
        set(labels),
        key=lambda label: (label != "pre", photometry["mjd"][labels == label].min()),
    )

    stacks = {}
    for label in ordered_labels:
        group = photometry[labels == label]
        stack = stack_in_wavelength(group, edges=edges)
        stack["flux_err"] *= error_scale
        stack["flux_scatter"] = stack["flux_scatter"].fillna(stack["flux_err"])
        stack["mjd"] = group["mjd"].mean()
        stack["phase_days"] = (
            np.nan if reference_mjd is None else group["mjd"].mean() - reference_mjd
        )
        stacks[label] = stack

    spans = ", ".join(
        f"{label} (n = {int(stack['n'].sum())})" for label, stack in stacks.items()
    )
    logger.info(f"Stacked {len(stacks)} blocks: {spans}")

    return stacks


def _calibrate_errors(photometry: pd.DataFrame, edges: np.ndarray | None) -> float:
    """
    Measure the error inflation factor, preferring epochs without the transient.

    Falls back to every epoch when the chosen phase does not cover enough
    wavelength bins to give a stable answer, since a factor measured from a
    handful of bins is worse than one contaminated by variability.

    :param photometry: Table of measurements, from measure_images
    :param edges: Wavelength bin edges in micron
    :return: Inflation factor
    """
    phase = "pre"
    if phase in set(photometry["phase"]):
        subset = photometry[photometry["phase"] == phase]
        n_bins = _populated_bins(subset, edges, CALIBRATION_MIN_PER_BIN)
        if n_bins >= MIN_CALIBRATION_BINS:
            logger.info(
                f"Calibrating errors on the {phase}-reference epochs, "
                f"where the transient is absent ({n_bins} bins)"
            )
            return error_inflation_factor(
                subset,
                edges=edges,
                min_per_bin=CALIBRATION_MIN_PER_BIN,
            )
        logger.info(
            f"Only {n_bins} usable {phase}-reference bins, "
            f"calibrating errors on every epoch instead"
        )

    return error_inflation_factor(photometry, edges=edges)


def _populated_bins(
    photometry: pd.DataFrame, edges: np.ndarray | None, min_per_bin: int
) -> int:
    """
    Count wavelength bins holding at least min_per_bin measurements.

    :param photometry: Table of measurements
    :param edges: Wavelength bin edges in micron
    :param min_per_bin: Threshold
    :return: Number of bins above the threshold
    """
    if edges is None:
        edges = default_wavelength_edges()
    if len(photometry) == 0:
        return 0

    counts = np.bincount(
        np.clip(np.digitize(photometry["wavelength"], edges) - 1, 0, len(edges) - 2),
        minlength=len(edges) - 1,
    )
    return int((counts >= min_per_bin).sum())


def difference_stacks(stacks: dict[str, pd.DataFrame]) -> pd.DataFrame | None:
    """
    Subtract the pre-reference baseline from each later block.

    Any transient flux shows up here, and a source that evolves shows a
    different difference in each block. The errors combined are the calibrated
    ones; 'flux_err_formal' and 'flux_scatter' are kept in the stacks if the raw
    or empirical versions are wanted instead.

    :param stacks: Mapping of block label to stack, from stack_by_block
    :return: Difference per wavelength bin per block, or None without a baseline
    """
    if "pre" not in stacks:
        logger.warning("No pre-reference baseline, so nothing to difference")
        return None

    pre = stacks["pre"]

    frames = []
    for label, post in stacks.items():
        if label == "pre":
            continue

        merged = pre.merge(post, on="wavelength", suffixes=("_pre", "_post"))
        merged["block"] = label
        merged["phase_days"] = post["phase_days"].iloc[0] if len(post) else np.nan
        merged["diff"] = merged["flux_post"] - merged["flux_pre"]
        merged["diff_err"] = np.hypot(merged["flux_err_pre"], merged["flux_err_post"])
        merged["sigma"] = merged["diff"] / merged["diff_err"]

        # If the differences really are noise, this is 1. Much above it means
        # the errors are still understated, or the block genuinely differs.
        rms_sigma = float(np.sqrt(np.mean(merged["sigma"] ** 2)))
        logger.info(f"{label}: rms(sigma) = {rms_sigma:.2f} over {len(merged)} bins")

        frames.append(merged)

    if len(frames) == 0:
        logger.warning("Only a pre-reference baseline, so nothing to difference")
        return None

    return pd.concat(frames, ignore_index=True)[
        [
            "block",
            "phase_days",
            "wavelength",
            "n_pre",
            "n_post",
            "flux_pre",
            "flux_post",
            "diff",
            "diff_err",
            "sigma",
        ]
    ]
