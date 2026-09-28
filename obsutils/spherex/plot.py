"""
Plots of SPHEREx coverage, light curves and stacked spectra.

Every function returns its figure, so callers can save it or adjust it further.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

logger = logging.getLogger(__name__)


# The baseline is drawn flat grey so the later blocks, which carry the colour
# ramp, read as the thing that changes.
BASELINE_COLOR = "0.45"
BLOCK_COLORMAP = "plasma"


def block_label(label: str, stack: pd.DataFrame) -> str:
    """
    Legend entry for one block, with its phase when that is known.

    :param label: Block label, from assign_blocks
    :param stack: That block's stack
    :return: Text for the legend
    """
    count = int(stack["n"].sum())
    phase = stack["phase_days"].iloc[0] if len(stack) else np.nan

    if not np.isfinite(phase):
        return f"{label} (n = {count})"

    return f"{label}, {phase:+.0f} d (n = {count})"


def block_colors(labels: list[str]) -> dict[str, str]:
    """
    Assign the baseline grey and spread a colour ramp over the later blocks.

    :param labels: Block labels in plotting order
    :return: Mapping of label to colour
    """
    later = [label for label in labels if label != "pre"]
    ramp = plt.get_cmap(BLOCK_COLORMAP)(np.linspace(0.05, 0.75, max(len(later), 1)))

    colors = {"pre": BASELINE_COLOR}
    colors.update({label: ramp[index] for index, label in enumerate(later)})
    return colors


def plot_coverage(
    images: pd.DataFrame, reference_mjd: float | None = None, title: str = ""
) -> Figure:
    """
    Show when each detector observed the position.

    :param images: Table of images, from query_images
    :param reference_mjd: Epoch to mark, e.g. discovery
    :param title: Name of the target, used in the title
    :return: The figure
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    for bandpass, group in images.groupby("bandpass"):
        ax.scatter(group["mjd"], [bandpass] * len(group), s=12)

    if reference_mjd is not None:
        ax.axvline(reference_mjd, color="k", ls="--", label=f"MJD {reference_mjd:.1f}")
        ax.legend()

    ax.set_xlabel("MJD")
    ax.set_title(f"{title} SPHEREx coverage ({len(images)} spectral images)".strip())
    fig.tight_layout()

    return fig


def plot_lightcurve(
    photometry: pd.DataFrame,
    reference_mjd: float | None = None,
    title: str = "",
) -> Figure:
    """
    Plot aperture flux against time, one panel per detector.

    Scatter within a panel is dominated by the fact that consecutive visits
    sample different wavelengths, not by variability, so points are coloured
    by the wavelength each exposure measured.

    :param photometry: Table of measurements, from measure_images
    :param reference_mjd: Epoch to mark, e.g. discovery
    :param title: Name of the target, used in the title
    :return: The figure
    """

    detectors = sorted(photometry["detector"].unique())
    n_rows = int(np.ceil(len(detectors) / 2))
    fig, axes = plt.subplots(
        n_rows, 2, figsize=(12, 3 * n_rows), sharex=True, squeeze=False
    )

    for ax, detector in zip(axes.ravel(), detectors):
        group = photometry[photometry["detector"] == detector]
        scatter = ax.scatter(
            group["mjd"],
            group["flux"] / 1e3,
            c=group["wavelength"],
            s=14,
            cmap="viridis",
        )
        ax.errorbar(
            group["mjd"],
            group["flux"] / 1e3,
            yerr=group["flux_err"] / 1e3,
            fmt="none",
            ecolor="grey",
            alpha=0.5,
            lw=0.8,
        )
        if reference_mjd is not None:
            ax.axvline(reference_mjd, color="k", ls="--", lw=1)
        ax.set_title(f"Detector {detector}")
        plt.colorbar(scatter, ax=ax, label="wavelength [um]")

    for ax in axes.ravel()[len(detectors) :]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("MJD")
    for ax in axes[:, 0]:
        ax.set_ylabel("PSF flux [mJy]")

    fig.suptitle(f"{title} SPHEREx PSF flux".strip())
    fig.tight_layout()

    return fig


def plot_sed(
    stacks: dict[str, pd.DataFrame],
    title: str = "",
    difference: pd.DataFrame | None = None,
    archival: pd.DataFrame | None = None,
) -> Figure:
    """
    Plot the wavelength-stacked spectrum, one line per block, with error bars.

    On a host-dominated source the phases lie on top of each other and the
    error bars are far smaller than the plotted range, so pass a difference
    table to add a residual panel: that is where the error bars can actually be
    judged against the epoch-to-epoch variation.

    :param stacks: Mapping of block label to stack, from stack_by_block
    :param title: Name of the target, used in the title
    :param difference: Difference table from difference_stacks, for the residual panel
    :param archival: Pre-SPHEREx estimate from archival_photometry, overlaid for
        comparison. It is the total catalogued flux in the PSF, so it lines up
        with a single-PSF measurement rather than with a joint fit that has
        removed companions.
    :return: The figure
    """
    if difference is None:
        fig, ax = plt.subplots(figsize=(10, 5))
        axes = [ax]
    else:
        fig, axes = plt.subplots(
            2,
            1,
            figsize=(10, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [2.5, 1]},
        )
        ax = axes[0]

    colors = block_colors(list(stacks))
    for label, stack in stacks.items():
        ax.errorbar(
            stack["wavelength"],
            stack["flux"] / 1e3,
            yerr=stack["flux_err"] / 1e3,
            color=colors[label],
            marker="o" if label == "pre" else "s",
            ms=4,
            ls="--" if label == "pre" else "-",
            capsize=3,
            label=block_label(label, stack),
        )

    if archival is not None and len(archival) > 0:
        ax.plot(
            archival["wavelength"],
            archival["flux"] / 1e3,
            ls=":",
            marker="D",
            ms=6,
            mfc="none",
            color="k",
            lw=1.2,
            label="estimated from archival photometry",
        )
        for _, row in archival.iterrows():
            ax.annotate(
                row["band"],
                (row["wavelength"], row["flux"] / 1e3),
                textcoords="offset points",
                xytext=(0, -14),
                ha="center",
                fontsize=7,
                color="k",
            )

    ax.set_xscale("log")
    ax.set_ylabel("PSF flux [mJy]")
    ax.set_title(f"{title} stacked SPHEREx SED".strip())
    ax.legend()

    if difference is not None:
        residual_ax = axes[1]
        for label, group in difference.groupby("block", sort=False):
            residual_ax.errorbar(
                group["wavelength"],
                group["diff"],
                yerr=group["diff_err"],
                fmt="o",
                ms=4,
                color=colors.get(label, "k"),
                capsize=3,
                label=label,
            )
        residual_ax.axhline(0, color="grey", lw=1)
        residual_ax.set_ylabel("block - pre [uJy]")
        rms_sigma = np.sqrt(np.mean(difference["sigma"] ** 2))
        residual_ax.set_title(
            f"residual, rms = {rms_sigma:.1f} sigma "
            f"(1 if the variation is just noise)",
            fontsize=10,
        )

    axes[-1].set_xlabel("observed wavelength [um]")
    fig.tight_layout()

    return fig


def plot_difference(difference: pd.DataFrame, title: str = "") -> Figure:
    """
    Plot each block's difference from the baseline, in flux and in sigma.

    The lower panel is the same points divided by their errors. If a block
    differs from the baseline only by noise, its points scatter within the
    shaded band.

    :param difference: Difference table, from difference_stacks
    :param title: Name of the target, used in the title
    :return: The figure
    """
    fig, axes = plt.subplots(
        2, 1, figsize=(10, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )

    labels = list(dict.fromkeys(difference["block"]))
    colors = block_colors(["pre"] + labels)

    for label, group in difference.groupby("block", sort=False):
        phase = group["phase_days"].iloc[0]
        legend = label if not np.isfinite(phase) else f"{label}, {phase:+.0f} d"
        axes[0].errorbar(
            group["wavelength"],
            group["diff"],
            yerr=group["diff_err"],
            fmt="o",
            ms=4,
            color=colors[label],
            capsize=3,
            label=legend,
        )
        axes[1].scatter(group["wavelength"], group["sigma"], color=colors[label], s=18)

    axes[0].axhline(0, color="grey", lw=1)
    axes[0].set_ylabel("block - pre flux [uJy]")
    axes[0].set_title(
        f"{title} flux difference from the pre-reference baseline".strip()
    )
    axes[0].legend(fontsize=9)

    axes[1].axhspan(-1, 1, color="grey", alpha=0.25, label="1 sigma")
    axes[1].axhline(0, color="grey", lw=1)
    axes[1].set_ylabel("sigma")
    axes[1].set_xlabel("observed wavelength [um]")
    axes[1].set_xscale("log")
    axes[1].legend(loc="upper right", fontsize=9)

    rms_sigma = np.sqrt(np.mean(difference["sigma"] ** 2))
    axes[1].set_title(f"rms = {rms_sigma:.2f} sigma", fontsize=10)

    fig.tight_layout()

    return fig
