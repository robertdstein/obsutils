"""
Forced PSF photometry on SPHEREx spectral images.

This is the method IRSA's Spectrophotometry Tool uses: the source position is
held fixed, the effective PSF (ePSF) is rendered at that sub-pixel position,
and only the flux amplitude is fitted. It is the right choice for a point
source, and it is far less sensitive to neighbours and structured background
than a fixed aperture.

The ePSF is a sub-pixel-resolved model of the pixel-convolved PSF, delivered as
a per-detector calibration product covering the array in zones. The models
stored inside the spectral images themselves are superseded for anything before
QR3, so the calibration products are used here instead, for every release.
"""

import logging
import re
from functools import lru_cache

import numpy as np
import requests
from astropy.io import fits
from astropy.table import Table
from photutils.psf import ImagePSF

from obsutils.spherex.constants import (
    EPSF_RELEASE,
    PIXEL_SCALE,
    PSF_STAMP,
    S3_BUCKET,
    detector_for_wavelength,
)

logger = logging.getLogger(__name__)

LISTING_TIMEOUT = 60  # seconds

# Minimum usable pixels for a fit to be attempted.
MIN_FIT_PIXELS = 10

# A component whose peak in the fitting region is below this contributes
# nothing, and would make the solve singular.
NEGLIGIBLE_COMPONENT = 1e-6

EPSF_EXTENSION = "EPSF"


@lru_cache(maxsize=None)
def epsf_urls() -> dict[int, str]:
    """
    Find the newest ePSF calibration file for each detector.

    Detectors are versioned independently, so the newest collection overall is
    not necessarily the newest for a given detector.

    :return: Mapping of detector number to file URL
    """
    response = requests.get(
        S3_BUCKET,
        params={"list-type": "2", "prefix": f"{EPSF_RELEASE}/epsf/"},
        timeout=LISTING_TIMEOUT,
    )
    response.raise_for_status()

    newest: dict[int, tuple[int, str]] = {}
    for key in re.findall(r"<Key>([^<]+)</Key>", response.text):
        match = re.search(r"epsf_D(\d)_spx_(cal-epsf-v(\d+)-[\d-]+)\.fits", key)
        if match is None:
            continue
        detector, version = int(match.group(1)), int(match.group(3))
        if detector not in newest or version > newest[detector][0]:
            newest[detector] = (version, f"{S3_BUCKET}/{key}")

    if len(newest) == 0:
        raise ValueError(f"No ePSF calibration products found in {EPSF_RELEASE}")

    return {detector: url for detector, (_, url) in newest.items()}


@lru_cache(maxsize=None)
def load_epsf(detector: int) -> tuple[Table, tuple[int, int]]:
    """
    Load the ePSF zone table for one detector.

    :param detector: Detector number, 1 to 6
    :return: Zone table, and the (y, x) oversampling of the stored models
    """
    url = epsf_urls()[detector]

    with fits.open(url, cache=True) as hdul:
        table = Table(hdul[EPSF_EXTENSION].data)
        header = hdul[EPSF_EXTENSION].header

    oversampling = (int(header.get("OVSMPY", 1)), int(header.get("OVSMPX", 1)))
    logger.info(
        f"Loaded ePSF for detector {detector} "
        f"({url.rsplit('/', maxsplit=1)[-1]}, {len(table)} zones)"
    )

    return table, oversampling


def render_epsf(
    detector: int,
    parent: tuple[float, float],
    position: tuple[float, float],
    stamp_size: int = PSF_STAMP,
) -> tuple[np.ndarray, int, int]:
    """
    Render the ePSF onto the native pixel grid at a sub-pixel position.

    The zone is chosen from the position on the parent array, since the PSF
    varies across the detector. The stamp is then evaluated at the source's
    actual sub-pixel phase and integrated onto native pixels, which matters
    because the SPHEREx PSF is undersampled.

    :param detector: Detector number, 1 to 6
    :param parent: Source (x, y) on the parent array, for choosing the zone
    :param position: Source (x, y) in the frame the stamp will be used in
    :param stamp_size: Width of the rendered stamp, in native pixels
    :return: Unit-sum PSF stamp, and the x and y index of its lower-left pixel
    """
    table, oversampling = load_epsf(detector)

    zone = int(
        np.argmin(np.hypot(parent[0] - table["XCENTER"], parent[1] - table["YCENTER"]))
    )
    model = ImagePSF(
        data=np.asarray(table["EPSF"][zone], dtype=float), oversampling=oversampling
    )

    x, y = position
    half = (stamp_size - 1) / 2.0
    x_ll, y_ll = int(np.floor(x - half)), int(np.floor(y - half))
    grid_x, grid_y = np.meshgrid(
        x_ll + np.arange(stamp_size), y_ll + np.arange(stamp_size), indexing="xy"
    )

    stamp = model.evaluate(grid_x, grid_y, flux=1.0, x_0=x, y_0=y)
    total = stamp.sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError(f"Degenerate ePSF stamp for detector {detector}")

    return stamp / total, x_ll, y_ll


def fit_psf(
    flux_image: np.ndarray,
    variance_image: np.ndarray,
    mask: np.ndarray,
    psfs: list[tuple[np.ndarray, int, int]],
    sky_order: int = 0,
) -> dict | None:
    """
    Fit one or more fixed-position PSFs plus a smooth background to a cutout.

    The model is linear in every amplitude, so this is a single weighted
    least-squares solve, with each pixel weighted by its inverse variance.
    The fitting region is set by the first PSF, which is the target; any
    further PSFs are companions fitted simultaneously so that their light is
    not absorbed into the target's amplitude.

    A flat background cannot represent an extended host galaxy. Raising
    sky_order adds a tilted plane (1) or a quadratic surface (2), which helps,
    but no smooth term removes host light that is concentrated on the scale of
    the PSF itself. On an extended host, read the difference between epochs
    rather than the absolute amplitude.

    :param flux_image: Cutout in flux units
    :param variance_image: Matching variance image, in those units squared
    :param mask: Boolean mask, True where a pixel is unusable
    :param psfs: (stamp, x_ll, y_ll) per source, from render_epsf, target first
    :param sky_order: Polynomial order of the background, 0, 1 or 2
    :return: Fitted fluxes, errors, background and diagnostics, or None
    """
    if len(psfs) == 0:
        raise ValueError("Need at least one PSF to fit")
    if sky_order not in (0, 1, 2):
        raise ValueError(f"sky_order must be 0, 1 or 2, got {sky_order}")

    n_y, n_x = flux_image.shape
    target, x_ll, y_ll = psfs[0]
    stamp_y, stamp_x = target.shape

    # The fitting region is the target's stamp, clipped to the cutout.
    x_start, x_stop = max(x_ll, 0), min(x_ll + stamp_x, n_x)
    y_start, y_stop = max(y_ll, 0), min(y_ll + stamp_y, n_y)
    if (x_stop <= x_start) or (y_stop <= y_start):
        return None

    data = flux_image[y_start:y_stop, x_start:x_stop]
    variance = variance_image[y_start:y_stop, x_start:x_stop]
    bad = mask[y_start:y_stop, x_start:x_stop]

    good = (~bad) & np.isfinite(data) & np.isfinite(variance) & (variance > 0)
    if good.sum() < MIN_FIT_PIXELS + 2 * sky_order:
        return None

    region = (x_start, x_stop, y_start, y_stop)
    components = [
        _place(stamp, stamp_x_ll, stamp_y_ll, region)
        for stamp, stamp_x_ll, stamp_y_ll in psfs
    ]

    # A companion far enough out that none of its PSF reaches the fitting
    # region contributes an all-zero column, which makes the solve singular.
    # Its amplitude is unconstrained here anyway, so drop it. The target is
    # index 0 and is always kept, so the reported flux stays comparable.
    contributing = [0] + [
        index
        for index in range(1, len(components))
        if components[index][good].max() > NEGLIGIBLE_COMPONENT
    ]
    components = [components[index] for index in contributing]

    # Background terms, in coordinates centred on the fitting region.
    grid_y, grid_x = np.mgrid[y_start:y_stop, x_start:x_stop]
    offset_x = grid_x[good] - 0.5 * (x_start + x_stop)
    offset_y = grid_y[good] - 0.5 * (y_start + y_stop)
    background = [np.ones(good.sum())]
    if sky_order >= 1:
        background += [offset_x, offset_y]
    if sky_order >= 2:
        background += [offset_x**2, offset_y**2, offset_x * offset_y]

    design = np.column_stack([c[good] for c in components] + background)
    weights = 1.0 / variance[good]
    weighted = design.T * weights

    try:
        covariance = np.linalg.inv(weighted @ design)
    except np.linalg.LinAlgError:
        return None

    parameters = covariance @ (weighted @ data[good])
    residual = data[good] - design @ parameters
    errors = np.sqrt(np.diag(covariance))

    result = {
        "flux": float(parameters[0]),
        "flux_err": float(errors[0]),
        "sky": float(parameters[len(components)]),
        "chi2": float(
            np.sum(residual**2 * weights) / max(good.sum() - design.shape[1], 1)
        ),
        "n_fit_pixels": int(good.sum()),
        "n_sources": len(components),
        # Below 1 when the stamp was clipped by an edge, so the fit saw only
        # part of the target's PSF.
        "psf_enclosed": float(components[0][good].sum()),
    }

    for fitted, original in enumerate(contributing):
        if original == 0:
            continue
        result[f"flux_companion{original}"] = float(parameters[fitted])
        result[f"flux_err_companion{original}"] = float(errors[fitted])
        # Strong anti-correlation means the fit cannot tell the two apart.
        result[f"corr_companion{original}"] = float(
            covariance[0, fitted]
            / np.sqrt(covariance[0, 0] * covariance[fitted, fitted])
        )

    return result


def _place(
    stamp: np.ndarray, x_ll: int, y_ll: int, region: tuple[int, int, int, int]
) -> np.ndarray:
    """
    Drop a rendered stamp into the fitting region, zero-padded.

    A companion can sit partly or wholly outside the region, so only the
    overlap is copied.

    :param stamp: Rendered PSF stamp
    :param x_ll: x index of the stamp's lower-left pixel in the cutout
    :param y_ll: y index of the stamp's lower-left pixel in the cutout
    :param region: Fitting region as (x_start, x_stop, y_start, y_stop)
    :return: Array shaped like the fitting region
    """
    x_start, x_stop, y_start, y_stop = region
    placed = np.zeros((y_stop - y_start, x_stop - x_start))

    x_from, x_to = max(x_ll, x_start), min(x_ll + stamp.shape[1], x_stop)
    y_from, y_to = max(y_ll, y_start), min(y_ll + stamp.shape[0], y_stop)
    if (x_to <= x_from) or (y_to <= y_from):
        return placed

    placed[y_from - y_start : y_to - y_start, x_from - x_start : x_to - x_start] = (
        stamp[y_from - y_ll : y_to - y_ll, x_from - x_ll : x_to - x_ll]
    )
    return placed


# Separations, in arcsec, at which the crosstalk curve is tabulated.
CROSSTALK_SAMPLES = np.concatenate([np.arange(0, 12, 0.5), np.arange(12, 41, 2.0)])

# Grids used to build the crosstalk curve, in native pixels.
CROSSTALK_STAMP = 15
CROSSTALK_CANVAS = 61


@lru_cache(maxsize=None)
def _crosstalk_curve(detector: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Tabulate how much of a neighbour's flux a forced fit assigns to the target.

    For a fit that solves only for amplitude at a fixed position, a source
    offset by s contributes the normalised overlap of the two PSFs. The real
    ePSF has far broader wings than a Gaussian, so this is much larger at
    separations beyond a resolution element than a Gaussian would suggest.

    :param detector: Detector number, 1 to 6
    :return: Separations in arcsec, and the corresponding weights
    """
    centre = CROSSTALK_CANVAS / 2.0

    def rendered(x: float) -> np.ndarray:
        stamp, x_ll, y_ll = render_epsf(
            detector, (1000, 1000), (x, centre), stamp_size=CROSSTALK_STAMP
        )
        canvas = (0, CROSSTALK_CANVAS, 0, CROSSTALK_CANVAS)
        return _place(stamp, x_ll, y_ll, canvas)

    reference = rendered(centre)
    norm = float(np.sum(reference**2))

    weights = [
        float(np.sum(reference * rendered(centre + s / PIXEL_SCALE)) / norm)
        for s in CROSSTALK_SAMPLES
    ]
    # The overlap must fall off; enforce that against sub-pixel rendering noise.
    weights = list(np.minimum.accumulate(np.clip(weights, 0.0, 1.0)))

    return tuple(CROSSTALK_SAMPLES), tuple(weights)


def psf_weight(separation: np.ndarray, wavelength: float) -> np.ndarray:
    """
    Fraction of a source's flux that a forced fit at the target picks up.

    :param separation: Offsets from the target in arcsec
    :param wavelength: Wavelength in micron, setting which PSF applies
    :return: Weight between 0 and 1 per source
    """
    samples, weights = _crosstalk_curve(detector_for_wavelength(wavelength))

    return np.interp(
        np.asarray(separation, dtype=float),
        samples,
        weights,
        left=weights[0],
        right=0.0,
    )
