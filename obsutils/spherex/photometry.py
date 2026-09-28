"""
Forced PSF photometry on SPHEREx spectral image cutouts.

For each cutout the spatial WCS gives the target's pixel position, the
alternative 'W' WCS gives the wavelength that exposure happened to sample
there, and the FLAGS bitmap marks unusable pixels. Surface brightness is
converted to a flux density with the per-detector solid angle pixel map, the
effective PSF is rendered at the target's sub-pixel position, and its amplitude
is fitted together with a smooth background.

This is the method IRSA's Spectrophotometry Tool uses. Against its published
measurements the fluxes here agree to a median ratio of 1.05 and the errors to
3%; see notebooks/README_spherex_validation.md.
"""

import concurrent.futures
import contextlib
import http.client
import logging
import time
import urllib.error

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.units import Unit
from astropy.utils.data import conf
from astropy.wcs import WCS

from obsutils.spherex.calibration import load_solid_angle_maps
from obsutils.spherex.constants import (
    BAD_FLAGS,
    DEBLEND_SKY_ORDER,
    FLAGS_EXTENSION,
    IMAGE_EXTENSION,
    PSF_STAMP,
    VARIANCE_EXTENSION,
)
from obsutils.spherex.psf import fit_psf, load_epsf, render_epsf

logger = logging.getLogger(__name__)

# SPHEREx MEFs are large enough that astropy's default read timeout is too short.
REMOTE_TIMEOUT = 120

MAX_RETRIES = 4
MAX_WORKERS = 8

# A fit whose stamp was clipped by the cutout edge saw only part of the PSF.
MIN_PSF_ENCLOSED = 0.99

# Errors IRSA intermittently returns when many cutouts are requested at once.
TRANSIENT_ERRORS = (
    TimeoutError,
    urllib.error.HTTPError,
    urllib.error.URLError,
    http.client.IncompleteRead,
    OSError,
)


@contextlib.contextmanager
def quiet_astropy():
    """
    Temporarily silence astropy's logger and lengthen its read timeout.

    Every SPHEREx image carries alternative WCSs, about which astropy logs at
    INFO for each one read. Restores the previous settings on exit.
    """
    astropy_logger = logging.getLogger("astropy")
    previous_level = astropy_logger.level
    previous_timeout = conf.remote_timeout

    astropy_logger.setLevel(logging.ERROR)
    conf.remote_timeout = REMOTE_TIMEOUT
    try:
        yield
    finally:
        astropy_logger.setLevel(previous_level)
        conf.remote_timeout = previous_timeout


def bad_pixel_mask(
    flux_image: np.ndarray, flags: np.ndarray, flags_header: fits.Header
) -> np.ndarray:
    """
    Build a mask of pixels that should not be used.

    The bit assignments are read from the FLAGS header itself, as MP_* keywords,
    rather than assumed, since they can change between releases.

    :param flux_image: Image data, used to mask non-finite pixels
    :param flags: FLAGS extension data
    :param flags_header: FLAGS extension header, holding the MP_* bit definitions
    :return: Boolean mask, True where a pixel is unusable
    """
    bits = {key[3:]: flags_header[key] for key in flags_header if key.startswith("MP_")}

    mask = ~np.isfinite(flux_image)
    for name in BAD_FLAGS:
        if name in bits:
            mask |= (flags & (1 << bits[name])) > 0

    return mask


def surface_brightness_to_flux(
    header: fits.Header,
    shape: tuple[int, int],
    solid_angle_map: np.ndarray,
    solid_angle_unit: Unit,
) -> np.ndarray:
    """
    Build the per-pixel factor converting the IMAGE extension to uJy.

    CRPIX1A/CRPIX2A give the cutout's position on the parent frame, so the
    matching part of the solid angle map can be looked up.

    :param header: IMAGE extension header of the cutout
    :param shape: Shape of the cutout
    :param solid_angle_map: Solid angle pixel map for this detector
    :param solid_angle_unit: Unit of the solid angle pixel map
    :return: Array of the same shape as the cutout, to multiply the image by
    """
    n_y, n_x = shape
    y_grid, x_grid = np.mgrid[0:n_y, 0:n_x]

    parent_x = np.clip(
        np.rint(1 + x_grid - header["CRPIX1A"]).astype(int),
        0,
        solid_angle_map.shape[1] - 1,
    )
    parent_y = np.clip(
        np.rint(1 + y_grid - header["CRPIX2A"]).astype(int),
        0,
        solid_angle_map.shape[0] - 1,
    )

    pixel_area = (solid_angle_map[parent_y, parent_x] * solid_angle_unit).to(
        u.arcsec**2
    )

    return (1 * Unit(header["BUNIT"])).to(u.uJy / u.arcsec**2).value * pixel_area.value


def measure_cutout(
    row: pd.Series, coord: SkyCoord, companions: list[SkyCoord] | None = None
) -> dict | None:
    """
    Measure a target in one SPHEREx spectral image cutout.

    :param row: Row from query_images, providing 'cutout_url' and 'release'
    :param coord: Position of the target
    :param companions: Nearby sources to fit simultaneously, so their light is
        not absorbed into the target's amplitude. Passing any also tilts the
        fitted background, which a flat one cannot do for an extended host.
    :return: Epoch, sampled wavelength, and fitted fluxes in uJy, or None if the
        target is not measurable on this image
    """
    sky_order = DEBLEND_SKY_ORDER if companions else 0
    with fits.open(row["cutout_url"], cache=True) as hdul:
        header = hdul[IMAGE_EXTENSION].header
        detector = int(header["DETECTOR"])
        image = hdul[IMAGE_EXTENSION].data.astype(float)
        variance = hdul[VARIANCE_EXTENSION].data.astype(float)

        spatial_wcs = WCS(header)
        x, y = spatial_wcs.world_to_pixel(coord)

        # Wavelength and bandwidth sampled at that position. The spectral WCS is
        # a WAVE-TAB lookup table in the WCS-WAVE extension, and the SIP
        # distortion of the spatial WCS does not apply to it.
        spectral_wcs = WCS(header, fobj=hdul, key="W")
        spectral_wcs.sip = None
        wavelength, bandwidth = spectral_wcs.pixel_to_world(x, y)

        solid_angle_map, solid_angle_unit = load_solid_angle_maps(row["release"])[
            detector
        ]
        scale = surface_brightness_to_flux(
            header, image.shape, solid_angle_map, solid_angle_unit
        )
        flux_image = image * scale
        variance_image = variance * scale**2

        mask = bad_pixel_mask(
            flux_image, hdul[FLAGS_EXTENSION].data, hdul[FLAGS_EXTENSION].header
        )

        # The ePSF varies across the array, so the zone is chosen from the
        # position on the parent frame rather than on the cutout. Companions get
        # a wider stamp, since they sit off to one side of the fitting region.
        psfs = []
        for index, position in enumerate([coord] + list(companions or [])):
            source_x, source_y = (
                (x, y) if index == 0 else spatial_wcs.world_to_pixel(position)
            )
            psfs.append(
                render_epsf(
                    detector,
                    (
                        1 + float(source_x) - header["CRPIX1A"],
                        1 + float(source_y) - header["CRPIX2A"],
                    ),
                    (float(source_x), float(source_y)),
                    stamp_size=PSF_STAMP if index == 0 else PSF_STAMP + 6,
                )
            )

    fit = fit_psf(flux_image, variance_image, mask, psfs, sky_order=sky_order)
    if fit is None:
        # Usually a cutout truncated at the edge of the array, leaving the
        # target outside it. Nothing to measure, so skip the image.
        logger.debug(f"No usable fit for {row['uri']}")
        return None

    stamp, x_ll, y_ll = psfs[0]
    footprint = mask[
        max(y_ll, 0) : y_ll + stamp.shape[0], max(x_ll, 0) : x_ll + stamp.shape[1]
    ]

    return {
        "mjd": header["MJD-OBS"],
        "detector": detector,
        "bandpass": row.get("bandpass"),
        "release": row["release"],
        "phase": row.get("phase"),
        "wavelength": wavelength.to(u.micron).value,
        "bandwidth": bandwidth.to(u.micron).value,
        "x": float(x),
        "y": float(y),
        "n_bad": float(footprint.sum()),
        "uri": row["uri"],
        **fit,
    }


def measure_images(
    images: pd.DataFrame,
    coord: SkyCoord,
    companions: list[SkyCoord] | None = None,
) -> pd.DataFrame:
    """
    Measure a target in every spectral image covering it.

    Cutouts are downloaded in parallel, since the work is dominated by IO, and
    the transient read errors IRSA returns are retried. Images that still could
    not be read are skipped, and logged.

    :param images: Table of images, from query_images
    :param coord: Position of the target
    :param companions: Nearby sources to fit simultaneously
    :return: One row per successfully measured image, sorted by epoch
    """
    start = time.time()

    # Preload the calibration so the worker threads share one cached copy
    # rather than racing to download it.
    for release in images["release"].unique():
        load_solid_angle_maps(release)
    for detector in sorted(images["bandpass"].str[-1].astype(int).unique()):
        load_epsf(detector)

    def measure(row: pd.Series) -> dict | None:
        for attempt in range(MAX_RETRIES):
            try:
                return measure_cutout(row, coord, companions=companions)
            except TRANSIENT_ERRORS as exc:
                if attempt == MAX_RETRIES - 1:
                    logger.warning(f"Giving up on {row['uri']}: {exc}")
                    return None
                time.sleep(5 * (attempt + 1))
        return None

    with quiet_astropy():
        with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            results = list(executor.map(measure, (row for _, row in images.iterrows())))

    records = [result for result in results if result is not None]
    logger.info(
        f"Measured {len(records)}/{len(images)} images in {time.time() - start:.0f} s"
    )
    if len(records) < len(images):
        logger.info(
            f"{len(images) - len(records)} images were skipped, either unreadable "
            f"or with the target off the edge of the cutout"
        )

    if len(records) == 0:
        return pd.DataFrame(records)

    return pd.DataFrame(records).sort_values("mjd").reset_index(drop=True)


def clean(photometry: pd.DataFrame) -> pd.DataFrame:
    """
    Drop unusable measurements.

    Flagged pixels are already excluded from the fit itself, so a measurement is
    dropped only when the stamp ran off the edge of the cutout, or the fit
    failed to return a finite flux.

    :param photometry: Table of measurements, from measure_images
    :return: Subset that passed
    """
    mask = np.isfinite(photometry["flux"]) & (photometry["flux_err"] > 0)
    mask &= photometry["psf_enclosed"] >= MIN_PSF_ENCLOSED

    logger.info(f"{mask.sum()}/{len(photometry)} measurements are usable")
    return photometry[mask].reset_index(drop=True)
