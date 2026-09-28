"""
Find SPHEREx spectral images covering a transient position.

SPHEREx images the sky through linear variable filters, so a single exposure
measures one wavelength at a given sky position, and the sampled wavelength
changes from visit to visit as the spacecraft repoints. There is therefore no
ready-made spectrum to download for a transient: the spectrum is assembled by
measuring the source in every spectral image that covers it.
"""

import logging

import numpy as np
import pandas as pd
import pyvo
from astropy.coordinates import SkyCoord

from obsutils.spherex.constants import (
    COMPANION_RADIUS,
    CUTOUT_SIZE,
    IRSA_ROOT,
    TAP_URL,
)
from obsutils.spherex.psf import psf_weight

# Catalogue searched for neighbours, and the limits applied to it. AllWISE is
# the only all-sky catalogue covering the SPHEREx wavelength range.
COMPANION_CATALOGUE = "allwise_p3as_psd"
SELF_MATCH_ARCSEC = 3.0
MAX_COMPANIONS = 5

# WISE W1 Vega zero point in Jy, and the wavelength its PSF overlap is
# evaluated at, for ranking how much each neighbour contaminates the target.
W1_ZERO_POINT_JY = 309.54
COMPANION_RANKING_WAVELENGTH = 3.353

logger = logging.getLogger(__name__)


def resolve_position(
    name: str | None = None,
    ra: float | None = None,
    dec: float | None = None,
) -> SkyCoord:
    """
    Resolve a transient to a sky position.

    Explicit coordinates take precedence. Otherwise the name is resolved
    through Sesame, which knows TNS names such as 'AT2025abcr'.

    :param name: Transient name to resolve
    :param ra: Right ascension in degrees, overriding the name
    :param dec: Declination in degrees, overriding the name
    :return: Position of the transient
    """
    if (ra is not None) & (dec is not None):
        return SkyCoord(ra, dec, unit="deg")

    if name is None:
        raise ValueError("Provide either a name or both ra and dec")

    coord = SkyCoord.from_name(name)
    logger.info(f"Resolved {name} to {coord.to_string('hmsdms')}")
    return coord


def query_images(coord: SkyCoord, reference_mjd: float | None = None) -> pd.DataFrame:
    """
    Find every SPHEREx Level 2 spectral image whose footprint contains a position.

    Uses TAP rather than the Simple Image Access service because TAP sees data
    as soon as it is ingested, while SIA lags by around a day.

    The returned 'cutout_url' points at a cutout of each image rather than the
    full 2040 x 2040 frame, which is what makes measuring hundreds of images
    practical.

    :param coord: Position to search
    :param reference_mjd: Epoch used to label rows as 'pre' or 'post', e.g. discovery
    :return: One row per spectral image
    """
    ra, dec = coord.ra.deg, coord.dec.deg

    query = f"""
    SELECT
        '{IRSA_ROOT}' || a.uri
            || '?center={ra},{dec}d&size={CUTOUT_SIZE}' AS cutout_url,
        a.uri AS uri,
        p.energy_bandpassname AS bandpass,
        p.time_bounds_lower AS mjd,
        p.time_exposure AS exptime,
        p.energy_bounds_lower AS wave_lower,
        p.energy_bounds_upper AS wave_upper,
        p.quality_flag AS quality_flag
    FROM spherex.artifact a
    JOIN spherex.plane p ON a.planeid = p.planeid
    WHERE 1 = CONTAINS(POINT('ICRS', {ra}, {dec}), p.poly)
      AND a.producttype = 'science'
    ORDER BY p.time_bounds_lower
    """

    service = pyvo.dal.TAPService(TAP_URL)
    images = service.search(query).to_table().to_pandas()

    if len(images) == 0:
        logger.warning(f"No SPHEREx images cover {coord.to_string('hmsdms')}")
        return images

    # The quick release an image belongs to is encoded in its artifact path,
    # and determines which calibration products apply to it.
    images["release"] = [
        uri.split("/spherex/")[1].split("/")[0] for uri in images["uri"]
    ]

    images["phase"] = label_phase(images["mjd"], reference_mjd)

    logger.info(
        f"Found {len(images)} SPHEREx images, "
        f"MJD {images['mjd'].min():.1f} to {images['mjd'].max():.1f}"
    )
    return images


def label_phase(mjd, reference_mjd: float | None) -> np.ndarray:
    """
    Label each epoch as before or after a reference epoch.

    :param mjd: Observation times
    :param reference_mjd: Reference epoch, e.g. discovery. If None, everything is 'all'.
    :return: Array of 'pre'/'post' labels, or 'all' if no reference epoch is given
    """
    if reference_mjd is None:
        return np.full(len(mjd), "all")
    return np.where(np.asarray(mjd) < reference_mjd, "pre", "post")


def summarise_coverage(images: pd.DataFrame) -> pd.DataFrame:
    """
    Summarise how many images cover the position, per detector and phase.

    :param images: Table of images, from query_images
    :return: Counts and time range per bandpass and phase
    """
    return (
        images.groupby(["bandpass", "phase"])
        .agg(n=("mjd", "size"), first_mjd=("mjd", "min"), last_mjd=("mjd", "max"))
        .unstack("phase")
    )


def find_companions(coord: SkyCoord) -> list[SkyCoord]:
    """
    Find the catalogued sources that most contaminate a target's photometry.

    At 6.15 arcsec/pixel a neighbour tens of arcsec away is only a pixel or two
    from the target, so its light lands in the same PSF. Fitting it explicitly
    keeps that light out of the target's amplitude.

    Sources are ranked by how much they actually contribute at the target
    position, which is their flux times the PSF overlap at their separation,
    not by brightness. A bright source 35 arcsec away contributes nothing, while
    a faint one at 10 arcsec contributes a great deal, so ranking on brightness
    fits the wrong ones. Anything closer than the PSF can resolve is the target
    itself and is left alone.

    This finds point sources. It does not solve the harder problem of an
    extended host, whose light no point-source model removes.

    :param coord: Position of the target
    :return: Companion positions, most contaminating first
    """
    ra, dec = coord.ra.deg, coord.dec.deg

    query = f"""
    SELECT ra, dec, w1mpro
    FROM {COMPANION_CATALOGUE}
    WHERE CONTAINS(POINT('ICRS', ra, dec),
                   CIRCLE('ICRS', {ra}, {dec}, {COMPANION_RADIUS / 3600})) = 1
    """
    table = pyvo.dal.TAPService(TAP_URL).search(query).to_table()

    if len(table) == 0:
        logger.info("No catalogued companions near the target")
        return []

    positions = SkyCoord(table["ra"], table["dec"], unit="deg")
    separation = coord.separation(positions).arcsec

    magnitude = np.asarray(table["w1mpro"], dtype=float)
    flux = W1_ZERO_POINT_JY * 10 ** (-magnitude / 2.5)
    contribution = flux * psf_weight(separation, COMPANION_RANKING_WAVELENGTH)
    contribution[~np.isfinite(contribution)] = 0.0

    # Anything this close cannot be separated from the target by the PSF, so it
    # is the target as far as the fit is concerned.
    contribution[separation < SELF_MATCH_ARCSEC] = 0.0

    order = np.argsort(-contribution)[:MAX_COMPANIONS]
    order = [index for index in order if contribution[index] > 0]

    companions = [positions[index] for index in order]
    separations = ", ".join(f"{separation[index]:.1f}" for index in order)
    logger.info(
        f"Fitting {len(companions)} companions within {COMPANION_RADIUS:.0f} arcsec, "
        f"at {separations} arcsec"
    )

    return companions
