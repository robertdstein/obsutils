"""
SPHEREx calibration products, fetched from the public IRSA S3 bucket.
"""

import logging
import re
from functools import lru_cache

import numpy as np
import requests
from astropy.io import fits
from astropy.units import Unit

from obsutils.spherex.constants import S3_BUCKET

logger = logging.getLogger(__name__)

N_DETECTORS = 6

LISTING_TIMEOUT = 60  # seconds


def latest_calibration_collection(release: str, dataset_type: str) -> str:
    """
    Find the most recent calibration collection for a release.

    The bucket is listed over plain HTTPS, so no S3 client is needed.

    :param release: Quick release, e.g. 'qr3'
    :param dataset_type: Calibration product, e.g. 'solid_angle_pixel_map'
    :return: Collection name
    """
    response = requests.get(
        S3_BUCKET,
        params={
            "list-type": "2",
            "delimiter": "/",
            "prefix": f"{release}/{dataset_type}/",
        },
        timeout=LISTING_TIMEOUT,
    )
    response.raise_for_status()

    pattern = rf"<Prefix>{release}/{dataset_type}/([^<]+?)/</Prefix>"
    collections = sorted(re.findall(pattern, response.text))

    if len(collections) == 0:
        raise ValueError(f"No {dataset_type} products found for release {release}")

    return collections[-1]


@lru_cache(maxsize=None)
def load_solid_angle_maps(
    release: str, collection: str | None = None
) -> dict[int, tuple[np.ndarray, Unit]]:
    """
    Download the solid angle pixel map for each detector.

    The IMAGE extension of a spectral image is a surface brightness in MJy/sr,
    so it must be multiplied by the per-pixel solid angle to become a flux
    density. The pixel scale is not quite constant across an array, so these
    maps are preferable to assuming a fixed pixel size everywhere.

    Results are cached for the lifetime of the process.

    :param release: Quick release, e.g. 'qr3'
    :param collection: Calibration collection. Default is the most recent.
    :return: Mapping of detector number to (solid angle map, its unit)
    """
    if collection is None:
        collection = latest_calibration_collection(release, "solid_angle_pixel_map")

    maps = {}

    for detector in range(1, N_DETECTORS + 1):
        url = (
            f"{S3_BUCKET}/{release}/solid_angle_pixel_map/{collection}/{detector}/"
            f"solid_angle_pixel_map_D{detector}_spx_{collection}.fits"
        )
        with fits.open(url, cache=True) as hdul:
            maps[detector] = (
                hdul["IMAGE"].data.astype(float),
                Unit(hdul["IMAGE"].header["BUNIT"]),
            )

    logger.info(f"Loaded solid angle pixel maps for {release} ({collection})")
    return maps
