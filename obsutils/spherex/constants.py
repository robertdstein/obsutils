"""
Constants describing the SPHEREx archive at IRSA.
"""

import numpy as np

# IRSA Table Access Protocol endpoint, used to find spectral images.
TAP_URL = "https://irsa.ipac.caltech.edu/TAP"

# Root of the IRSA data server, prepended to the artifact URIs returned by TAP.
IRSA_ROOT = "https://irsa.ipac.caltech.edu/"

# Public S3 bucket holding the SPHEREx calibration products.
S3_BUCKET = "https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com"

# Nominal SPHEREx pixel scale, in arcsec.
PIXEL_SCALE = 6.15

# Wavelength range covered by the six detectors, in micron.
WAVELENGTH_MIN = 0.74
WAVELENGTH_MAX = 5.05

# Wavelength range of each detector, in micron, from the SPHEREx
# Explanatory Supplement.
DETECTOR_WAVELENGTHS = (
    (1, 0.75, 1.09),
    (2, 1.10, 1.62),
    (3, 1.63, 2.41),
    (4, 2.42, 3.82),
    (5, 3.83, 4.41),
    (6, 4.42, 5.00),
)

# SPHEREx completes a full-sky survey every six months, so a position is
# revisited on that cadence, in passes that span a few weeks. Measurements
# separated by more than a quarter of that period are from different passes:
# in practice the gaps within a pass run to about 40 days and those between
# passes are never shorter than 120, so the split is insensitive to the exact
# fraction chosen.
SURVEY_PERIOD_DAYS = 182.6
BLOCK_GAP_DAYS = SURVEY_PERIOD_DAYS / 4

# Cutout half-size in degrees, comfortably larger than the fitting stamp.
CUTOUT_SIZE = 0.05

# Width of the PSF fitting stamp, in native pixels. Odd values keep the source
# centred. The ePSF support is under 7 pixels across.
PSF_STAMP = 11

# Radius in arcsec within which catalogued neighbours are fitted alongside the
# target when deblending.
COMPANION_RADIUS = 40.0

# Background polynomial order used when deblending. A flat background cannot
# follow an extended host across the fitting stamp; a tilted plane can.
DEBLEND_SKY_ORDER = 1

# Release holding the ePSF calibration products. QR3 is the only one that has
# them; earlier releases ship a superseded PSF, so these are used for all data.
EPSF_RELEASE = "qr3"

# FLAGS bits treated as unusable. Full definitions are in Table 8 of the
# SPHEREx Explanatory Supplement, and in each FLAGS header as MP_* keywords.
BAD_FLAGS = (
    "NONFUNC",
    "HOT",
    "DARKFRAME",
    "NONLINEAR",
    "SATURATE",
    "MISSING",
    "TRANSIENT",
    "COSMICRAY",
    "MASKED",
)

# Names of the MEF extensions used by the photometry.
IMAGE_EXTENSION = "IMAGE"
VARIANCE_EXTENSION = "VARIANCE"
FLAGS_EXTENSION = "FLAGS"


def detector_for_wavelength(wavelength: float) -> int:
    """
    Which detector covers a wavelength.

    :param wavelength: Wavelength in micron
    :return: Detector number, 1 to 6
    """
    for detector, lower, upper in DETECTOR_WAVELENGTHS:
        if lower <= wavelength <= upper:
            return detector

    return 1 if wavelength < WAVELENGTH_MIN else 6


def default_wavelength_edges(n_bins: int = 30) -> np.ndarray:
    """
    Logarithmically-spaced wavelength bin edges spanning the SPHEREx range.

    :param n_bins: Number of bins
    :return: Array of n_bins + 1 bin edges, in micron
    """
    return np.geomspace(WAVELENGTH_MIN, WAVELENGTH_MAX, n_bins + 1)
