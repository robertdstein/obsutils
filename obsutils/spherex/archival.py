"""
Predict the pre-SPHEREx flux at a position from archival catalogues.

A forced measurement at a position picks up every source close enough to share
the PSF, so its baseline is rarely zero. Summing catalogued sources, each
weighted by how much of its light a PSF centred on the target would collect,
gives an independent estimate of that baseline from data taken years before
SPHEREx launched. Overlaid on a SPHEREx SED it shows at a glance whether
anything has changed.

Extended sources need their extended photometry. A galaxy appears in the
point-source catalogues too, but only its core does, so using the point-source
magnitude leaves most of the galaxy out: the 2MASS point-source entry for the
AT2025abcr host is 4.1 mJy at Ks against 10.5 mJy for the whole galaxy. The
2MASS extended source catalogue is queried as well, and its magnitudes take
precedence where a source appears in both, and the WISE elliptical aperture
magnitudes are used for sources WISE flags as extended. That leaves 2MASS and
WISE covering 1.2 to 4.6 micron, which is most of the SPHEREx range.

Pan-STARRS is deliberately not used: neither of its magnitudes works for a
galaxy, the PSF ones measuring only the core and the Kron ones double counting,
since Pan-STARRS deblends an extended source into fragments that each carry a
near-total Kron magnitude.

Compare it against a single-PSF measurement, not a joint fit. The estimate is
the total catalogued flux falling in the PSF, which is what one PSF fitted at
the position measures. A fit that also solves for companions has deliberately
removed some of that, so the two are no longer the same quantity, and excluding
the companions from the estimate does not fix it: their catalogue entries carry
the extended light that physically remains at the target position.

The estimate is still approximate, which is why plots label it "estimated from
archival photometry": the catalogues are incomplete in crowded fields, an
extended galaxy is not really a point source however it is measured, and the
archival epochs span decades. On AT2025abcr it agrees with the single-PSF
SPHEREx photometry to 2-29% across 1.2 to 4.6 micron. In the crowded Galactic
plane field of V1935 Cen it over-predicts by a factor of two in the near
infrared and by much more at W1 and W2.
"""

import logging

import numpy as np
import pandas as pd
import pyvo
from astropy.coordinates import SkyCoord

from obsutils.spherex.constants import TAP_URL, WAVELENGTH_MAX, WAVELENGTH_MIN
from obsutils.spherex.psf import psf_weight

logger = logging.getLogger(__name__)

# Archival bands inside the SPHEREx range. Zero points are in Jy: 3631 for the
# AB magnitudes of Pan-STARRS, and the Vega zero points for 2MASS and WISE.
BANDS = (
    ("2MASS J", "2mass", "j_m", 1.235, 1594.0),
    ("2MASS H", "2mass", "h_m", 1.662, 1024.0),
    ("2MASS Ks", "2mass", "k_m", 2.159, 666.7),
    ("WISE W1", "wise", "w1mpro", 3.353, 309.54),
    ("WISE W2", "wise", "w2mpro", 4.603, 171.787),
)

# Extended-source equivalents, used where a source is resolved.
XSC_COLUMNS = {"j_m": "j_m_k20fe", "h_m": "h_m_k20fe", "k_m": "k_m_k20fe"}
WISE_EXTENDED_COLUMNS = {"w1mpro": "w1gmag", "w2mpro": "w2gmag"}

# A source this close in both catalogues is the same object.
CROSS_MATCH_ARCSEC = 5.0

# Catalogue search radius in arcsec. By 30 arcsec the PSF crosstalk weight
# is already below 0.01, so nothing further contributes.
SEARCH_RADIUS = 30.0


def _query_irsa(coord: SkyCoord, catalogue: str, columns: list[str]) -> pd.DataFrame:
    """
    Fetch positions and magnitudes from one IRSA catalogue.

    :param coord: Position of the target
    :param catalogue: IRSA table name
    :param columns: Magnitude columns to retrieve
    :return: Positions and magnitudes, empty if the query fails
    """
    ra, dec = coord.ra.deg, coord.dec.deg
    selected = ", ".join(["ra", "dec"] + columns)

    query = f"""
    SELECT {selected}
    FROM {catalogue}
    WHERE CONTAINS(POINT('ICRS', ra, dec),
                   CIRCLE('ICRS', {ra}, {dec}, {SEARCH_RADIUS / 3600})) = 1
    """
    try:
        return pyvo.dal.TAPService(TAP_URL).search(query).to_table().to_pandas()
    except Exception as exc:  # noqa: BLE001 - any archive failure is non-fatal
        logger.warning(f"Could not query {catalogue}: {exc}")
        return pd.DataFrame()


def _build_2mass(coord: SkyCoord) -> pd.DataFrame:
    """
    Combine the 2MASS point and extended source catalogues.

    Where a galaxy appears in both, the extended magnitudes replace the
    point-source ones, which measure only its core.

    :param coord: Position of the target
    :return: Positions and J, H, Ks magnitudes
    """
    point = _query_irsa(coord, "fp_psc", ["j_m", "h_m", "k_m"])
    extended = _query_irsa(coord, "fp_xsc", list(XSC_COLUMNS.values()))

    if len(extended) == 0:
        return point

    extended = extended.rename(columns={v: k for k, v in XSC_COLUMNS.items()})

    if len(point) > 0:
        # Drop point-source entries that are really the extended sources.
        point_positions = SkyCoord(
            point["ra"].to_numpy(), point["dec"].to_numpy(), unit="deg"
        )
        duplicate = np.zeros(len(point), dtype=bool)
        for _, row in extended.iterrows():
            centre = SkyCoord(row["ra"], row["dec"], unit="deg")
            duplicate |= centre.separation(point_positions).arcsec < CROSS_MATCH_ARCSEC
        if duplicate.any():
            logger.info(
                f"{duplicate.sum()} 2MASS point source(s) replaced by their "
                f"extended-source photometry"
            )
        point = point[~duplicate]

    return pd.concat([point, extended], ignore_index=True)


def _build_wise(coord: SkyCoord) -> pd.DataFrame:
    """
    Fetch AllWISE photometry, preferring aperture magnitudes where extended.

    :param coord: Position of the target
    :return: Positions and W1, W2 magnitudes
    """
    columns = list(WISE_EXTENDED_COLUMNS) + list(WISE_EXTENDED_COLUMNS.values())
    table = _query_irsa(coord, "allwise_p3as_psd", columns + ["ext_flg"])
    if len(table) == 0:
        return table

    extended = pd.to_numeric(table["ext_flg"], errors="coerce").fillna(0) > 0
    for profile, aperture in WISE_EXTENDED_COLUMNS.items():
        aperture_mag = pd.to_numeric(table[aperture], errors="coerce")
        use = extended & aperture_mag.notna()
        table.loc[use, profile] = aperture_mag[use]

    if extended.any():
        logger.info(
            f"{int(extended.sum())} WISE source(s) using elliptical aperture magnitudes"
        )

    return table


def archival_photometry(coord: SkyCoord) -> pd.DataFrame:
    """
    Estimate the pre-SPHEREx flux a forced fit at a position would have measured.

    Every catalogued source within the search radius contributes its flux times
    the fraction of its light a PSF centred on the target would collect, so a
    bright star a few arcsec away counts for most of its flux and one 15 arcsec
    away counts for almost none.

    :param coord: Position of the target
    :return: One row per band, with the summed flux in uJy
    """
    tables = {"2mass": _build_2mass(coord), "wise": _build_wise(coord)}

    rows = []
    for band, catalogue, column, wavelength, zero_point in BANDS:
        if not WAVELENGTH_MIN <= wavelength <= WAVELENGTH_MAX:
            continue

        table = tables.get(catalogue, pd.DataFrame())
        if (len(table) == 0) or (column not in table.columns):
            continue

        magnitude = pd.to_numeric(table[column], errors="coerce").to_numpy()
        usable = np.isfinite(magnitude)
        if usable.sum() == 0:
            continue

        positions = SkyCoord(
            table["ra"].to_numpy()[usable], table["dec"].to_numpy()[usable], unit="deg"
        )
        separation = coord.separation(positions).arcsec
        flux = zero_point * 10 ** (-magnitude[usable] / 2.5) * 1e6  # uJy
        weight = psf_weight(separation, wavelength)

        rows.append(
            {
                "band": band,
                "wavelength": wavelength,
                "flux": float(np.sum(flux * weight)),
                "flux_brightest": float(np.max(flux * weight)),
                "n_sources": int(usable.sum()),
                "nearest_arcsec": float(separation.min()),
                "catalogue": catalogue,
            }
        )

    estimate = pd.DataFrame(rows).sort_values("wavelength").reset_index(drop=True)

    if len(estimate) == 0:
        logger.warning("No archival photometry found near the target")
    else:
        lowest = estimate["flux"].min() / 1e3
        highest = estimate["flux"].max() / 1e3
        logger.info(
            f"Archival estimate in {len(estimate)} bands, "
            f"{lowest:.2f} to {highest:.2f} mJy"
        )

    return estimate
