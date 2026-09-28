"""
End-to-end SPHEREx forced photometry for a transient.

Can be used as a library::

    from obsutils.spherex import run_photometry
    result = run_photometry("AT2025abcr", reference_mjd=60961.28)

or from the command line::

    spherex-phot AT2025abcr --reference-mjd 60961.28 --deblend --plot

Output always goes to <OUTPUT_DIR>/<name>/spherex, with OUTPUT_DIR read from
.env.
"""

import argparse
import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from dotenv import load_dotenv

from obsutils.spherex.archival import archival_photometry
from obsutils.spherex.constants import default_wavelength_edges
from obsutils.spherex.photometry import clean, measure_images
from obsutils.spherex.plot import (
    plot_coverage,
    plot_difference,
    plot_lightcurve,
    plot_sed,
)
from obsutils.spherex.query import (
    find_companions,
    query_images,
    resolve_position,
    summarise_coverage,
)
from obsutils.spherex.stack import difference_stacks, stack_by_block

load_dotenv()

logger = logging.getLogger(__name__)


@dataclass
class SpherexResult:
    """
    Everything produced for one target.

    :param name: Name of the target
    :param coord: Position of the target
    :param reference_mjd: Epoch used to split the stacks, if any
    :param images: One row per spectral image covering the position
    :param photometry: One row per measured image
    :param archival: Pre-SPHEREx flux estimated from archival catalogues
    :param stacks: Wavelength-stacked spectrum per block, baseline first
    :param difference: Each block minus the baseline, if there is a baseline
    :param error_scale: Factor the formal errors were rescaled by
    """

    name: str
    coord: SkyCoord
    reference_mjd: float | None
    images: pd.DataFrame
    photometry: pd.DataFrame
    stacks: dict[str, pd.DataFrame] = field(default_factory=dict)
    archival: pd.DataFrame | None = None
    difference: pd.DataFrame | None = None
    output_dir: Path | None = None
    error_scale: float = 1.0

    @property
    def coverage(self) -> pd.DataFrame:
        """
        Images per detector and phase.

        :return: Summary table
        """
        return summarise_coverage(self.images)


def default_name(coord: SkyCoord) -> str:
    """
    Build an IAU-style name for a position, for targets given only as coordinates.

    :param coord: Position of the target
    :return: Name of the form J014655.40-152215.7
    """
    ra = coord.ra.to_string(unit=u.hourangle, sep="", precision=2, pad=True)
    dec = coord.dec.to_string(sep="", precision=1, alwayssign=True, pad=True)
    return f"J{ra}{dec}"


def safe_name(name: str) -> str:
    """
    Sanitise a target name for use in a file path.

    :param name: Name of the target
    :return: Name with characters awkward in a path replaced
    """
    return re.sub(r"[^A-Za-z0-9.+-]+", "_", name).strip("_")


def get_output_dir(name: str) -> Path:
    """
    Build the output directory for a target, creating it if needed.

    Always <OUTPUT_DIR>/<name>/spherex, reading OUTPUT_DIR from .env as the
    other obsutils notebooks do.

    :param name: Name of the target
    :return: The directory
    """
    base = Path(os.getenv("OUTPUT_DIR", str(Path.home()) + "/obsutils"))
    output_dir = (base / safe_name(name) / "spherex").expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir


def run_photometry(
    name: str | None = None,
    ra: float | None = None,
    dec: float | None = None,
    reference_mjd: float | None = None,
    deblend: bool = False,
    wavelength_edges: np.ndarray | None = None,
) -> SpherexResult:
    """
    Measure a transient in every SPHEREx spectral image covering it.

    :param name: Name of the target, resolved through Sesame
    :param ra: Right ascension in degrees, overriding the name
    :param dec: Declination in degrees, overriding the name
    :param reference_mjd: Epoch splitting the stacks, e.g. discovery
    :param deblend: Fit catalogued neighbours alongside the target and tilt the
        background, for a target blended with a host or a crowded field
    :param wavelength_edges: Wavelength bin edges in micron
    :return: Images, photometry, stacks and difference spectrum
    """
    coord = resolve_position(name=name, ra=ra, dec=dec)
    label = name if name is not None else default_name(coord)

    if wavelength_edges is None:
        wavelength_edges = default_wavelength_edges()

    images = query_images(coord, reference_mjd=reference_mjd)

    if len(images) == 0:
        return SpherexResult(
            name=label,
            coord=coord,
            reference_mjd=reference_mjd,
            images=images,
            photometry=pd.DataFrame(),
        )

    companions = find_companions(coord) if deblend else None
    photometry = measure_images(images, coord, companions=companions)

    measured = clean(photometry) if len(photometry) > 0 else photometry

    stacks, difference, error_scale = {}, None, 1.0
    if len(measured) > 0:
        stacks = stack_by_block(
            measured, edges=wavelength_edges, reference_mjd=reference_mjd
        )
        if len(stacks) > 0:
            first = next(iter(stacks.values()))
            if len(first) > 0:
                error_scale = float(
                    (first["flux_err"] / first["flux_err_formal"]).median()
                )
        difference = difference_stacks(stacks)

    # The overlay shows what was at this position before SPHEREx, so it includes
    # the companions even when they were fitted out of the target's amplitude.
    archival_estimate = archival_photometry(coord)

    result = SpherexResult(
        name=label,
        coord=coord,
        reference_mjd=reference_mjd,
        images=images,
        photometry=photometry,
        stacks=stacks,
        archival=archival_estimate,
        difference=difference,
        error_scale=error_scale,
    )

    result.output_dir = save_result(result)

    return result


def save_result(result: SpherexResult) -> Path:
    """
    Write the photometry, stacks and difference spectrum to CSV.

    :param result: Result to write
    :return: The directory written to
    """
    directory = get_output_dir(result.name)
    stem = safe_name(result.name)

    result.photometry.to_csv(directory / f"{stem}_spherex_photometry.csv", index=False)

    if len(result.stacks) > 0:
        stacked = pd.concat(
            [stack.assign(block=label) for label, stack in result.stacks.items()],
            ignore_index=True,
        )
        stacked.to_csv(directory / f"{stem}_spherex_stacks.csv", index=False)

    if result.archival is not None:
        result.archival.to_csv(directory / f"{stem}_spherex_archival.csv", index=False)

    if result.difference is not None:
        result.difference.to_csv(
            directory / f"{stem}_spherex_difference.csv", index=False
        )

    logger.info(f"Wrote SPHEREx products to {directory}")
    return directory


def save_plots(result: SpherexResult) -> Path:
    """
    Write the coverage, light curve, SED and difference plots to PNG.

    :param result: Result to plot
    :return: The directory written to
    """
    directory = get_output_dir(result.name)

    figures = {
        "coverage": plot_coverage(
            result.images, reference_mjd=result.reference_mjd, title=result.name
        ),
        "lightcurve": plot_lightcurve(
            result.photometry,
            reference_mjd=result.reference_mjd,
            title=result.name,
        ),
        "sed": plot_sed(
            result.stacks,
            title=result.name,
            difference=result.difference,
            archival=result.archival,
        ),
    }
    if result.difference is not None:
        figures["difference"] = plot_difference(result.difference, title=result.name)

    stem = safe_name(result.name)
    for kind, figure in figures.items():
        figure.savefig(directory / f"{stem}_spherex_{kind}.png", dpi=150)
        plt.close(figure)

    logger.info(f"Wrote SPHEREx plots to {directory}")
    return directory


def block_summary(difference: pd.DataFrame) -> pd.DataFrame:
    """
    One line per block: its mean offset from the baseline and its scatter.

    :param difference: Difference table, from difference_stacks
    :return: Summary table
    """
    rows = []
    for label, group in difference.groupby("block", sort=False):
        weights = 1.0 / group["diff_err"] ** 2
        offset = float(np.sum(group["diff"] * weights) / np.sum(weights))
        error = float(1.0 / np.sqrt(np.sum(weights)))
        rows.append(
            {
                "block": label,
                "phase_days": round(float(group["phase_days"].iloc[0]), 1),
                "bins": len(group),
                "offset_uJy": round(offset, 1),
                "offset_err": round(error, 1),
                "sigma": round(offset / error, 1),
                "rms_sigma": round(float(np.sqrt(np.mean(group["sigma"] ** 2))), 2),
            }
        )
    return pd.DataFrame(rows)


def main(args: list[str] | None = None):
    """
    Command line entry point.

    :param args: Argument list, defaulting to sys.argv
    :return: None
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "name", nargs="?", help="transient name, resolved through Sesame"
    )
    parser.add_argument("--ra", type=float, help="right ascension in degrees")
    parser.add_argument("--dec", type=float, help="declination in degrees")
    parser.add_argument(
        "--reference-mjd", type=float, help="epoch splitting the stacks, e.g. discovery"
    )
    parser.add_argument(
        "--deblend",
        action="store_true",
        help="fit catalogued neighbours alongside the target, for a blended "
        "or crowded position",
    )
    parser.add_argument("--plot", action="store_true", help="also write plots")
    parsed = parser.parse_args(args)

    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s %(name)s: %(message)s"
    )

    result = run_photometry(
        name=parsed.name,
        ra=parsed.ra,
        dec=parsed.dec,
        reference_mjd=parsed.reference_mjd,
        deblend=parsed.deblend,
    )

    if len(result.photometry) == 0:
        logger.warning(f"No SPHEREx photometry for {result.name}")
        return

    print(result.coverage)
    if result.difference is not None:
        print(block_summary(result.difference).to_string(index=False))

    if parsed.plot:
        save_plots(result)


if __name__ == "__main__":
    main()
