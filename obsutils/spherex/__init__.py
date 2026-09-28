"""
Query and measure SPHEREx spectral images for a transient.

The usual entry point is run_photometry, which finds every spectral image
covering a position, measures the target in each by forced PSF fitting, and
stacks the result into a spectrum either side of a reference epoch.
"""

from obsutils.spherex.archival import archival_photometry
from obsutils.spherex.constants import default_wavelength_edges
from obsutils.spherex.photometry import clean, measure_cutout, measure_images
from obsutils.spherex.pipeline import (
    SpherexResult,
    run_photometry,
    save_plots,
    save_result,
)
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
from obsutils.spherex.stack import (
    assign_blocks,
    difference_stacks,
    error_inflation_factor,
    stack_by_block,
    stack_in_wavelength,
)

__all__ = [
    "SpherexResult",
    "archival_photometry",
    "assign_blocks",
    "clean",
    "default_wavelength_edges",
    "difference_stacks",
    "error_inflation_factor",
    "find_companions",
    "measure_cutout",
    "measure_images",
    "plot_coverage",
    "plot_difference",
    "plot_lightcurve",
    "plot_sed",
    "query_images",
    "resolve_position",
    "run_photometry",
    "save_plots",
    "save_result",
    "stack_by_block",
    "stack_in_wavelength",
    "summarise_coverage",
]
