"""RePrint: DNA Repair FootPrint computed from COSMIC mutational signatures."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("reprint-tool")
except PackageNotFoundError:  # running from a source checkout without installation
    __version__ = "0.0.0+unknown"

from reprint_tool.analyze import (
    calculate_cosine,
    calculate_rmse,
    create_heatmap_with_custom_sim,
)
from reprint_tool.core import reprint
from reprint_tool.plot_static import generate_all_reprints, generate_single_reprint

__all__ = [
    "__version__",
    "calculate_cosine",
    "calculate_rmse",
    "create_heatmap_with_custom_sim",
    "generate_all_reprints",
    "generate_single_reprint",
    "reprint",
]
