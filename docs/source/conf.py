"""Sphinx configuration with executable Plotly gallery examples."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from MomentRT import __version__

project = "MomentRT"

release = __version__

extensions = ["sphinx.ext.mathjax", "sphinx_gallery.gen_gallery"]

html_theme = "alabaster"

html_static_path = ["_static"]

exclude_patterns = ["_build"]

sphinx_gallery_conf = {
    "examples_dirs": "../examples",
    "gallery_dirs": "auto_examples",
    "filename_pattern": r"\.py$",
    "image_scrapers": ("gallery_config.plotly_scraper",),
    "capture_repr": (),
    "abort_on_example_error": True,
    "download_all_examples": False,
    "backreferences_dir": None,
}
