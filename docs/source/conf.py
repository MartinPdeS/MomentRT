"""Sphinx configuration for MomentRT documentation."""

import os
import sys
from pathlib import Path

from martinpdes_docs import configure

package_name = "MomentRT"
docs_source = Path(__file__).resolve().parent
sys.path.insert(0, str(docs_source))
os.environ["MOMENTRT_SPHINX_GALLERY"] = "1"

globals().update(
    configure(
        project=package_name,
        package=package_name,
        version="auto",
        author="Martin Poinsinet de Sivry-Houle",
        repository=f"MartinPdeS/{package_name}",
        logo="../images/logo-light.svg",
        favicon="../images/favicon.ico",
        html_theme_options={
            "logo": {
                "image_light": "../images/logo-light.svg",
                "image_dark": "../images/logo-dark.svg",
            },
        },
        prefer_installed_package=True,
        gallery_sections=["fundamentals", "finite_samples"],
        gallery_default_thumb=str(
            docs_source / "_static" / "images" / "sample-geometry.png"
        ),
        gallery_options={
            "image_scrapers": ("matplotlib", "gallery_config.plotly_scraper"),
            "abort_on_example_error": True,
            "capture_repr": (),
        },
        logo_text="",
        logo_alt_text="MomentRT deterministic light transport",
        html_css_files=["layout.css"],
        autodoc_default_options={
            "members": True,
            "member-order": "bysource",
        },
    )
)

# Class pages own their member documentation; avoid numpydoc-generated links
# to separate method pages that are not built by this reference.
numpydoc_show_class_members = False
