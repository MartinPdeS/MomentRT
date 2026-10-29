"""Capture volume figures produced by the Sphinx-Gallery documentation build.

Run after ``make docs`` with Playwright and Chromium installed. The captured
images are committed so documentation builds do not need a browser.
"""

import argparse
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent.parent

# Each entry identifies an existing Plotly volume embedded by its example.
VOLUMES = {
    "absorption_sweep": (1, "gallery_volumes/absorption_sweep.png"),
    "infinite_diffusion": (2, "gallery_volumes/infinite_diffusion.png"),
    "compare_p3": (1, "p3-fluence.png"),
    "compare_sp3": (1, "sp3-fluence.png"),
    "detector_convergence": (2, "gallery_volumes/detector_convergence.png"),
    "mesh_refinement": (1, "gallery_volumes/mesh_refinement.png"),
    "reflectance_sweep": (1, "gallery_volumes/reflectance_sweep.png"),
    "seven_fiber_skin": (2, "seven-fiber-skin.png"),
    "sp3_layered_anisotropy": (2, "sp3_layered_anisotropy.png"),
    "sp3_spectral_anisotropy": (2, "sp3_spectral_anisotropy.png"),
    "spectral_detectors": (3, "gallery_volumes/spectral_detectors.png"),
    "structured_sample": (2, "gallery_volumes/structured_sample.png"),
    "surface_detector": (2, "gallery_volumes/surface_detector.png"),
    "surface_diffusion": (2, "gallery_volumes/surface_diffusion.png"),
    "uniform_surface_source": (1, "gallery_volumes/uniform_surface_source.png"),
}


def main() -> None:
    """Render each generated volume in a browser and save a gallery preview."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--browser",
        help="Path to Chromium or Chrome; omit to use Playwright's Chromium.",
    )

    arguments = parser.parse_args()

    output = ROOT / "docs/source/_static/images"

    output.mkdir(parents=True, exist_ok=True)

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=arguments.browser,
            headless=True,
        )

        page = browser.new_page(viewport={"width": 1000, "height": 750})

        for name, (number, filename) in VOLUMES.items():
            html = ROOT / f"docs/_build/html/_static/gallery/{name}-{number}.html"

            if not html.exists():
                raise FileNotFoundError(f"Missing {html}; run make docs first.")

            page.goto(html.as_uri(), wait_until="load")

            page.wait_for_selector(".main-svg", timeout=30000)

            page.wait_for_function(
                "document.querySelector('.js-plotly-plot')?._fullLayout?.scene?._scene?.glplot",
                timeout=30000,
            )

            destination = output / filename

            destination.parent.mkdir(parents=True, exist_ok=True)

            page.screenshot(path=str(destination))

            print(f"Captured {name}: {destination}")

        browser.close()


if __name__ == "__main__":
    main()
