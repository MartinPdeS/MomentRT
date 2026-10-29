"""Capture geometry figures produced by the Sphinx-Gallery documentation build.

Run after ``make docs`` with Playwright and Chromium installed. The captured
images are committed so documentation builds do not need a browser.
"""

import argparse
import re
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent.parent

# Each example ends with its physical geometry, or observation window.
EXAMPLES = tuple(sorted((ROOT / "docs/examples").glob("*/*.py")))


def main() -> None:
    """Render each generated geometry in a browser and save a gallery preview."""
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

        for example in EXAMPLES:
            name = example.stem

            rst = ROOT / "docs/source/gallery" / example.parent.name / f"{name}.rst"

            figures = re.findall(
                rf"_static/gallery/{re.escape(name)}-(\d+)\.html",
                rst.read_text(),
            )

            if not figures:
                raise ValueError(f"No Plotly figures in {rst}; run make docs first.")

            html = ROOT / f"docs/_build/html/_static/gallery/{name}-{figures[-1]}.html"

            if not html.exists():
                raise FileNotFoundError(f"Missing {html}; run make docs first.")

            page.goto(html.as_uri(), wait_until="load")

            page.wait_for_selector(".main-svg", timeout=30000)

            page.wait_for_function(
                "document.querySelector('.js-plotly-plot')?._fullLayout?.scene?._scene?.glplot",
                timeout=30000,
            )

            destination = output / "gallery_systems" / f"{name}.png"

            destination.parent.mkdir(parents=True, exist_ok=True)

            page.screenshot(path=str(destination))

            print(f"Captured {name}: {destination}")

        browser.close()


if __name__ == "__main__":
    main()
