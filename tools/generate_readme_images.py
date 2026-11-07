"""Capture README figures from the executed skin contact gallery example.

Run after ``make docs`` with Playwright and Chromium installed. The generated
images are committed so documentation and package builds do not need a browser.
"""

from __future__ import annotations
import argparse
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    """Capture the system, mesh, and fluence views from the built example."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--browser",
        help="Path to Chromium or Chrome; omit to use Playwright's Chromium.",
    )
    arguments = parser.parse_args()
    gallery = ROOT / "docs/_build/html/_static/gallery"
    destination = ROOT / "docs/images"
    figures = (
        ("skin_contact-1.html", "readme-skin-system.png"),
        ("skin_contact-2.html", "readme-skin-mesh.png"),
        ("skin_contact-3.html", "readme-skin-fluence.png"),
    )
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=arguments.browser,
            headless=True,
        )
        page = browser.new_page(viewport={"width": 1000, "height": 750})
        for source, output in figures:
            html = gallery / source
            if not html.is_file():
                raise FileNotFoundError(f"Missing {html}; run make docs first.")
            page.goto(html.as_uri(), wait_until="load")
            page.wait_for_selector(".main-svg", timeout=30000)
            page.wait_for_function(
                "document.querySelector('.js-plotly-plot')?._fullLayout?.scene?._scene?.glplot",
                timeout=30000,
            )
            path = destination / output
            page.screenshot(path=str(path))
            print(f"Captured {path}")
        browser.close()


if __name__ == "__main__":
    main()
