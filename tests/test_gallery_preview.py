"""Gallery preview selection remains stable when later figures are added."""

import importlib.util
from pathlib import Path

import pytest

plotly = pytest.importorskip("plotly.graph_objects")


def test_explicit_preview_is_not_replaced_by_a_later_figure(tmp_path):
    path = Path(__file__).resolve().parents[1] / "docs/source/gallery_config.py"

    spec = importlib.util.spec_from_file_location("gallery_config", path)

    module = importlib.util.module_from_spec(spec)

    spec.loader.exec_module(module)

    example = tmp_path / "gallery" / "example.py"

    example.parent.mkdir()

    selected = plotly.Figure()

    variables = {
        "target_file": str(example),
        "example_globals": {"gallery_system_figure": selected},
    }

    config = {"src_dir": str(tmp_path)}

    module.plotly_scraper(None, variables, config)

    manifest = tmp_path / "_static/gallery/example-preview.txt"

    assert manifest.read_text() == "example-1.html"

    variables["example_globals"]["later_figure"] = plotly.Figure()

    module.plotly_scraper(None, variables, config)

    assert manifest.read_text() == "example-1.html"
