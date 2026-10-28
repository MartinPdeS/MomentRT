"""Validate release metadata against the published package configuration."""

import json
import re
import tomllib
from pathlib import Path

import MomentRT


def test_release_metadata_and_license_are_aligned():
    root = Path.cwd()

    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]

    citation = (root / "CITATION.cff").read_text()

    recipe = (root / "conda.recipe" / "meta.yaml").read_text()

    zenodo = json.loads((root / ".zenodo.json").read_text())

    citation_version = re.search(r'^version:\s*"([^"]+)"\s*$', citation, re.MULTILINE)

    recipe_version = re.search(r'^  version:\s*"([^"]+)"\s*$', recipe, re.MULTILINE)

    assert citation_version is not None

    assert recipe_version is not None

    assert project["version"] == MomentRT.__version__

    assert citation_version.group(1) == project["version"]

    assert recipe_version.group(1) == project["version"]

    assert zenodo["version"] == project["version"]

    assert project["license"] == "MIT"

    assert re.search(r'^license:\s*"MIT"\s*$', citation, re.MULTILINE)

    assert re.search(r"^  license:\s*MIT\s*$", recipe, re.MULTILINE)

    assert zenodo["license"] == "mit"

    assert zenodo["creators"][0]["name"] == "Poinsinet de Sivry-Houle, Martin"

    assert "MIT License" in (root / "LICENSE").read_text()
