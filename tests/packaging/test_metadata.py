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


def test_plotly_is_a_core_dependency():
    root = Path.cwd()

    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]

    recipe = (root / "conda.recipe" / "meta.yaml").read_text()

    assert "plotly>=6,<7" in project["dependencies"]

    run_requirements = recipe.split("  run:\n", 1)[1].split("\ntest:", 1)[0]

    assert "    - plotly >=6,<7" in run_requirements


def test_typedunit_and_pint_requirements_match_conda_recipe():
    """Both distributions use the supported TypedUnit and Pint versions."""
    root = Path.cwd()
    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]
    recipe = (root / "conda.recipe" / "meta.yaml").read_text()
    run_requirements = recipe.split("  run:\n", 1)[1].split("\ntest:", 1)[0]

    assert "TypedUnit>=0.3.0,<0.4" in project["dependencies"]
    assert "pint>=0.24,<0.26" in project["dependencies"]
    assert "    - typedunit >=0.3.0,<0.4" in run_requirements
    assert "    - pint >=0.24,<0.26" in run_requirements
