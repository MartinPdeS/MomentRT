"""Verify and install a pure Python wheel in an isolated environment."""

import argparse
import subprocess
import sys
import tempfile
import tomllib
import venv
from email import message_from_bytes
from pathlib import Path
from zipfile import ZipFile

REQUIRED_MODULES = (
    "MomentRT/__init__.py",
    "MomentRT/py.typed",
    "MomentRT/_version.py",
    "MomentRT/boundary.py",
    "MomentRT/convergence.py",
    "MomentRT/detector.py",
    "MomentRT/diffusion.py",
    "MomentRT/fem.py",
    "MomentRT/fields.py",
    "MomentRT/geometry.py",
    "MomentRT/material.py",
    "MomentRT/medium.py",
    "MomentRT/mesh.py",
    "MomentRT/p3.py",
    "MomentRT/plotting.py",
    "MomentRT/results.py",
    "MomentRT/sections.py",
    "MomentRT/simulation.py",
    "MomentRT/source.py",
    "MomentRT/sp3.py",
    "MomentRT/units.py",
    "MomentRT/validation.py",
)


def check_archive(*, wheel):
    if not wheel.name.endswith("-py3-none-any.whl"):
        raise ValueError("Expected a platform-independent Python wheel.")

    root = Path(__file__).resolve().parent.parent

    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]

    with ZipFile(wheel) as archive:
        names = set(archive.namelist())

        missing = set(REQUIRED_MODULES) - names

        if missing:
            raise ValueError(f"Missing package modules: {sorted(missing)}")

        metadata_files = [
            name for name in names if name.endswith(".dist-info/METADATA")
        ]

        if len(metadata_files) != 1:
            raise ValueError("Expected exactly one wheel metadata file.")

        metadata = message_from_bytes(archive.read(metadata_files[0]))

        if (
            metadata["Name"].lower() != project["name"].lower()
            or metadata["Version"] != project["version"]
        ):
            raise ValueError("Wheel name or version disagrees with pyproject.toml.")

        if metadata["License-Expression"] != "MIT" or not any(
            name.endswith(".dist-info/licenses/LICENSE") for name in names
        ):
            raise ValueError("Wheel must include MIT license metadata and text.")


def check_install(*, wheel):
    with tempfile.TemporaryDirectory(prefix="momentrt-wheel-") as directory:
        temporary = Path(directory)

        environment = temporary / "venv"

        venv.EnvBuilder(with_pip=True).create(environment)

        python = environment / (
            "Scripts/python.exe" if sys.platform == "win32" else "bin/python"
        )

        subprocess.run(
            [str(python), "-m", "pip", "install", str(wheel)],
            cwd=temporary,
            check=True,
        )

        command = (
            "from pathlib import Path; import sys; import MomentRT; "
            "location = Path(MomentRT.__file__).resolve(); "
            "assert location.is_relative_to(Path(sys.prefix).resolve()), location; "
            "assert all(hasattr(MomentRT, name) for name in "
            "('SurfaceDetector', 'SpectralMaterial', 'SpectralSimulationResult', "
            "'DetectorConvergenceResult')); "
            "print(f'Imported MomentRT {MomentRT.__version__} from {location}')"
        )

        subprocess.run(
            [str(python), "-c", command],
            cwd=temporary,
            check=True,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("wheel", type=Path)

    arguments = parser.parse_args()

    wheel = arguments.wheel.resolve()

    check_archive(wheel=wheel)

    check_install(wheel=wheel)

    print(f"Verified platform-independent wheel: {wheel.name}")


if __name__ == "__main__":
    main()
