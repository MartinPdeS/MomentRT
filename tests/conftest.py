"""Keep example-driven tests noninteractive on all platforms."""

import os

os.environ["MOMENTRT_SPHINX_GALLERY"] = "1"
os.environ["MPLBACKEND"] = "Agg"
