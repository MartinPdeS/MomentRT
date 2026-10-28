.DEFAULT_GOAL := help

PYTHON ?= .venv/bin/python
BUILD_DIR ?= build
TAG_VERSION ?= $(or $(VERSION),$(filter v%,$(MAKECMDGOALS)))
RELEASE_KIND := $(filter major minor patch,$(MAKECMDGOALS))

.PHONY: help configure build install uninstall quick rebuild editable clean quality typecheck stubcheck test check package-check release-check tag release major minor patch

help:
	@echo "MomentRT development commands"
	@echo ""
	@echo "  make configure             Check whether native configuration is needed"
	@echo "  make build                 Build source and wheel distributions"
	@echo "  make install               Install the package"
	@echo "  make uninstall             Uninstall the package"
	@echo "  make quick                 Install the package"
	@echo "  make rebuild               Clean, build, and install the package"
	@echo "  make editable              Install an editable package"
	@echo "  make clean                 Remove build artifacts"
	@echo "  make quality               Run static checks"
	@echo "  make typecheck             Report Python API typing support"
	@echo "  make stubcheck             Report native binding stub support"
	@echo "  make test                  Run the test suite"
	@echo "  make check                 Run quality and tests"
	@echo "  make docs                  Build the executable example gallery"
	@echo "  make package-check         Build source and wheel distributions"
	@echo "  make release-check         Check version metadata consistency"
	@echo "  make tag VERSION=vX.Y.Z    Create a release commit and annotated tag"
	@echo "  make release patch         Create and push the next patch release"
	@echo "  make release minor         Create and push the next minor release"
	@echo "  make release major         Create and push the next major release"

configure:
	@echo "MomentRT has no native build to configure."

build package-check:
	$(PYTHON) -m build

install:
	$(PYTHON) -m pip install .

uninstall:
	$(PYTHON) -m pip uninstall -y MomentRT

quick: install

rebuild:
	$(MAKE) clean
	$(MAKE) build
	$(MAKE) install

editable:
	$(PYTHON) -m pip install -e .

clean:
	rm -rf $(BUILD_DIR) dist

quality:
	$(PYTHON) -m ruff check MomentRT tests docs/examples tools
	$(PYTHON) -m ruff format --check MomentRT tests docs/examples tools

typecheck:
	@echo "MomentRT does not currently provide a typed Python API."

stubcheck:
	@echo "MomentRT has no native binding stubs to check."

test:
	$(PYTHON) -m pytest

check: quality test

release-check:
	$(PYTHON) -m martinpdes_maintenance.cli check $(if $(VERSION),--version $(VERSION),)

ifneq ($(filter tag,$(MAKECMDGOALS)),)
ifneq ($(strip $(TAG_VERSION)),)
.PHONY: $(TAG_VERSION)
$(TAG_VERSION):
	@:
endif
endif

tag:
	$(PYTHON) -m martinpdes_maintenance.cli tag "$(TAG_VERSION)"

release:
	@test "$(words $(RELEASE_KIND))" -eq 1 || { echo "usage: make release [patch|minor|major]" >&2; exit 2; }
	@set -eu; release_tag="$$($(PYTHON) -m martinpdes_maintenance.cli next $(RELEASE_KIND))"; \
	$(PYTHON) -m martinpdes_maintenance.cli tag "$$release_tag"; \
	git push origin HEAD "refs/tags/$$release_tag"

major minor patch:
	@:

.PHONY: docs
docs:
	$(PYTHON) -m sphinx -W --keep-going -b html docs/source docs/_build/html
