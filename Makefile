PYTHON ?= python3
PY_ENV ?= .venv
PYTHON_BIN := $(PY_ENV)/bin/python
JULIA ?= julia

.PHONY: help setup-python setup-julia analysis bootstrap test test-python test-julia validate-python handbook

help:
	@echo "make setup-python  Install Python dependencies and tests"
	@echo "make analysis      Verify public data; reproduce coefficients and figure"
	@echo "make bootstrap     Independently rerun the fixed resampling"
	@echo "make test-python   Run the Python tests"
	@echo "make test-julia    Run Julia tests and small command-line workflows"
	@echo "make test          Run both language test suites"
	@echo "make handbook      Build the Typst research handbook"

setup-python:
	"$(PYTHON)" -m venv "$(PY_ENV)"
	"$(PYTHON_BIN)" -m pip install --upgrade pip
	"$(PYTHON_BIN)" -m pip install -e './square_glauber_python[test,accel]'

# Julia dependencies are standard libraries; this checks the package loads.
setup-julia:
	"$(JULIA)" --startup-file=no --project=aztec -e 'using AztecDiamond; println("AztecDiamond ready")'

analysis:
	OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=tmp/matplotlib "$(PYTHON_BIN)" analysis/run.py

bootstrap:
	OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=tmp/matplotlib "$(PYTHON_BIN)" analysis/run.py --bootstrap 10000 --output output/recomputed

test-python:
	OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=tmp/matplotlib PYTHONPATH=square_glauber_python/src \
		"$(PYTHON_BIN)" -m pytest square_glauber_python/tests

validate-python:
	MPLCONFIGDIR=tmp/matplotlib PYTHONPATH=square_glauber_python/src "$(PYTHON_BIN)" -m square_glauber.cli validate-exact

test-julia:
	"$(JULIA)" --startup-file=no --project=aztec aztec/test/runtests.jl
	PATH="$(dir $(JULIA)):$$PATH" sh aztec/test/smoke_workflows.sh

test: test-python test-julia

handbook:
	mkdir -p output/pdf
	typst compile --root . docs/research/handbook.typ output/pdf/research-handbook.pdf
