PYTHON ?= python3.13
PY_ENV ?= square_glauber_python/.venv
PYTHON_BIN := $(PY_ENV)/bin/python
JULIA ?= julia

.PHONY: help setup-python setup-julia demo-aztec demo-glauber figures test-python validate-python test-julia test

help:
	@echo "make setup-julia                     Instantiate the Julia project"
	@echo "make setup-python PYTHON=python3.13  Create the Python environment"
	@echo "make demo-aztec                      Generate a small Aztec tiling"
	@echo "make demo-glauber                    Run a small Glauber diagnostic"
	@echo "make figures                         Rebuild README illustrations"
	@echo "make test-python                     Run square-grid Python tests"
	@echo "make validate-python                 Run exact tiny-grid validation"
	@echo "make test-julia                      Run the Julia/Aztec test suite"
	@echo "make test                            Run both project test suites"

setup-python:
	$(PYTHON) -m venv $(PY_ENV)
	$(PYTHON_BIN) -m pip install --upgrade pip
	$(PYTHON_BIN) -m pip install -e './square_glauber_python[test,accel]'

setup-julia:
	$(JULIA) --project=aztec -e 'using Pkg; Pkg.instantiate()'

demo-aztec:
	$(JULIA) --project=aztec aztec/scripts/run_gamma_disordered.jl \
		--order 24 --seed 20260920 --output-dir aztec/output/demo

demo-glauber:
	MPLCONFIGDIR=$(PY_ENV)/matplotlib-cache XDG_CACHE_HOME=$(PY_ENV)/cache \
		PYTHONPATH=square_glauber_python/src $(PYTHON_BIN) -m square_glauber.cli mixing \
		--L 8 --model uniform --sweeps 500 --record-every 5 --seed 20260920 \
		--output square_glauber_python/outputs/demo

figures:
	MPLCONFIGDIR=$(PY_ENV)/matplotlib-cache XDG_CACHE_HOME=$(PY_ENV)/cache \
		PYTHONPATH=square_glauber_python/src $(PYTHON_BIN) docs/scripts/render_figures.py

test-python:
	MPLCONFIGDIR=$(PY_ENV)/matplotlib-cache XDG_CACHE_HOME=$(PY_ENV)/cache \
		PYTHONPATH=square_glauber_python/src $(PYTHON_BIN) -m pytest square_glauber_python/tests

validate-python:
	MPLCONFIGDIR=$(PY_ENV)/matplotlib-cache XDG_CACHE_HOME=$(PY_ENV)/cache \
		PYTHONPATH=square_glauber_python/src $(PYTHON_BIN) -m square_glauber.cli validate-exact

test-julia:
	$(JULIA) --project=aztec -e 'using Pkg; Pkg.test()'
	sh aztec/test/smoke_workflows.sh

test: test-python test-julia
