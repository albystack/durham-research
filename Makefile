PYTHON ?= python3.13
PY_ENV ?= square_glauber_python/.venv
PYTHON_BIN := $(PY_ENV)/bin/python

.PHONY: help setup-python test-python validate-python test-julia test

help:
	@echo "make setup-python PYTHON=python3.13  Create the Python environment"
	@echo "make test-python                     Run square-grid Python tests"
	@echo "make validate-python                 Run exact tiny-grid validation"
	@echo "make test-julia                      Run the Julia/Aztec test suite"
	@echo "make test                            Run both project test suites"

setup-python:
	$(PYTHON) -m venv $(PY_ENV)
	$(PYTHON_BIN) -m pip install --upgrade pip
	$(PYTHON_BIN) -m pip install -e './square_glauber_python[test,accel]'

test-python:
	MPLCONFIGDIR=$(PY_ENV)/matplotlib-cache XDG_CACHE_HOME=$(PY_ENV)/cache \
		PYTHONPATH=square_glauber_python/src $(PYTHON_BIN) -m pytest square_glauber_python/tests

validate-python:
	MPLCONFIGDIR=$(PY_ENV)/matplotlib-cache XDG_CACHE_HOME=$(PY_ENV)/cache \
		PYTHONPATH=square_glauber_python/src $(PYTHON_BIN) -m square_glauber.cli validate-exact

test-julia:
	julia --project=aztec -e 'using Pkg; Pkg.test()'
	sh aztec/test/smoke_workflows.sh

test: test-python test-julia
