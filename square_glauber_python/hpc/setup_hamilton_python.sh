#!/bin/bash
set -euo pipefail

repository_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$repository_root"

if [[ -z "${PYTHON_BIN:-}" ]]; then
  module load python/3.12.6
  PYTHON_BIN=python3
fi
"$PYTHON_BIN" -c 'import sys; assert sys.version_info >= (3, 11), sys.version'

if [[ ! -x .hamilton_python/bin/python ]]; then
  "$PYTHON_BIN" -m venv .hamilton_python
fi

mkdir -p .hamilton_python/matplotlib_cache
export MPLCONFIGDIR="$repository_root/.hamilton_python/matplotlib_cache"

.hamilton_python/bin/python -m pip install --upgrade pip setuptools wheel
.hamilton_python/bin/python -m pip install '.[test,accel]'
PYTHONPATH=src .hamilton_python/bin/python -c \
  'from square_glauber.accelerated import require_numba; require_numba(); print("Numba acceleration available")'
PYTHONPATH=src .hamilton_python/bin/python -m pytest

echo "Hamilton Python environment is ready at $repository_root/.hamilton_python"
