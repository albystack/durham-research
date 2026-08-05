#!/bin/bash
set -euo pipefail

source /etc/profile
module purge
module load julia/1.10.4

repository_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$repository_root"
mkdir -p .hamilton_env

# The repository Manifest was generated with Julia 1.12.6. Hamilton currently
# supplies Julia 1.10.4, so use a separate local environment and develop the
# package into it. This preserves aztec/Manifest.toml unchanged.
julia --project=.hamilton_env --startup-file=no -e '
    using Pkg
    Pkg.develop(path=abspath("aztec"))
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.precompile()
    Pkg.test("AztecDiamond")
'

echo "Hamilton Julia 1.10 environment is ready at $repository_root/.hamilton_env"
