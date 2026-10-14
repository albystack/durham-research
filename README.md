# Research

- **aztec/**: `aztec.jl` and `graph.pdf`.
- **square/**: `square.py` and `graph.pdf`.
- **old folders/**: the previous project files, including data, analysis, and notes. Git history is preserved in this repository.

The code files consolidate the simulation code used for the plotted results. The graphs use the completed campaign data; the default runs below are small examples and do not reproduce the full campaign or create the graphs.

## Run a small example

From this folder, with Julia 1.10+ and Python 3.11+ installed:

```sh
julia aztec/aztec.jl --size 8 --environments 2 --output /tmp/aztec-example
python3 -m pip install numpy
python3 square/square.py --size 8 --environments 2 --output /tmp/square-example
```

Choose fresh output directories for repeat runs. The model uses independent Gamma edge families with shapes 0.2 and 0.25 (scale 1). The plotted fits cover L = 16-128 and do not establish an asymptotic scaling law.

AI assistance covered implementation, testing, analysis, and cleanup. Original campaign source, data, and analysis remain in `old folders/research-original-2026-10-14/`.
