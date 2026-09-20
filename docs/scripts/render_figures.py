"""Regenerate README illustrations from the repository's samplers.

Run `make figures` from the repository root. These are illustrative samples,
not measurements used for any scaling claim.
"""

from pathlib import Path
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch
import numpy as np

from square_glauber.environment import sample_environment
from square_glauber.glauber import run_attempts
from square_glauber.matching import DimerState


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "docs" / "figures"
COLORS = ("#247C85", "#76B6BA", "#E2AD57", "#C5694D")
INK = "#233746"


def save(fig, name):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT / f"{name}.svg", facecolor="white", metadata={"Date": None})
    svg = OUTPUT / f"{name}.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    fig.savefig(OUTPUT / f"{name}.png", dpi=180, facecolor="white")
    plt.close(fig)


def samples():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "tiling.csv"
        subprocess.run(
            ["julia", "--project=aztec", "--startup-file=no", "-e", """
            using AztecDiamond, Random, DelimitedFiles
            t = sample_tiling(Xoshiro(20260920), [fill(0.5, k, k) for k in 1:28])
            @assert validate_tiling(t).valid
            writedlm(ARGS[1], Int.(t), ',')
            """, str(path)], cwd=ROOT, check=True,
        )
        tiling = np.loadtxt(path, delimiter=",", dtype=int)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    fig.subplots_adjust(left=.04, right=.96, bottom=.10, top=.84, wspace=.22)
    ax = axes[0]
    side = tiling.shape[0]
    for row, col in zip(*np.nonzero(tiling)):
        i, j = row + 1, col + 1
        x, y = j - i, side + 1 - i - j
        horizontal = i % 2 == j % 2
        dx, dy = (2, 1) if horizontal else (1, 2)
        corners = np.array([(x-dx, y-dy), (x+dx, y-dy),
                            (x+dx, y+dy), (x-dx, y+dy)], dtype=float)
        # Keep the native diamond orientation so the two domains are distinct.
        points = corners
        color = COLORS[0 if i % 2 and j % 2 else
                       1 if i % 2 else 2 if not j % 2 else 3]
        ax.add_patch(Polygon(points, facecolor=color, edgecolor="white", linewidth=.22))
    extent = side + 3
    ax.set(xlim=(-extent, extent), ylim=(-extent, extent))
    ax.set_title("Aztec diamond", loc="left", color=INK, fontsize=19, pad=20)
    ax.text(.0, -.06, "Uniform domino shuffling · order 28", transform=ax.transAxes,
            fontsize=11, color=INK)

    L = 20
    state = DimerState.all_horizontal(L)
    environment = sample_environment(L, np.random.default_rng(20260920), "uniform")
    run_attempts(state, environment, np.random.default_rng(20260921), 500 * (L-1)**2)
    state.validate()
    ax = axes[1]
    for row, col in zip(*np.nonzero(state.horizontal_occupied)):
        ax.add_patch(Rectangle((col-.5, row-.5), 2, 1, facecolor=COLORS[0],
                               edgecolor="white", linewidth=.65))
    for row, col in zip(*np.nonzero(state.vertical_occupied)):
        ax.add_patch(Rectangle((col-.5, row-.5), 1, 2, facecolor=COLORS[3],
                               edgecolor="white", linewidth=.65))
    ax.set(xlim=(-1, L), ylim=(L, -1))
    ax.set_title("Glauber dynamics", loc="left", color=INK, fontsize=19, pad=20)
    ax.text(.0, -.06, "20 × 20 grid · after 500 uniform-weight sweeps", transform=ax.transAxes,
            fontsize=11, color=INK)
    for ax in axes:
        ax.set_aspect("equal")
        ax.axis("off")
    save(fig, "samples")


def flip():
    fig, ax = plt.subplots(figsize=(10, 2.5))
    fig.subplots_adjust(left=.02, right=.98, bottom=.05, top=.95)
    ax.set(xlim=(0, 10), ylim=(-.4, 2.35), aspect="equal")
    ax.axis("off")
    for y in (0, 1):
        ax.add_patch(Rectangle((.5, y), 2, 1, facecolor=COLORS[0],
                               edgecolor="white", linewidth=2))
    for x in (7.5, 8.5):
        ax.add_patch(Rectangle((x, 0), 1, 2, facecolor=COLORS[3],
                               edgecolor="white", linewidth=2))
    ax.add_patch(FancyArrowPatch((3.1, 1), (6.9, 1), arrowstyle="<->",
                                 mutation_scale=22, linewidth=1.5, color=INK))
    ax.text(5, 1.6, "Choose a face. Redraw its two dimers.", ha="center", fontsize=12, color=INK)
    ax.text(5, .35, "Equal weights: each orientation has probability ½.",
            ha="center", fontsize=10, color=INK)
    save(fig, "local-flip")


if __name__ == "__main__":
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                         "svg.hashsalt": "dimer-research"})
    samples()
    flip()
