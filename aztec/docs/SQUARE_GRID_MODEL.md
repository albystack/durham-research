# Square-grid paired Temperley-dimer model

Status: implementation contract for the square-grid experiment requested by Professor Sunil Chhita.

## Scientific objective

For each random environment, sample two conditionally independent spanning trees, map both trees to generalized Temperley perfect matchings, and measure spatial dimer-height increments from the same fixed central cut. Estimate

\[
T_L(r)=\frac12\operatorname{Var}(\Delta H_1(r)-\Delta H_2(r))
\]

and

\[
D_L(r)=\operatorname{Cov}(\Delta H_1(r),\Delta H_2(r)).
\]

The first quantity isolates conditional sampling noise and the second isolates the shared-environment contribution.

## Finite graph and boundary

For integer `L >= 1`, the interior vertices are

\[
\Lambda_L=\{(x,y)\in\mathbb Z^2: \max(|x|,|y|)<L\}.
\]

A step from an interior vertex to a site with `max(abs(x),abs(y)) = L` is wired to one root. The geometric direction of every boundary edge is retained. This gives `V=(2L-1)^2+1` vertices and `E=4L(2L-1)` labelled primal edges.

The planar embedding is represented by the `4L^2` unit cells with lower-left corners in `[-L,L-1]^2`. The cell `(-L,-L)` is the distinguished outer face. The other `4L^2-1` cells are bounded face nodes in the Temperley graph.

## Environment laws

### Baseline

All four outgoing weights equal one.

### Directed site-i.i.d. model

At every interior site, north/east/south/west weights are independent. Opposite orientations of one geometric edge are generally different. The initial disorder law was mean-one Gamma:

\[
w_{x,d}\sim\operatorname{Gamma}(k,1/k),\qquad k=0.5.
\]

The robustness panel also supports Gamma `k=1`, mean-one lognormal weights,
and `Uniform(0,2)` weights. Overall deterministic rescaling of every edge
weight leaves the tree law unchanged.

### Undirected conductance comparison

Each unoriented geometric edge receives one conductance shared by both orientations. This is implemented as a separate model and is never pooled with the directed model.

### Scope relative to the random-bond dimer literature

The directed model is the square-grid spanning-tree environment studied in the
earlier LERW experiments and is the model requested in the supervisor email.
The Kenyon--Propp--Wilson correspondence makes it a genuine weighted dimer
measure on the associated Temperley graph: a directed primal weight is carried
by the link from a primal-vertex white node to the corresponding edge-node,
while the face-to-edge links have unit weight.

This structured disorder is **not** the canonical random-bond domino model in
which every dimer edge receives an independent random energy. Consequently, a
null squared-log result here is evidence about this directed
spanning-tree/Temperley law only. It neither contradicts nor tests the full
i.i.d. random-bond square-lattice dimer model used in the Zeng--Leath--Hwa and
Cardy--Ostlund super-roughness literature. The undirected-conductance model is
also a separate structured law and is reported separately.

All environment weights are deterministic functions of the public environment seed and a geometric site or edge key. Query order and Julia thread scheduling cannot change them.

## Weighted tree law

Every tree is oriented toward the wired root. For fixed environment `omega`,

\[
\mathbb P_\omega(T)\propto\prod_{(x\to y)\in T}w_{x\to y}(\omega).
\]

Wilson's algorithm uses the row-normalised transition probabilities. Since every rooted tree contains one outgoing edge from each non-root vertex, row normalisation multiplies all tree weights by the same environment-dependent constant.

## Exact generalized Temperley matching

For every primal edge, create one black edge-node. White nodes are all non-root primal vertices and all non-outer faces.

Given a rooted tree:

1. each primal vertex is matched to the black node of its outgoing tree edge;
2. non-tree primal edges form the complementary dual tree;
3. orient that dual tree toward the outer face;
4. each bounded face is matched to the black node of its outgoing dual-tree edge.

The code validates that every primal edge-node has exactly one owner and that the matching contains `(2L-1)^2` vertex matches and `(2L)^2-1` face matches.

## Exact spatial dimer-height increment

The Temperley graph is embedded at half-grid spacing:

- primal vertex nodes at even/even coordinates;
- horizontal edge-nodes at odd/even coordinates;
- vertical edge-nodes at even/odd coordinates;
- face nodes at odd/odd coordinates.

For separation `r`, define

\[
x_{\rm left}=-\lfloor r/2\rfloor,
\qquad x_{\rm right}=x_{\rm left}+r.
\]

The height is integrated along the horizontal dual segment at scaled height `1/2`, from the left probe to the right probe. It alternately crosses:

- a vertex-to-vertical-edge link, contributing `+1` when matched;
- an upper-face-to-horizontal-edge link, contributing `-1` when matched.

All crossed white nodes have degree four. There are equally many upward and downward crossings, so the deterministic `1/4` reference-flow contribution cancels exactly. The returned integer is therefore the exact dimer-height increment in the chosen convention, not a branch-winding proxy.

## Relative separations

Fractions are taken relative to the full square side `2L`:

```text
1/32, 1/16, 1/8, 1/4
```

Thus the integer separations are

\[
r_L(\rho)=\operatorname{round}(2L\rho).
\]

This gives distinct separations `1,2,4,8` already at `L=16`.

## Replica and bootstrap contract

For each `(model,L,sample_id)`, derive separate deterministic seeds for:

- the environment;
- replica 1 Wilson walks;
- replica 2 Wilson walks.

Both replicas share only the environment. Every separation for one replica is evaluated on the same complete tree and matching. Statistical resampling must keep both replicas and all fractions from one environment together.

The covariance-aware pooled sensitivity analysis estimates the joint `4 x 4`
covariance of the four fraction-level statistics at each order using common
environment-bootstrap indices. Block generalized least squares uses these
within-order covariances and treats different orders as independent blocks.

## Initial pilot

The frozen pilot configuration is:

```text
L = 16, 32, 64, 128, 256
baseline environments = 200, 200, 100, 50, 20
directed Gamma(k=0.5) environments = the same counts
replicas = 2
fractions = 1/32, 1/16, 1/8, 1/4
```

This is a correctness and detectability pilot. It is not an asymptotic proof.

## Primary files

- `src/SquareGrid.jl`
- `scripts/run_square_grid_campaign.jl`
- `scripts/merge_square_grid_batches.jl`
- `test/square_grid_runtests.jl`
- `configs/square_grid_smoke.csv`
- `configs/square_grid_pilot.csv`
- `../hpc/square_grid_smoke.slurm`
- `../hpc/square_grid_pilot.slurm`
- `../hpc/analyze_square_grid.slurm`

## References

1. R. W. Kenyon, J. G. Propp and D. B. Wilson, *Trees and Matchings*, Electronic Journal of Combinatorics 7 (2000), R25.
2. D. B. Wilson, *Generating Random Spanning Trees More Quickly than the Cover Time*, STOC 1996.
