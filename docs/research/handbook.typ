// Typst source with repository-local figures; no network packages required.
// Build from the repository root:
// typst compile --root . docs/research/handbook.typ output/pdf/research-handbook.pdf
#set document(title: "Random environments, trees and dimers: research handbook", author: "Alberto Rescigno", date: datetime(year: 2026, month: 10, day: 10))
#set page(paper: "a4", margin: (x: 23mm, y: 22mm), numbering: "1", header: context { if counter(page).get().first() > 1 { text(size: 8pt, fill: rgb("667085"))[Alberto Rescigno · Research reference · 10 October 2026] } })
#set text(font: "New Computer Modern", size: 10.5pt)
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): it => { v(0.8em); text(fill: rgb("164E63"), it); v(0.3em) }
#show raw: set text(size: 8.5pt)
#show table: set par(justify: false)
#let Var = math.op("Var")
#let Cov = math.op("Cov")
#let Aztec = math.op("Aztec")
#let square = math.op("square")
#let note(body) = block(fill: rgb("EFF6F8"), inset: 10pt, radius: 3pt, width: 100%, body)
#let source(path) = block(text(size: 8pt, fill: rgb("52616B"))[Repository source: #raw(path)])

#align(center)[
  #v(25mm)
  #text(size: 27pt, weight: "bold", fill: rgb("164E63"))[Random environments,\ trees and dimers]
  #v(8mm)
  #text(size: 17pt)[Research handbook]
  #v(8mm)
  Alberto Rescigno\ Durham University\ Supervisor: Sunil Chhita
  #v(10mm)
  Living reference · 10 October 2026
]
#v(15mm)
#note[Purpose: collect the definitions, mathematical reasoning, model conventions, computational methods, experimental history, numerical evidence and open questions needed to resume this project without relying on memory.]
This is a synthesis of the available repository and the July–September correspondence supplied in conversation. It is not a complete archive of every discussion or attachment: the referenced screenshots and handwritten notes have not all been supplied. Historical numerical claims are labelled by provenance. No theorem is inferred from a numerical fit. The completed corrected Aztec refinement is recorded separately from the still-unresolved domain contrast.

#pagebreak()
#block[
  #set text(size: 9pt)
  #outline(title: [Contents], indent: auto, depth: 2)
]
#pagebreak()

= The question and the current position
#note[Read the final chapter for the completed 258,048-pair experiment. Both Aztec and ordinary square geometries show positive fitted squared-log contributions over sizes 16–128. Higher-cutoff fits remain inconclusive. Earlier chapters preserve the dated experiments, sampler correction and control investigation.]
We study fluctuations in lattice models with random environments. The central question is whether the disorder-induced part of a height or winding fluctuation grows approximately like $log L$ or has a positive contribution proportional to $(log L)^2$, and why experiments on the Aztec diamond and on square regions differ.

The project developed through uniform spanning trees (USTs), loop-erased random walks (LERWs), random transition environments, Aztec domino shuffling, a Temperley representation, direct dimer Glauber dynamics, finite-volume Kasteleyn checks and finally coupling from the past (CFTP).

Historically, the strongest reported positive numerical evidence was in *spatial height-increment disorder covariance* on the Gamma-disordered Aztec diamond. Its interpretation is now qualified by the final audit-update chapter. Central height, a two-replica difference and spatial disorder covariance are different observables; success or failure for one does not transfer automatically to another.

The September all-edge Gamma square experiment did not show a robust positive quadratic-log coefficient. The completed October campaign uses the original two-family Aztec Gamma law instead. It shows positive fitted curvature in both geometries over sizes 16–128, with uncertainty at higher cutoffs. These are different models and observables; preserve both dated findings.

== Sunil's latest guidance
The supplied reply endorses the understanding of weighted heat-bath dynamics and continuation with CFTP. Any independent mean-one edge weights are acceptable as a model, but the weights used in the positive Aztec experiment are preferred for the comparison. The working explanation is a special role for Aztec boundary conditions. This is a hypothesis to investigate, not an established conclusion of the existing data.

Large inverse-Kasteleyn calculations are not the preferred main direction. Existing finite-volume calculations remain useful validation. Sunil also requests an honest account of AI contributions. A contribution record now documents confirmed assistance; historical attribution still requires Alberto’s account. This handbook does not reconstruct authorship from guesses.

Authenticated Hamilton access was used for the bounded validation, feasibility and recovery jobs. The corrected Aztec refinement, matched two-domain campaign and separate all-one confirmation are complete. The fresh control primary contrast is $0.081$ with interval $[-0.282,0.437]$ and permutation $p=0.6615$: the original anomaly was not reproduced at this precision. The completed larger campaign gives Aztec $c=0.17994$ and ordinary square $c=0.15882$, with positive primary intervals in both. See the final chapter for the finite-size qualifications and scientific review.

= Dictionary and notation
#table(columns: (1fr, 3fr), inset: 6pt, stroke: 0.4pt,
  table.header([Symbol / term], [Meaning]),
  [$G=(V,E)$], [Finite graph with vertices and undirected edges.],
  [$omega$], [One frozen realisation of all random edge weights.],
  [$M$], [Perfect matching: every vertex is incident to exactly one selected edge.],
  [$w_e$], [Positive weight of edge $e$, fixed while replicas are sampled.],
  [$Z_omega$], [Partition function: sum of matching weights.],
  [$h(f)$], [Height at face $f$, with a specified common boundary reference.],
  [$X$], [Observable: central height, spatial increment, winding, etc.],
  [$X_1,X_2$], [Two conditionally independent copies given the same environment.],
  [$T$], [Expected conditional variance, also called connected or thermal component.],
  [$D$], [Variance of the conditional mean, also called disorder component.],
  [$rho$], [Requested relative separation; actual integer endpoints must also be stored.],
  [$c$], [Coefficient of the squared logarithm in a finite-size fitted model.],
  [Quenched], [Environment held fixed; conditional distribution given that environment.],
  [Annealed], [Average over both environments and samples from their conditional laws.],
  [Replica], [A sample or trajectory sharing weights with another but using independent sampling randomness.],
  [Censored], [A requested CFTP stream has reached its cap without certification. It is not an equilibrium sample.],
)

== Size conventions that must not be mixed
Python square code uses an even $L$ for an $L times L$ grid of *vertices*, with $(L-1)^2$ bounded faces and $2L(L-1)$ edges. Julia's retained square module uses a size parameter for a $2L times 2L$ vertex grid. Aztec order $n$ has $n(n+1)$ dominoes and a $2n times 2n$ Boolean encoding; that encoding is not an ordinary square board.

One cannot pool data merely because both tables contain a column named `L` or `order`. Record physical lattice spacing, domain extent, endpoint positions and number of steps along a probe. A constant size rescaling changes intercept and linear-log coefficient but not the quadratic coefficient in an unrestricted quadratic in the log; different windows, rounded separations and different observables create additional differences.

= Random walks, spanning trees and the original experiments
== Walks, loop erasure and Wilson's algorithm
A nearest-neighbour walk moves along graph edges with transition probabilities $P(x,y)$. A simple random walk chooses uniformly among neighbours. Chronological loop erasure removes the segment between repeated visits as each loop is formed, leaving a simple path.

A spanning tree is a connected acyclic subgraph containing every vertex. Wilson's algorithm grows a tree by loop-erasing walks until they hit an existing rooted tree. Its target depends on the transition kernel: symmetric conductances and a directed transition environment must not be conflated. Ordinary boundaries retain the finite graph; wired boundaries identify appropriate boundary vertices with a root.

The original winding statistic counts left turns minus right turns along the retained path, in quarter-turn units. Angles in radians require multiplication by $pi/2$, and variances then by $(pi/2)^2$. Verify the specific archived routine before converting historical values.

== Fixed and refreshed environments
In a fixed site environment, weights determine transition probabilities after normalisation. Making only some directions random can introduce drift. The early construction with north/east weights one and south/west random weights was not a symmetric all-edge disorder model. Drawing each outgoing direction from the same law gives distributional directional symmetry, but does not by itself make a realised directed walk reversible.

Refreshing weights at each step changes the stochastic model. It is not a faster way to sample a quenched model. In particular, the shared-environment replica decomposition needs an explicitly defined common random environment; it cannot be imported unchanged into a model whose weights are continually refreshed.

== Temperley correspondence: why the model matters
Temperley-type constructions relate suitable spanning trees and perfect matchings on associated bipartite graphs. The graph, boundary conditions and weight transformation are part of that correspondence. Independent random transition or tree weights need not become independent identically distributed dimer-edge weights. The large archived square/Temperley null results are valuable for their defined model, but do not establish the same claim for the direct random-bond dimer model.
#source("archive/experiments/square_grid/; archive/random_walk/; archive/history/docs/RESULTS.md")

= Weighted dimers and domino tilings
A domino tiling can be represented as a perfect matching of the adjacency graph of board cells. For positive weights the finite-volume Gibbs law is
$ pi_omega(M) = frac(1,Z_omega) product_(e in M) w_e, quad Z_omega = sum_(M) product_(e in M) w_e. $
The partition function is finite and positive whenever at least one matching exists. Randomness enters twice: first draw $omega$, then draw $M$ under $pi_omega$. Shared weights do not mean shared update coins.

== Gamma distribution and conventions
For shape $k>0$ and scale $theta>0$, a Gamma variable has mean $k theta$ and variance $k theta^2$. Scale is the reciprocal of rate. Always store which convention is used.

The September square environment assigns every undirected edge an independent Gamma variable with shape $k=1$ and scale $1/k=1$. Its general all-edge sampler uses scale $1/k$, hence mean one for every shape.

The original Aztec experiment uses independent families $a_(i,j) tilde Gamma(0.2,1)$ and $b_(i,j) tilde Gamma(0.25,1)$, with the other two families one. Their means are 0.2 and 0.25. Small shape gives substantial mass near zero; strictly positive mathematical weights can still present floating-point and computational difficulties.

== Gauge transformations and normalisation
If each edge $(u,v)$ is replaced by $w'_(u,v)=g_u g_v w_(u,v)$ with fixed positive vertex factors, every matching weight acquires the same factor $product_(v in V) g_v$. The Gibbs probabilities are unchanged. This is a gauge transformation.

Multiplying selected orientation families by arbitrary constants is not automatically such a transformation. Check alternating-cycle products, or construct the vertex factors explicitly. Scaling both Gamma families or normalising them separately should not be assumed equivalent. The new comparison preserves the original scale-one law.
#source("square_glauber_python/src/square_glauber/environment.py; aztec/src/AztecDiamond.jl")

= Height functions and observables
== The Python height convention
Vertices have zero-based array coordinates $(i,j)$ and are white when $i+j$ is even. Orient primal edges from white to black. Crossing such an edge from its geometric left to its right adds $+1$ if unoccupied and $-3$ if occupied. Reversing the crossing changes the sign. Around an interior vertex the matching constraint makes the circulation zero.

The code fixes exterior height zero and enters bounded face $(0,0)$ across the top-left horizontal edge. It then integrates over interior dual edges, checking that alternative paths give the same integer value. Boundary vertices have degree two or three, so identifying all exterior sectors and imposing the degree-four rule everywhere would be inconsistent. The fixed cut is part of the observable definition.

Do not subtract a different interior reference for each sample when measuring central height: that changes the random variable. A spatial difference cancels additive constants naturally.

== Local flips and height order
A flippable square has two opposite occupied edges. Replacing them by the other pair changes the selected face height by four in magnitude. All other bounded-face heights remain unchanged. Legal normalised heights form a partially ordered set under pointwise comparison. On the simply connected rectangle there are true minimum and maximum tilings.

All-horizontal and all-vertical tilings are convenient initial states, but are not generally those extrema. The code constructs bounds by exhausting directed downward or upward flips. Exact enumeration already distinguishes the convenient starts from true extrema on a $4 times 4$ grid.

== Central height and spatial increments
The Python central bounded face is $((L-2)/2,(L-2)/2)$. The September increment is rightward from that face. Nominal offsets are $L/32,L/16,L/8,L/4$, rounded half upward; zeros, duplicates and out-of-domain offsets are omitted.

Aztec probes use two face columns separated by an integer $r$, centred as closely as parity permits on the central row. In the Julia staggered table, `left = fld(n-r,2)+1`, `right = left+r`, `row = n+1`. Table coordinates need a geometric conversion before being compared with square coordinates. A separation parameter can mean total endpoint separation, not half-separation.

For the matched comparison use symmetric endpoints and record both coordinates, actual displacement and direction. Small-grid probe collapse must be visible in metadata. Never count duplicate observables as independent measurements.
#source("square_glauber_python/src/square_glauber/height.py; square_glauber_python/src/square_glauber/geometry.py; aztec/src/AztecDiamond.jl")

= Glauber heat bath and its target
For clockwise top/right/bottom/left weights $(a,b,c,d)$ on a flippable face, the two local matching weights are $a c$ and $b d$. Conditional on the outside matching,
$ p_H = frac(a c,a c+b d), quad p_V=frac(b d,a c+b d). $
Choose uniformly from *all* bounded faces, including nonflippable faces. If flippable, redraw an orientation with these probabilities. Redrawing the current orientation is a valid unchanged update. Otherwise do nothing. One sweep means $(L-1)^2$ attempted updates, not accepted flips.

Detailed balance follows because two matchings differing on one face have weight ratio $a c/(b d)$ and the reverse heat-bath probabilities give exactly that ratio. Positive weights, flip connectivity and self-loops supply irreducibility and aperiodicity on this finite rectangle.

Evaluate probabilities through log weight ratios to avoid product underflow. Accelerated code must consume the same supplied face and uniform arrays as the reference implementation. Speed must not change the stochastic kernel.

A long ordinary trajectory is not automatically an equilibrium sample. Autocorrelation reduces effective sample size; agreement of two starts is a diagnostic, not a proof of mixing. Parallel tempering maintains a larger joint stationary law but still needs effective communication between temperatures. The historical exchange-clock bug preserved individual kernels while preventing intended exchanges, demonstrating why correct stationarity alone is insufficient.
#source("square_glauber_python/src/square_glauber/glauber.py; accelerated.py; aztec/src/GlauberSquareGrid.jl")

= Coupling from the past
== Mathematical mechanism
Choose random update maps indexed by negative time. Start every possible state at time $-T$ and apply the same maps up to zero. If every trajectory agrees, their common state is a sample from the stationary distribution. Monotonicity allows the implementation to track only the true minimum and maximum: all other states remain between them.

The heat-bath threshold at a selected face depends only on frozen weights. A common face and common uniform variable therefore define the coupling. The local monotonicity audit enumerates both parities, admissible local occupancy masks and relevant height gaps. The saved mathematical contract records 496 admissible local cases preserving order.

If the bounds have not coalesced, prepend an older independent block, restart at the extrema and replay *all existing later maps unchanged*. Merely running forward until two chains meet is not CFTP. Redrawing the old maps after failure also breaks the construction.

== What certification does and does not mean
A returned state has no arbitrary burn-in truncation, conditional on the implemented random-map model. The implementation still uses pseudorandom generators and floating-point probabilities. Certification is not an asymptotic scaling theorem, nor a guarantee that a runtime-censored collection represents the environment distribution without bias.

Two separately seeded streams share one frozen environment. A capped stream returns no state. Preserve its environment, seed, initial horizon and random-block convention. Extend the same history; do not replace it with an easier seed or environment. Power-of-two horizons are the safe standard for repeatable doubling. Arbitrary final partial blocks require care when extending a history.

Runtime and coalescence horizons are diagnostics. A horizon is not an independent sample and is not automatically a mixing-time estimate. Report the censored fraction by size and model; runtime quantiles among completed streams are conditional on completion.
#source("square_glauber_python/PERFECT_SAMPLING.md; src/square_glauber/perfect.py; restart.py; extrema.py; order_audit.py")

= Two replicas and variance decomposition
Let $mu_omega = E[X | omega]$ and $v_omega = Var(X | omega)$. Conditional independence gives
$ E[X_1 X_2 | omega] = mu_omega^2. $
Taking expectations and subtracting the common squared mean gives
$ D = Cov(X_1,X_2) = Var(mu_omega). $
Also $E[X_1-X_2|omega]=0$ and $Var(X_1-X_2|omega)=2v_omega$, so
$ T=frac(1,2) Var(X_1-X_2)=E[v_omega]. $
The law of total variance gives
$ Var(X)=T+D. $
These identities apply to a height, increment or other square-integrable observable under the stated pairing. They do not require Gaussianity.

With $N$ environment pairs, use the sample covariance with denominator $N-1$ for $D$ and half the sample variance of differences for $T$. Finite-sample covariance can be negative even though the population disorder variance is nonnegative. Do not clip negative estimates before fitting. The sum of these two estimators equals the average of the two replica sample variances algebraically.

An all-one environment has no random conditional mean, so $D=0$. Two independent configurations still have nonzero conditional fluctuations and finite-sample covariance noise. Replica correlation and marginal means/variances are useful diagnostics, not replacements for the decomposition.

= Statistical analysis and interpretation
== Whole-environment bootstrap
Resample independent environments within each size/model group. Keep both replicas and every spatial observable of each selected environment together. Recompute covariances and then refit. Resampling replicas separately destroys the shared-environment effect; treating different separations as independent underestimates uncertainty.

The frozen September plan uses 10,000 bootstrap replicates and seed 2026090902. Its primary observables are the original rightward spatial increments, with central height secondary. New symmetric analyses are labelled exploratory for old data and prospectively specified for new data; they do not retroactively replace that plan.

== Logarithmic models
Compare
$ D(L)=a+b log L $
with
$ D(L)=a+b log L+c(log L)^2. $
For pooled spatial analysis each relative separation has its own intercept and log slope, with one common $c$. The sign of $c$ is unconstrained. Use actual separations when they are the declared independent variable. Inspect residuals, conditioning, leverage and fixed lower-size cutoffs.

A power fit $C(log L)^p$ is a different family; it cannot be applied to negative covariance estimates and is not interchangeable with an additive quadratic. A positive point estimate, bootstrap fraction above zero, confidence interval, BIC comparison and held-out prediction error answer different questions. A bootstrap fraction is not a Bayesian posterior probability.

Unweighted SSE-derived AIC/BIC in the archived scripts are descriptive curve-comparison metrics, not a complete likelihood for the correlated raw measurements. Block GLS can account for covariance across separations, but requires stable covariance estimation. Near-singular designs and too few sizes should fail clearly rather than return spurious precision.

== Comparing domains
A coefficient significant in one model and nonsignificant in another does not prove the coefficients differ. After matching conventions, estimate $Delta c=c_(Aztec)-c_(square)$ with independent whole-environment bootstrap samples from the two campaigns. Report individual and difference intervals, common physical windows where possible and explicitly non-overlapping windows otherwise.

Changing a height scale by a factor $s$ scales its covariance and fitted $c$ by $s^2$. A sign reversal alone does not change covariance. Probe direction, local slope, frozen regions and distance to the boundary may matter even after the edge law is matched. A finite-volume contrast alone does not isolate boundary conditions as its unique cause.

== Controls, multiplicity and missingness
Individual 95% intervals will sometimes exclude zero under the null. Here the intervals share measurements, so a binomial calculation assuming independence is inappropriate. Use a joint test preserving each replica's vector across separations, such as permuting entire second-replica vectors within each size and evaluating a maximum standardised covariance. This tests pairing dependence under the all-one model; it does not prove the single-replica law is correct.

Check duplicate identifiers, seeds, state hashes, full probe coverage and aggregation joins. Equal configurations can occur legitimately, especially on tiny grids; a duplicated seed is a different issue.

Omitting a slow environment is informative missingness. Standard bootstraps of completed environments do not account for it. Preserve unresolved cases and report support bounds or explicitly limited conclusions. Do not impute a missing covariance contribution as zero.

= Experimental history and numerical ledger
Unless stated otherwise, this section summarises retained historical reports; it does not claim all raw experiments were rerun on 29 September. The authoritative dated files remain unchanged.

== April to July: trees and walks
The April UST implementation recorded 37 passing tests, ordinary and wired boundaries, seeded walks and structural checks. Early asymmetric random environments produced strong northeast drift: at side 256, mean boundary hit coordinates were approximately $(0.894L,0.889L)$ for Gamma shape 1 and $(0.916L,0.918L)$ for shape 0.5. This motivated more balanced models.

The strict one-walk-per-environment campaign recorded 379,000 walks/environments across 15 distributions, sizes 16 to 4096, with baseline and Gamma extended to 8192. Fitted log-power exponents were approximately 1.01 to 1.16; BIC preferred ordinary log in 14 of 15 cases. Paired LERW experiments used 379,000 environments and 758,000 walks; reported exponents ranged 1.036 to 1.187, with largest bootstrap upper bound 1.260. These are winding results for their particular models.

Refreshed-weight pilots used about 46,000 runs through size 1024. Reported exponents were 1.064, 1.123, 1.069 and 1.105 for baseline, Gamma, lognormal and Pareto cases. A shape-0.5 extension reached 8192; the path-length exponent was 1.2518 with interval [1.2456, 1.2573], consistent with the familiar five-quarters benchmark. Its winding exponent 1.134 had interval [0.957, 1.265].
#source("archive/history/docs/RESULTS.md, sections 1-4; archive/random_walk/results/RESULTS_INDEX.md")

== Aztec: from central height to spatial disorder covariance
The retained central-height dataset contains 35,536 rows, orders 16 through 1300. The paired centre dataset has 28,304 independent environments. Metadata records shapes 0.2 and 0.25, scale one and deterministic seed derivation. The historical email central-height interval contains a transcription error, “[90.894, 1.032]”; do not silently reuse a guessed correction.

For the two-tiling central-height difference, the ledger reports log-power exponent 0.814 with interval [0.700, 0.902]. This isolates the conditional component and therefore does not rule out disorder curvature.

The original spatial analysis uses 6,900 Gamma environments and 5,500 all-one pairs, with orders 128, 192, 256, 384, 512, 700, 900, 1100 and 1300. Four symmetric relative separations share a common quadratic coefficient.
#table(columns: (2.2fr, 1fr, 1.5fr), inset: 6pt, stroke: 0.4pt,
 [Fit], [$c$], [95% bootstrap interval],
 [Disorder, unweighted], [0.5843], [[0.1957, 0.9520]],
 [Disorder, precision weighted], [0.4948], [[0.1932, 0.7977]],
 [Conditional, unweighted], [-0.2013], [[-0.5410, 0.1411]],
 [Uniform marginal control], [0.00882], [[-0.2623, 0.2835]],
)
The original unweighted disorder comparison gives BIC difference 10.848 in favour of the quadratic extension; held-out two-largest-order RMSE improves from 0.886 to 0.588. Individual separation coefficients are approximately 0.566, 0.338, 1.045 and 0.389, but the four probes are not four independent experiments.

Later block-GLS broad-window coefficients are 0.5019 for the original law and 0.6077 for stronger disorder, with intervals excluding zero. Weak and symmetric laws are less convincing. Raising the lower cutoff to 384 or 512 substantially widens intervals and weakens the original positive conclusion. This is parameter-sensitive finite-size evidence, not universality or proof.
#source("aztec/data/DATASETS.md; archive/history/results/aztec/spatial/spatial_analysis_report.md; archive/history/results/aztec/hamilton_20260809/spatial_gls_campaign_comparison.md")

== August: square Temperley and direct dimers
The archived structured Temperley campaign reached size 6144 across five disorder laws. Fourteen of fifteen reported window/model intervals included zero, and fourteen BIC comparisons preferred ordinary log. Those are not automatically the independent-edge dimer law.

The direct Julia dimer campaign retained 960 Gamma(shape 0.5, scale 2) environments and 352 controls. Its size labels 2 to 20 correspond to vertex sides 4 to 40. The full-range central disorder coefficient was 0.108 with interval [-0.545, 0.794]. Some larger-size trajectories showed mixing concerns.

A Kasteleyn replay of the same 960 environments removed the MCMC mixing assumption for central moments. The full-range disorder coefficient was 0.130 with interval [-0.552, 0.785]. One ill-conditioned environment needed a 256-bit fallback; maximum retained solve residual was reported as $7.39 times 10^(-12)$. The aggregate null conclusion persisted. This does not independently settle the spatial-increment question.
#source("aztec/data/glauber_square_grid_20260822/DATASET.md; aztec/data/glauber_square_grid_kasteleyn_20260822/DATASET.md")

== September: certified square campaign
The frozen Python design has sizes 8, 12, 16, 20, 24, 28 and 32; 18,200 Gamma(shape 1) environments; and 5,600 all-one pairs. The saved September snapshot contains 23,799 of 23,800 requested pairs. The missing case is `gamma-L32-env001015`; both streams remained unresolved through 1,073,741,824 sweeps.

Gamma central disorder covariance rises from about 0.714 at side 8 to 1.414 at side 32. The full-range linear-log slope interval is [0.404, 0.718]. The central quadratic coefficient is 0.0446 with interval [-0.3444, 0.4254]; the pooled spatial coefficient is -0.2244 with interval [-0.6370, 0.1855]. Central coefficients become -0.4354 and -1.5176 at lower cutoffs 16 and 20. These are unstable finite-size curvature estimates.

The missing pair is one of 1,200 requested Gamma pairs at side 32. Extremal central heights range from -33 to 31. The saved conservative support calculation allows completed side-32 covariance from 0.5593 to 2.2684 and all-size central curvature approximately -0.8191 to 0.9077. These are worst-case bounds, not probable outcomes.

The 29 September table review independently confirmed the pair count and found 6 of 30 individual all-one covariance intervals excluding zero. The affected entries are centre at 24; spatial at (16, 1/8), (20, 1/4), (24, 1/4), (28, 1/8), (32, 1/4). Their dependence requires a joint audit.
#source("archive/history/square_glauber_python/reports/provisional_20260915/REPORT.md; input/center_measurements.csv; results/bootstrap_component_intervals.csv")

= Kasteleyn calculations as an independent check
For a finite planar bipartite graph, a correctly signed or complex-phased bipartite Kasteleyn matrix encodes the weighted matching partition function. With a consistent convention an edge probability is $K_(b,w) (K^(-1))_(w,b)$, and joint probabilities follow determinant identities. A height difference is an affine signed sum of edge occupations along a dual path, so conditional means and variances can be calculated from occupation moments.

Random weights do not prevent a finite-volume calculation for a fixed environment. They do create conditioning and cost issues; interpreting an averaged random inverse or an asymptotic kernel is a further problem. Numerical residuals, positivity/realness checks and higher precision fallbacks matter. Keep exact enumeration as the independent smallest-grid reference.

= Software and reproducibility map
#table(columns: (1.2fr, 2.3fr), inset: 6pt, stroke: 0.4pt,
 [Location], [Purpose],
 [`aztec/src/AztecDiamond.jl`], [Gamma weights, reduction, shuffling, height observables and rendering.],
 [`aztec/src/GlauberSquareGrid.jl`], [Retained Julia direct-dimer and Kasteleyn reference.],
 [`aztec/scripts/`], [Resumable campaigns, merging, analysis and plotting.],
 [`square_glauber_python/src/`], [Environment, matching, height, heat bath, extrema, CFTP and campaign logic.],
 [`square_glauber_python/tests/`], [Exact, mathematical, restart and workflow tests.],
 [`docs/METHODS.md`], [Short cross-project conventions.],
 [`archive/RELOCATION.csv`], [Mapping from old file locations to archived locations.],
 [`archive/local_runs/`], [Ignored historical local outputs; availability must be checked.],
)

Record revision, dependency versions, command, configuration, environment identifier, environment seed, replica seeds and hashes. Persist weights before sampling and certify each replica separately. Publish paired rows only when both states are certified. Preserve failures and censored records alongside successful observations.

Use unique seed namespaces for new campaigns, including model and parameter settings. Never reuse old IDs to disguise a new law. Preserve existing frozen manifests and analysis plans. Generated output is not a replacement for source data.

For historical reproduction, archived commands can still name old paths. Resolve them with relocation tables. The September snapshot is a local dated copy, not a live Hamilton status query. There is no claim that external jobs have been checked in this session.

= The matched-comparison programme
+ Define a single full-lattice edge-family convention, map it explicitly to the Aztec encoding and restrict it to square regions. Preserve scale one and original shapes.
+ Define symmetric physical probes, with coordinates, directions and height scales. Check both exact small configurations and the conversion between implementations.
+ Audit all-one controls jointly; check seed uniqueness, pairing joins and bootstrap blocks. Reanalyse old states only when their certified state files are available.
+ Use a small predetermined feasibility pilot. Retain every attempted environment and both stream outcomes. Measure completion rate, runtime tails and horizon tails by size and law.
+ Choose any production size/sample budget after feasibility, then freeze its primary observable, cutoffs, fit method, uncertainty method and stopping rules before examining its scaling result.
+ Compare domain coefficients directly, with appropriately blocked independent resampling and explicit physical windows. Add a small-order Aztec bridge if useful.
+ Log into Hamilton with the user's involvement before large runs. AI contribution reconstruction follows this implementation stage.

The all-one control separates disorder from thermal noise. A deterministic two-family control with weights 0.2 and 0.25 additionally separates the random environment from the bias of its family means. These answer different questions and should have distinct model names.

A local pilot cannot establish the large-size law. If censoring becomes prohibitive, report that computational limitation and reconsider the method; do not condition the scientific sample on quick coalescence.

= Questions to revisit and common mistakes
- Are the two compared edge laws identical after the stated coordinate mapping?
- Does the height observable retain a common reference and the same fluctuation scale?
- Is $r$ total separation, half-separation, table-column offset or physical lattice distance?
- Are probes in comparable parts of the liquid region, with comparable local slopes?
- Are all requested environments represented, including censored cases?
- Does the reported uncertainty preserve replica and spatial dependence?
- Are apparent curvature differences stable to the prospectively chosen windows?
- Is a claim about a finite interval being mistakenly described as an asymptotic theorem?
- Has a weight refresh, active-face selection or restart shortcut changed the model?
- Is an old historical status being mistaken for the current external campaign status?

= Reading and source trail
Primary mathematical references:
- Duits and Van Peski, _The Gamma-disordered Aztec diamond_, arXiv:2512.03033 (2025). Definition 1.1 and Figure 1 specify the biased model; equations 1.22 and 1.23 describe the reduction/creation convention used by the Julia implementation. #link("https://arxiv.org/abs/2512.03033")[arXiv record].
- Propp and Wilson, _Exact sampling with coupled Markov chains and applications to statistical mechanics_ (1996). #link("https://www.math.cmu.edu/~af1p/Teaching/MCC17/Papers/propp_wilson.pdf")[Paper].
- Desreux and Rémila, _An optimal algorithm to generate tilings_ (2005). #link("https://www.sciencedirect.com/science/article/pii/S1570866705000092")[Publisher record].

Start with the root README and methods notes. Then read the Aztec implementation guide, the Python project guide, the perfect-sampling contract, the frozen paired-covariance analysis plan, the September provisional report and the historical results ledger. Their exact locations are listed in the software map and repository source notes. The historical ledger explicitly warns that many early numbers were distilled from correspondence rather than re-derived from raw data.

== Maintenance rule
Update this source when a definition, model or verified result changes. Record date, provenance and limitations; retain the old dated evidence. Rebuild the PDF with Typst. Do not edit historical reports to make them appear prospectively consistent with new decisions. Supporting implementation and validation results for the matched comparison are documented separately in `docs/history/MATCHED_COMPARISON.md` and should be incorporated here when verified.

#pagebreak()
= Implementation and audit update: 29 September 2026
#note[*This section supersedes the earlier interpretation of the historical Aztec benchmark.* The newly implemented independent distribution audit finds an array-orientation mismatch in the legacy Gamma sampler. Historical coefficients remain a record of what was computed, but should not currently be presented as validated evidence for the stated independent red/yellow edge-weight law.]

== What the new exact check tests
The old tests check the recurrence, reproducibility and geometric validity of tilings. All can pass while the probability distribution on those tilings is wrong for the intended edge environment. The new audit enumerates every possible creation-coin path through orders two and three, adds the probabilities of paths reaching the same tiling and compares the result with independently computed products of occupied-edge weights.

The paper's red edge corresponds to matrix entry `(2i,2j-1)`, yellow to `(2i,2j)`, and the two other entries in that block have weight one. For two fixed non-symmetric test environments:
#table(columns: (0.7fr, 0.7fr, 1.5fr, 1.5fr), inset: 6pt, stroke: 0.4pt,
 [Order], [Tilings], [Legacy maximum probability error], [Corrected maximum error],
 [2], [8], [0.09305996], [$5.55 times 10^(-17)$],
 [3], [64], [0.04069997], [$4.16 times 10^(-17)$],
)
The corresponding legacy total variation discrepancies are 0.15661 and 0.25658. The generic weighted-shuffling calculation agrees with the occupied-edge Gibbs law to roundoff. No Monte Carlo uncertainty is involved in these finite enumerations.

== Corrected orientation adapter
The paper recurrence and this implementation's tiling array have opposite first-axis orientations. Reverse the rows of both input family arrays, apply the existing Gamma reduction and reverse the rows of each resulting creation table. Apply the same conversion to packed creation-coin tables for the efficient paired route. This gives the intended edge law for the new comparison. Corrected probability tables agree with the generic weighted reduction through orders one to eight in the regression tests.

The historical functions and archived results are preserved. The corrected bridge is a distinct path with explicit provenance. The exact mismatch means the old quenched Gamma samples cannot simply be relabelled as correct. This audit does not predict the corrected large-size scaling; a fresh, separately specified assessment is required before using Aztec as the positive benchmark.
#source("aztec/scripts/comparison_geometry.jl; audit_comparison_mapping.jl; docs/validation/aztec_mapping_audit.json")

== Physical probes: an additional correction to the proposed comparison
The Julia table location `(R,C)` has physical coordinates, in units of one primal edge,
$ x=C-1-floor((R-1)/2), quad y=n+1-C-ceil((R-1)/2). $
Thus a central-table-row displacement of $r$ columns is a lattice displacement $(r,-r)$, of length $sqrt(2) r$. It is diagonal in the ordinary square lattice, even though it appears horizontal in the rotated Aztec display. The new square primary observable therefore uses a symmetric diagonal probe, with row/column displacement $(r,r)$. Symmetric horizontal probes are supplementary.

A separate occupied-edge dual-path calculation on every enumerated order-two and order-three tiling confirms unit fluctuation scale. The overall height sign depends on the chosen colour/order convention; covariance is invariant under a common sign reversal. The bridge uses even orders. The new tables store endpoints, direction, actual separation and height units explicitly.

== Prospective square weight placement
With black vertices at odd row-plus-column parity and rows increasing downward, assign independent $a tilde Gamma(0.2,1)$ to edges north from black vertices and $b tilde Gamma(0.25,1)$ to edges east from black vertices. All remaining edges have weight one. This is the paper's red/yellow assignment after undoing its drawing rotation, with a fixed translation choosing the square colour phase.

The 2-by-2 heat-bath probability is $b/(a+b)$ for horizontal occupation. Exact square enumeration checks stationarity and detailed balance for sides two and four under random, deterministic-mean and all-one versions of this placement. Separate parameter-aware seed namespaces leave the September campaign unchanged.

== Joint control audit outcome
The frozen September tables pass deterministic seed reconstruction, global stream-seed uniqueness and centre/spatial join checks. A joint test of the 30 all-one covariance observables uses 9,999 permutations, moving entire second-replica vectors within each size. This preserves within-replica spatial dependence.

For centred replica vectors $u,v$ of length $N$, the sample covariance under a random permutation has null variance
$ frac((sum_i u_i^2)(sum_i v_i^2),(N-1)^3). $
The test uses the largest absolute standardised covariance over all sizes/probes. The observed maximum is 2.26643 and the plus-one randomisation p-value is 0.4827. This audit does not find joint evidence against zero control covariance. It does not certify the single-replica sampler or prove absence of every possible defect.
#source("square_glauber_python/src/square_glauber/control_audit.py; docs/validation/september_control_audit.json")

== Local feasibility and state recovery
A predetermined small pilot requested four environments at each of sides 4, 8 and 12 for each of the three models, with a 4096-sweep cap: 36 environments and 72 streams. All 24 deterministic-control environment pairs certified. The random model certified only one pair at side four and no pairs at sides eight or twelve; 21 random-model streams are censored and retain their original histories. This is a computational feasibility result, not evidence about curvature. No scaling fit is performed on the completed random subset.

The corrected Aztec bridge generated twenty paired environments, four at each of orders 8, 16, 24, 32 and 48. It saves weights and both tilings. This verifies a reproducible workflow, not a positive large-size result.

A search of local historical states found 200 complete pairs that match the September snapshot's hashes, seeds and centre heights. Symmetric probes have been remeasured from those states without new sampling. Another 23,599 snapshot pairs have no locally available verified pair in that search. The available subset is explicitly exploratory; its missing-state inventory is retained for Hamilton recovery.

== Current handoff and remaining decisions
The handbook, corrected model/probe implementation, control audit, remeasurement tool, small pilots and prospective direct coefficient comparison are now the local preparation stage. A prepared Hamilton feasibility design uses sides 8, 12, 16, 24 and 32, sixteen environments per model and a 65,536-sweep cap. It is not a production scaling design. On 30 September, Hamilton validation job 18852141 passed all 108 Python tests, the Julia suite and the independent exact distribution audit in 59 seconds. The uploaded source snapshot was verified against its recorded file hashes. Capped feasibility job 18852143 and saved-state remeasurement job 18852144 then completed successfully, in 15 minutes 57 seconds and 50 minutes 58 seconds respectively. Both used one CPU and a 4 GB memory allocation. The dated job records and logs are retained in `docs/validation/hamilton_20260930/`.

Recover the original state observables, review certification tails and decide whether to extend identical capped histories. The historical campaign and its unresolved pair remain untouched. The September centre and spatial input hashes agree between Hamilton and the archived local snapshot. Before any large campaign, freeze its model, corrected sampler version, physical windows, probe coverage, sample counts, uncertainty plan and stopping rules. Reassessing the corrected Aztec scaling now takes priority over interpreting a boundary-effect contrast against the historical positive fit.

The comparison analysis requires complete authenticated tables, matching laws/directions/height scales and declared fit windows. It refuses incomplete probe vectors and rank-deficient designs; it independently resamples environment blocks in each domain and reports an interval for $c_(Aztec)-c_(square)$. No such scientific contrast is reported from the present tiny pilot.

The follow-up chapter below and `docs/AI_CONTRIBUTIONS.md` record subsequent analysis and confirmed AI assistance, with historical attribution still awaiting Alberto’s account. For commands and the full handoff, consult `docs/history/MATCHED_COMPARISON.md`.


#pagebreak()
= Completed follow-up assessments: 30 September 2026
#note[*Production is on hold.* All five requested areas have been assessed or documented. The evidence does not yet support launching a matched production comparison. Final historical AI attribution and human verification remain to be confirmed by the author.]

== Hamilton outcomes and CFTP diagnosis
Both bounded Hamilton jobs completed successfully. All 160 deterministic-control pairs certified, but none of the 80 random-weight pairs did so within 65,536 sweeps. Only two of 160 random streams certified. The original 23,799 completed September pairs were all recovered and remeasured; the unresolved 23,800th planned pair remains unresolved.

The follow-up regenerated all 240 pilot environments from recorded seeds. Across 30,160 faces in the random environments, 2,809 have minority heat-bath probabilities below $1/65536$. The smallest minority probability is approximately $3.86 times 10^(-33)$. Nine horizontal probabilities, in six environments, round to exactly one in the current floating-point implementation. The mathematical perfect-sampling theorem does not eliminate numerical transition error: coalescence certifies the implemented maps, subject to the correctness and arithmetic assumptions of that implementation.

The first four side-eight environments were selected by index, with both original replica streams. None of the eight streams coalesced by 262,144 sweeps. Median disagreeing-face counts fell from 49 at horizon one to 25.5 at 65,536 and 16.5 at 262,144. Final maximum height gaps are four or eight. Hundreds of thousands of actual flips still occurred. These selected environments have no probabilities rounded to one, so the numerical issue alone does not explain their slow coupling.

These findings are consistent with disorder-related bottlenecks, but do not identify a unique cause or measure the equilibrium mixing time. A coupling may coalesce slowly even when its individual chains relax more quickly. Increasing the cap fourfold did not resolve this fixed cohort. Original checkpoints were not altered; diagnostic endpoints are not certified equilibrium samples.
#source("docs/validation/followup_20260930/weight_diagnostics.csv; cftp_replay.csv; cftp_endpoint_maps.csv")

== Recovered September uncertainty analysis
The exploratory analysis resamples whole environments 10,000 times within size and model. Both replicas and all centre, original rightward and recovered symmetric observables move together. Intervals below are pointwise percentile intervals; they do not include uncertainty from the unresolved environment or provide simultaneous coverage over the collection of fits.

The historical model has independent mean-one Gamma(1) weights on every edge, and remains distinct from the new matched two-family model. At the symmetric diagonal probe with requested separation fraction one quarter, the fitted disorder covariance curvature changes substantially with the size cutoff:
#table(columns: (1fr, 1fr, 2fr), inset: 5pt, stroke: 0.4pt,
 [Square sizes], [Squared-log coefficient], [Pointwise 95% interval],
 [8–32], [0.9955], [0.4093, 1.6003],
 [12–32], [-0.3801], [-1.6387, 0.9152],
 [16–32], [-3.2509], [-6.0980, -0.3325],
)
The positive full-window value therefore does not establish robust asymptotic positive curvature. The centre coefficient is 0.0446 with interval [-0.3311, 0.4218], reproducing the earlier point estimate. The paired symmetric-horizontal minus original-rightward curvature contrast over sizes eight through 32 is 0.6516 with interval [-0.0735, 1.3774]. All fitted uniform-control curvature intervals include zero.

There are 152 moment rows, 24 accepted exploratory fits and 18 refused designs. The planned pooled small-probe fits lack adequate coverage or distinct physical separations after rounding; absent probes were not filled. Supplementary single-quarter-probe and centre fits were specified in the analysis code before execution, though not separately listed in the initial plan JSON. No fit here is a matched-law Aztec–square contrast.
#source("docs/validation/followup_20260930/recovered_analysis.json; recovered_intervals.csv; analysis_amendments.json")

== Corrected Aztec reassessment
The corrected orientation adapter generated 384 paired environments: 64 at each of orders 16, 32, 48, 64, 96 and 128, with fresh master seed 2026093002. Weights and both tilings are retained. This is a bounded pilot using the mathematical weighted-shuffling algorithm in floating point; the tiny exact-distribution checks validate their enumerated cases, not arbitrary-precision sampling at every size.

A coverage audit withdrew the initial three-probe full-window fit because rounding removes the one-sixteenth probe at order 16. An explicitly exploratory two-probe subset was added after initial results were inspected. The amendment is recorded rather than presented as prospective confirmation.
#table(columns: (1.5fr, 0.9fr, 0.9fr, 1.8fr), inset: 5pt, stroke: 0.4pt,
 [Probe fractions], [Orders], [Coefficient], [Pointwise 95% interval],
 [1/16, 1/8, 1/4], [32–128], [-0.4056], [-3.3970, 2.6176],
 [1/16, 1/8, 1/4], [48–128], [0.2222], [-7.0879, 7.1973],
 [1/8, 1/4], [16–128], [0.6586], [-0.7506, 2.0800],
 [1/8, 1/4], [32–128], [1.2169], [-2.5525, 5.0765],
 [1/8, 1/4], [48–128], [5.0796], [-4.1941, 14.0347],
)
All intervals include zero. The corrected pilot neither confirms nor excludes a positive squared-log term. The historical finite-environment orientation mismatch does not, by itself, determine the effect on disorder-averaged scaling; no equivalence theorem for the old and corrected annealed observables has been established here.
#source("docs/validation/followup_20260930/aztec_analysis.json; aztec/output/corrected_reassessment_20260930/")

#pagebreak()
== Production decision
Do not launch production with the current setup. First validate numerical rare-event handling without silently clipping weights or changing the law; then establish completion costs on fixed retained square environments. Any revised random map needs a versioned stream/checkpoint identity. The corrected Aztec benchmark also needs substantially greater precision, complete physical-probe coverage and fixed cutoff sensitivities.

For planning only, the two-probe Aztec full-window interval has half-width about 1.4153 with 64 environments per size. A naive $1/sqrt(N)$ extrapolation to half-width 0.25 suggests approximately 2,052 environments per size in the same design. This is not an approved allocation, power calculation or guarantee: pilot variance is uncertain, the precision target needs scientific justification, and the direct contrast also requires square uncertainty. Square production sizes and counts cannot yet be fixed responsibly from zero completed random pairs.
#source("docs/history/PRODUCTION_DECISION.md; docs/history/FOLLOWUP_20260930.md")

== AI assistance and author review
Codex assisted with the handbook, model and coordinate mapping, audit implementation, corrected adapter, tests, bounded HPC execution, state recovery, statistical analysis and interpretation. These contributions are recorded in `docs/AI_CONTRIBUTIONS.md`. Alberto supplied the correspondence, research objective, instructions, access and authorisation. Earlier personal mathematical work and historical AI involvement have not been inferred from repository author fields.

Automated tests and independent computational formulations are evidence about tested properties, not independent human review; many tests were themselves AI-assisted. The disclosure draft records confirmed assistance and marks historical attribution and human mathematical verification as unconfirmed. It must not be converted into a claim that every result has been personally checked until that review occurs. No disclosure or email has been sent.


#pagebreak()
= Numerical remediation: 1 October 2026
#note[*Original-input small-width sampling is verified.* The direct integer sampler completed the fixed side-eight and side-twelve cohort using the original Hamilton archives. Large-width feasibility remains open. Hamilton validation job 18872074 passed. The corrected Aztec refinement completed as array 18897605; its results are in the next chapter.]

== Exactness for stored edge weights
Every finite positive Float64 weight is a dyadic rational. Multiplying all edges by a common power of two makes them integers. Since every matching has the same number of edges, this multiplies every matching weight by the same constant and preserves its Gibbs probability. A common integer divisor may also be removed.

For integer weights, generate an unbiased integer below a positive bound by drawing enough random bits and rejecting out-of-range values. Comparing that integer with an exact numerator implements a rational Bernoulli coin without floating-point division. Very rare outcomes remain possible. This statement assumes ideal independent bits; the implementation uses a reproducible pseudorandom bit stream, and exactness is conditional on the stored weights, not an exact continuous Gamma generator.

The new reference heat-bath map uses the rational face probability $a c/(a c+b d)$, a uniformly selected face, and one coin shared by the lower and upper configurations even on null faces. It preserves the same local orientation maps. It is not silently substituted into the old CFTP checkpoints; the stream and algorithm versions differ.

== Direct frontier partition sampling
Scan the square's vertices in row-major order. A bit mask records whether each vertex at the scanning boundary has already been matched. At the next unmatched vertex, any new partner must lie to its right or below. Let $Z(s)$ be the total integer weight of completions from a position-and-mask state $s$. Then
$ Z(s)=sum_(e " allowed") W_e Z(s_e). $
Occupied vertices are skipped deterministically; the final empty mask has value one and other final masks value zero. Cache these values and select edge $e$ with probability $W_e Z(s_e)/Z(s)$. The successive conditional factors telescope to the matching weight divided by the total partition function. This samples the Gibbs law directly, without mixing or coalescence.

The implementation deliberately limits the width to twelve and the number of cached states to two million. Its cost grows exponentially with width. It is currently an independent small-grid oracle and practical small-width sampler, not an established large-grid replacement for CFTP.
#table(columns: (1.7fr, 0.6fr, 0.8fr, 1.2fr, 1.2fr), inset: 5pt, stroke: 0.4pt,
 [Model], [Side], [Pairs], [Total setup time], [States per environment],
 [Matched random Gamma], [8], [16/16], [0.0495 s], [3,361],
 [Matched random Gamma], [12], [4/4], [0.3572 s], [104,284],
 [Two control models], [4,8,12], [6/6], [See record], [At most 104,284],
)
These are local timings for the original-input rerun, excluding interpreter startup. All 240 original Hamilton archives were downloaded and verified against the pilot hashes. Comparing their arrays with local seeded regeneration found 176 identical environments and 64 with tiny numerical differences, all in the random Gamma arm (at most 9 ULP; maximum relative difference $1.225 times 10^(-15)$). The precise cause is not isolated; the discrepancy is not solely archive packaging. Earlier regenerated diagnostics retain this qualification.

The rerun loaded original files directly: all 24 overlapping input archives match the pilot hashes; the two side-four controls are new. All 52 saved tilings satisfy the matching constraints, and all 224 saved observable values agree with recomputation. Replica seeds use the new `integer-frontier-v1` namespace. These samples do not resume or replace old CFTP streams. No scaling fit was drawn from this small cohort.

Validation compares every tiny-grid matching mass with independent rational Gibbs products, verifies exact detailed balance for every flippable four-by-four matching, tests scripted rare outcomes and integer rejection boundaries, and checks common-map height order and deterministic replay.
#source("square_glauber_python/src/square_glauber/integer_sampling.py; docs/validation/numerical_20261001/integer_original_feasibility_summary.json")

== Robust Aztec creation coins
The new `gamma-interval-coins-v1` algorithm uses outward neighbouring Float64 values to bound $b/(a+b)$. A 53-bit uniform prefix gives an interval containing the underlying uniform number. If that interval lies entirely above or below the probability interval, the coin is settled. Otherwise, convert the represented weights to exact rationals and extend the same random prefix in 64-bit blocks until integer comparisons settle the coin. Do not redraw the entire uniform after detecting a near-threshold case, since that would condition the draw incorrectly.

This fixes rounding in creation coins for the represented reduced weights. The reduction still uses Float64. Three fixed environments were separately reduced with 256- and 512-bit arithmetic. At order 128, the maximum creation-probability discrepancy of the Float64 reduction relative to the 512-bit calculation was $5.52 times 10^(-15)$; twice the sum over all potential coins was $1.30 times 10^(-10)$. The old coin formula rounded 277 probabilities to endpoints in that environment. The 256-versus-512-bit discrepancy was below $4.1 times 10^(-76)$.

These are measured discrepancies for the audited environments, not a universal interval proof. The exact coin tests include rare tails, subnormal inputs and sums overflowing Float64. Corrected orientation is preserved; the historical functions remain unchanged for reproducibility.
#source("aztec/scripts/robust_gamma_coins.jl; audit_robust_gamma.jl; docs/validation/numerical_20261001/aztec_numeric_audit.json")

== Prepared precision refinement and supervisor update
The next corrected Aztec design fixes 2,048 paired environments at each order 16, 32, 48, 64, 96 and 128. It uses only the one-eighth and one-quarter probes, which have full coverage at every size, with windows 16–128, 32–128 and 48–128. Master seed is 2026100102; the planned 10,000-resample bootstrap seed is 2026100103. This improves precision but does not guarantee a positive signal or the target interval half-width.

The 24 integration pairs passed the generation-to-analysis pipeline; they are not refinement results. The Hamilton array allows two concurrent tasks, each one CPU, 4 GB and two hours. Analysis checks counts, hashes, source versions, seed uniqueness and geometry. All 135 source hashes were verified on Hamilton. Validation job 18872074 completed in one minute: all 120 Python tests, 21,970 Julia assertions, the numerical audit and benchmark passed. Refinement array 18897605 was submitted on 2 October and completed. The next chapter records the frozen analysis, checked on 3 October.

The large matched comparison remains on hold pending wider-square feasibility. The unsent supervisor draft and AI contribution ledger record the audit, small-grid results and proposed refinement.
#source("docs/history/NUMERICAL_REMEDIATION_20261001.md; aztec/hpc/robust_refinement.slurm; docs/AI_CONTRIBUTIONS.md")


#pagebreak()
= Corrected Aztec results: 3 October 2026
#note[*All 12,288 paired environments completed.* The primary squared-log coefficient is $c=0.3568$, with a 95% bootstrap interval $[0.0723,0.6327]$. This supports positive finite-window curvature; it does not establish an asymptotic law or a boundary contrast with squares.]

Hamilton array 18897605 completed without errors or missing cases: orders 16, 32, 48, 64, 96 and 128 each have 2,048 paired environments. Model, orientation, sampler and probes were fixed beforehand. Raw weights and tilings remain on Hamilton; compact tables and provenance were downloaded.

The fit is $D_p(d)=a_p+b_p log d+c(log d)^2$: each probe has its own intercept and log slope, with a common quadratic coefficient. It uses ordinary least squares on the covariance estimates. The 10,000 bootstrap resamples move complete environment vectors, including both replicas and probes, with independent resampling across sizes. Seed is 2026100103. Intervals below are pointwise and correlated across windows.

#table(columns: (1.2fr, 1fr, 1.6fr), inset: 5pt, stroke: 0.4pt,
 [Aztec order window], [Coefficient $c$], [95% bootstrap interval],
 [16–128 (primary)], [0.3568], [0.0723 to 0.6327],
 [32–128], [0.7362], [0.0049 to 1.4770],
 [48–128], [1.9256], [0.1376 to 3.7377],
)

#image("../figures/refinement_20261003.png", width: 100%)

All intervals exclude zero, but estimates drift upward and uncertainty widens as smaller orders are removed. The middle interval only narrowly excludes zero. These correlated checks do not establish a stable asymptotic coefficient. The primary half-width is 0.2802, above the planning target 0.25; no sample-count extension was triggered.

Verification covered archive and table hashes, full probe coverage, all seed derivations, source identity and independent recomputation of all 36 moment estimates. Tilings were validated during generation; this follow-up did not revalidate all raw tilings. Float64 reduction remains an approximation, with the earlier three-environment high-precision audit providing numerical evidence rather than a global proof.

*Next scientific step:* resolve wider-square sampling under the same weight law, then freeze physically matched comparison windows and estimate the direct domain contrast. The small-width integer sampler is limited to $L<=12$; the capped matched CFTP pilot cannot establish absence of square curvature. The supervisor draft has been updated but not sent. Human review and historical AI attribution remain outstanding.
#source("docs/history/AZTEC_REFINEMENT_20261003.md; docs/validation/refinement_20261003/analysis.json; verification.json")


#pagebreak()
= Wider-square feasibility: 4 October 2026
#note[*All 18 retained environments now have verified pairs.* The first stage completed sides 14 and 16. A separately recorded follow-up completed the same six side-18 inputs after they reached the original state cap. This establishes cohort feasibility, not square curvature.]

The new `integer-frontier-budgeted-v1` sampler retains the exact conditional-partition law for stored dyadic weights. At each black vertex $b$, scale incident edges by a rational factor $g_b$ to obtain smaller integers. Every matching uses one edge at each black vertex, so its weight acquires the same factor $product_b g_b$. Gibbs probabilities are unchanged; the input Gamma law is preserved.

Already occupied vertices are skipped, remaining states use a packed integer key, and an explicit stack replaces recursion. Sampling uses exact integer conditional probabilities. The new representation changes random-bit consumption, so it has a separate stream namespace. The old frontier and CFTP code is unchanged. Width dependence remains exponential, with a hard implementation limit of side 18.

A local side-12 comparison confirmed exact partition equality after gauge correction: cached states fell from 104,284 to 52,341 and partition bit length from 6,430 to 4,311. All 132 Python tests passed locally and on Hamilton (validation job 18902368). Tests cover rational matching masses, conditional probabilities, gauge invariance, rare weights, budgets and output verification.

#table(columns: (0.5fr, 0.7fr, 1.2fr, 1.2fr, 0.9fr), inset: 5pt, stroke: 0.4pt,
 [Side], [Final pairs], [Cached states], [Gamma worker time], [Max. RSS],
 [14], [6/6], [268,759], [1.47–1.73 s], [0.236 GiB],
 [16], [6/6], [1,333,766], [7.53–7.66 s], [0.798 GiB],
 [18], [6/6], [6,448,885], [39.31–39.75 s], [4.356 GiB],
)

Each size has four random Gamma pairs, one deterministic-mean control pair and one uniform pair. Side 16 uses the first four retained random Hamilton environments and original controls, verified by archive hashes. Sides 14 and 18 use new saved environments with master seed 2026100301.

First-stage array 18902369 returned twelve complete pairs and six side-18 outcomes capped at three million states. These outcomes remain preserved. The capped cases used at most 1.247 GiB and about 17 seconds, so they had not tested the existing 120-second and 6 GiB soft limits. After explicit resource review, a separate six-case plan and array 18902387 reused the exact same inputs and seeds, raising only the state cap to eight million. All six completed. Both stages allowed two concurrent workers, with four-minute, one-CPU, 8 GiB scheduler ceilings per worker.

All scheduler tasks completed with exit code zero and empty error logs. Full input and sample archives were downloaded. Verification checked hashes and original identities, all 36 sampled perfect matchings and all 72 saved observable values. No failed environment was replaced. The six capped first-stage outcomes are not relabelled as completed; the follow-up is a separate record.

*What remains:* cached-state counts grew by roughly fivefold per increase of two in side length. This small cohort does not establish large-ensemble runtime tails or justify a broad square scaling campaign. No curvature fit is made: the size window is narrow, random-pair counts are small and rounded probe distances repeat. A scalable exact square sampler must be validated against these integer oracles before fixing a large matched-domain comparison. Generalized shuffling with forbidden edges is one candidate; it has not yet been implemented. No further campaign was submitted.
#source("docs/history/SQUARE_FEASIBILITY_20261004.md; evidence directory: docs/validation/frontier_scaling_20261003/")

#pagebreak()
= Matched shuffling comparison: 4–5 October 2026
#note[*The scalable square sampler has passed local validation.* Hamilton validation passed and the six-size comparison is submitted as array 18907659, with dependent analysis job 18907660. The current evidence is validation, not a new square scaling result. This chapter supersedes the preceding chapter's statement that embedded shuffling is unimplemented.]

== Preserve the square law inside a sampling host
Let the target square have $L times L$ vertices. Embed it in an Aztec diamond of order $L-1$, preserve every target edge weight, set all crossing edges to zero and exterior edges to one. Check an explicit perfect matching of the exterior. Every allowed host matching splits into independent square and exterior matchings, giving $Z_("host")=Z_("square") Z_("exterior")$. Projecting a host Gibbs sample onto the square therefore preserves its original Gibbs law and boundary conditions. The diamond is an algorithmic device, not a replacement physical domain.

Weighted shuffling reduces an active cell with clockwise edge weights $a,b,c,d$ using $a c+b d$ and opposite-edge weights divided by this factor. On growth, move a lone old edge to its opposite, remove an old opposite pair, and fill an old empty cell by a horizontal or vertical pair with probabilities $a c/(a c+b d)$ and $b d/(a c+b d)$. Creation cells are selected from the old matching; refilling every newly uncovered cell after deletion would be wrong.

This implementation supports zero edges only while each reduction cell has positive double product. The positive-weight distribution identity then extends continuously to the forbidden crossing edges. Singular zero cells, positive-product underflow and nonfinite arithmetic cause an error; no epsilon is introduced. The general singular-cell algorithm is not implemented. The construction follows #link("https://www.unige.ch/math/folks/velenik/papers/JRV06err.pdf")[Janvresse, de la Rue and Velenik, “A Note on Domino Shuffling”, section 2].

== Exact reference and numerical approximation
Rational arithmetic gives an exact small-grid reference for represented input weights, conditional on ideal random bits. Every matching probability agrees with an independent Gibbs enumeration at sides two and four for the random law and both controls. An arbitrary all-edge-weight side-four case also agrees with the integer frontier oracle. The full uniform side-six distribution contains 6,728 matchings with probability $1/6728$ each.

The scalable sampler uses Float64 reductions and exact rational creation coins for the represented reduced weights. Those coins preserve rare outcomes, but do not remove reduction rounding. The fixed audit compares every potential creation probability with 80- and 160-digit reductions at sides 8, 16, 32 and 64. At side 64, the maximum discrepancy is $1.79 times 10^(-15)$ and twice the sum across potential coins is $7.94 times 10^(-12)$. These are finite-cohort measurements, not a universal precision proof. Arithmetic-operation count is cubic in side length; exact rational bit complexity is not bounded by that statement.

All 146 Python tests pass locally. A separate local pipeline check completed four paired environments in every one of 36 domain/model/size cells. These 144 pairs are validation-only and excluded from inference. Saved weights, matchings and observables have dedicated revalidation tools. Hamilton validation subsequently passed before production submission.
#source("square_glauber_python/src/square_glauber/embedded_shuffle.py; docs/validation/embedded_shuffle_20261004/numeric_audit.json")

== Fixed comparison, without a predetermined result
Both geometries use independent north-from-black $a$ weights with Gamma shape 0.2, east-from-black $b$ weights with Gamma shape 0.25, scale one for both, and remaining weights one. These families have means 0.2 and 0.25: normalising them to mean one would change the experiment. The ordinary square's two-family law is now the same law as the Aztec arm. Earlier all-edge Gamma(shape one) results remain separate.

#table(columns: (1.1fr, 2.2fr), inset: 5pt, stroke: 0.4pt,
 [Choice], [Prespecified design],
 [Domains and sizes], [Square and fresh corrected Aztec: 8, 16, 24, 32, 48, 64.],
 [Random data], [2,048 paired environments per size and domain.],
 [Controls], [Deterministic family means and all-one; 256 pairs per size, domain and model.],
 [Total], [30,720 paired environments, including controls.],
 [Probes], [Symmetric diagonal increments at fractions 1/8 and 1/4; common physical distances and height scale one.],
 [Primary window], [8–64; sensitivity windows 16–64 and 24–64.],
 [Resampling], [10,000 whole-environment bootstrap draws, independent between domains, models and sizes.],
 [Resource ceiling], [36 tasks; two concurrent; one CPU, 4 GiB and two hours per task.],
)

Fit $D_p(d)=a_p+b_p log d+c(log d)^2$, allowing probe-specific intercepts and log slopes. The primary comparison is the direct difference $Delta c=c_(Aztec)-c_(square)$, with its own bootstrap interval. Separate significance decisions do not test that contrast. A joint 9,999-permutation control audit moves each complete second-replica probe vector within its size/model/domain group.

A square interval containing zero does not demonstrate absence of a squared-log contribution. A positive domain contrast would support a difference over the sampled window, not establish an asymptotic boundary theorem. Report every fixed cutoff window and any control discrepancy. Missing environments or probes block fitting; difficult environments are never silently replaced.

The production master seed is 2026100407; validation uses 2026100408. Source, input, sample and table checksums are retained. After renewed SSH access, all 155 source hashes were verified on Hamilton. Validation job 18907657 completed in 2 minutes 56 seconds with exit code zero: all 146 Python tests, 21,970 Julia assertions, both numerical audits and 144 validation pairs passed. Reloading the raw files verified 288 matchings and 576 observable values. The validation time request was reduced to the test partition’s 15-minute limit; scientific code and production budgets were unchanged.

Production array 18907659 is submitted. Verification and analysis job 18907660 depends on successful completion of every sampling task. It checks all raw records before the fixed bootstrap analysis; failures block automatic analysis. The side-64 Gamma validation cell took 5.14 seconds for four pairs and used about 348 MB, within the planned production limits. Results remain pending. Submitted jobs continue independently of the SSH connection. The AI ledger records confirmed assistance; historical attribution and independent human mathematical review remain outstanding.
#source("docs/history/MATCHED_SHUFFLE_20261004.md; square_glauber_python/configs/matched_shuffle_20261004.json; docs/AI_CONTRIBUTIONS.md")

#pagebreak()
= Completed matched comparison: 5 October 2026
#note[*The requested matched-law experiment is complete; the proposed scaling contrast is not established.* All 30,720 paired environments, raw verification and the fixed analysis completed. The primary difference is $Delta c=-0.046$, with 95% interval $[-0.366,0.282]$. A control-curvature discrepancy also requires follow-up.]

Hamilton array 18907659 completed all 36 sampling tasks, and job 18907660 completed raw verification and analysis. All exits were zero and error logs were empty. No environment was omitted or replaced. The data comprise 2,048 random pairs and 256 pairs for each of two deterministic controls at each of six sizes in each geometry. The slowest sampling cell took 44 minutes 17 seconds; post-processing took 22 minutes 9 seconds.

== Primary result and fixed cutoff checks
Use the prespecified fit $D_p(d)=a_p+b_p log d+c(log d)^2$, with probe-specific intercepts and slopes. The 10,000 bootstrap draws resample complete paired-environment vectors; intervals are pointwise. The direct comparison is $Delta c=c_(Aztec)-c_(square)$.

#table(columns: (0.9fr, 1.35fr, 1.35fr, 1.45fr), inset: 4pt, stroke: 0.4pt,
 [Window], [Aztec $c$ [95% CI]], [Square $c$ [95% CI]], [Difference [95% CI]],
 [8–64, primary], [0.238\ [0.010, 0.467]], [0.284\ [0.050, 0.517]], [−0.046\ [−0.366, 0.282]],
 [16–64], [0.227\ [−0.387, 0.837]], [−0.197\ [−0.806, 0.413]], [0.424\ [−0.434, 1.294]],
 [24–64], [−1.391\ [−2.847, 0.045]], [0.532\ [−0.909, 1.978]], [−1.922\ [−4.022, 0.129]],
)

Both primary individual intervals exclude zero; the interval for their difference does not. None of the three fixed contrast intervals excludes zero. Coefficients change substantially when small sizes are removed. This neither proves different asymptotic laws nor proves equal laws. The earlier corrected Aztec result over 16–128 remains a separate, wider-window result and cannot replace the matched 8–64 comparison.

#image("../figures/matched_shuffle_20261005.png", width: 95%)

#pagebreak()
== The control discrepancy
For a deterministic environment, the true disorder covariance is zero. The prespecified joint permutation test of the individual control covariance coordinates gives $p=0.5882$. Nevertheless, the all-one primary curvature contrast is $1.161$ with interval $[0.451,1.853]$: Aztec $c=0.642$ and square $c=-0.519$. This should not be dismissed because a different test gives a large p-value.

A separate *post hoc* diagnostic permutes the complete second-replica probe vector within each control/domain/size group. It calculates the same curvature contrasts, divides each by its exact permutation-null standard deviation and takes the maximum absolute statistic across the six contrasts (two control models, three windows). With 9,999 permutations and seed 2026100514, the maximum-statistic p-value is 0.0093. The uniform primary contrast has pointwise permutation p=0.0021 and six-contrast adjusted p=0.0093.

The statistic was chosen after the anomaly was observed. Adjustment across its six contrasts does not account for the broader choice of diagnostic; this is exploratory evidence requiring investigation, not a new prespecified test or proof of a sampler bug. It does not alter the frozen Gamma analysis. Maximum individual covariance and curvature contrast are different statistics, so their p-values are not logically contradictory. Chance, dependence or provenance errors, and numerical or analysis problems must be distinguished by focused follow-up.

== Verification, records and next steps
Hamilton verification reloaded 61,440 matchings, recomputed all 122,880 observable values and checked 138,240 files. The compact archive was downloaded with checksum verification; full raw files remain on Hamilton and were not all downloaded a second time. Local checks confirm all counts, seed identities, uniqueness, physical probe conventions and table hashes. Independent direct calculations reproduce all 216 moment estimates and 18 fitted coefficients, with maximum differences below $1.8 times 10^(-15)$ and $2.3 times 10^(-14)$ respectively. The original bootstrap was retained.

The requested comparison has now been performed under the same Gamma law and matched observables. The correct next step is to report this result to Sunil, including uncertainty and the control anomaly; the email remains an unsent draft. Investigate the control discrepancy and review physical probe positions, parity and finite-size sensitivity before fixing any broader campaign. Do not extend the sample until an expected pattern appears. Historical AI attribution and independent human mathematical review remain outstanding. No additional sampling campaign has been launched.
#source("docs/history/MATCHED_RESULTS_20261005.md; docs/validation/matched_shuffle_20261005_results/; square_glauber_python/outputs/matched_shuffle_20261005/production/analysis/")

#pagebreak()
= Control investigation and confirmation: 6–7 October 2026
#note[*The independent confirmation is complete.* All 12,288 fresh pairs passed verification. The primary contrast is $0.081$ with interval $[-0.282,0.437]$ and permutation $p=0.6615$. The original anomaly was not reproduced at this precision; this does not prove sampler correctness. All 38 selected original pairs replayed exactly, and the completed audit identified no pairing or geometry defect.]

== What deterministic controls test
With deterministic weights and independent random inputs $U_1,U_2$, let $X_1=f(U_1)$ and $X_2=g(U_2)$. Independence gives $EE[X_1 X_2]=EE[X_1] EE[X_2]$, hence $Cov(X_1,X_2)=0$. This remains true if either sampler has an imperfect marginal distribution. A deterministic probe offset, common height-scale error or marginal sampling bias alone cannot produce nonzero expected cross-replica covariance. Finite-sample chance, dependence, incorrect record association and analysis errors must be considered. Conversely, a zero covariance control would not by itself prove the marginal Gibbs law.

The original uniform primary contrast remains 1.160615. Its original interval and both the prespecified and exploratory p-values remain unchanged. New diagnostics do not replace that result.

== Saved-data and source audit
Both physical probe endpoints, expressed relative to the domain centre, agree exactly in the two geometries at all six sizes. The one-eighth probes at sizes 8 and 24 have odd separation and common midpoint $(-0.5,+0.5)$; the other selected probes have midpoint $(0,0)$. Thus parity is a documented finite-size feature, not an endpoint mismatch. It may affect Gamma fits, but cannot explain a nonzero expected deterministic-control covariance under independent replicas.

All 108,000 role seeds across old and new scientific/validation cohorts are distinct. This is a provenance check, not a proof about the pseudorandom generators. The source uses separate seeded local random streams for square replicas and distinct consecutive random draws with separate choice arrays for Aztec replicas. Under ideal independent draws, these constructions give independent replicas.

The original all-one contrast receives positive realised contributions from each size: 0.171, 0.162, 0.296, 0.040, 0.027 and 0.464 at sizes 8, 16, 24, 32, 48 and 64 respectively. Removing each whole pair in turn, solely as an influence diagnostic, gives contrasts between 1.111913 and 1.194972. The largest individual change is 0.048702. No observation was actually removed from the experiment.

An independent finite check of the exact permutation covariance formula exhausted all 24 pairings of a four-pair, two-probe example and agreed to $1.78 times 10^(-15)$, including cross-probe terms. The existing anomaly is therefore not explained by a detected formula error or a single extreme pair. These checks do not identify its cause.

#pagebreak()
== Exact replay and a fixed confirmation
Hamilton job 18919384 replayed pairs 1, 128 and 256 at each size in both geometries, plus the most influential size-64 square and Aztec pairs: 38 pairs total. Both saved matching arrays agreed exactly. Square replicas were regenerated in reverse order, with unchanged reduction tables; Aztec choice arrays were distinct. The job also passed the new 48-pair smoke pipeline, raw verification of 96 matchings and recomputation of 192 values. It completed in 57 seconds with exit zero and empty error logs.

The confirmation uses the unchanged verified sampler source, all-one weights only, and 1,024 fresh pairs per size/domain: 12,288 pairs total. Its sole primary test is a two-sided whole-vector permutation test of the 8–64 Aztec-minus-square curvature against zero, with 9,999 permutations. It retains 10,000 bootstrap draws and descriptive cutoff fits at 16–64 and 24–64. Production master seed is 2026100601; validation is 2026100602; bootstrap seeds are 2026100603 and 2026100604; permutation seed is 2026100605.

The plan was fixed before fresh sampling. Report the result regardless of sign or significance, keep it separate from the original data and do not extend it adaptively. A fresh rejection flags evidence requiring investigation; non-rejection means the earlier anomaly was not reproduced at this precision, not that correctness has been proved. Twelve tasks allow at most two concurrent workers, one CPU, 4 GiB and 45 minutes each. Verification and analysis require success of the whole sampling array. No further Gamma scaling campaign has been launched.
#source("docs/history/CONTROL_INVESTIGATION_20261006.md; docs/validation/control_investigation_20261006/")

#pagebreak()
== Completed confirmation and decision
All twelve tasks in array 18919405 completed successfully. Job 18919413 verified the raw data and analysed the fixed cohort in 24 minutes 31 seconds; all thirteen error logs are empty. The longest sampling task took 21 minutes 36 seconds. The scheduler records completion on 6 October at 14:03:21; the evidence was retrieved and locally verified on 7 October. No pair was dropped, replaced or added.

#table(columns: (0.9fr, 1.35fr, 1.35fr, 1.45fr), inset: 4pt, stroke: 0.4pt,
 [Window], [Aztec $c$ [95% CI]], [Square $c$ [95% CI]], [Difference [95% CI]],
 [8–64, primary], [0.039\ [−0.213, 0.290]], [−0.042\ [−0.297, 0.208]], [0.081\ [−0.282, 0.437]],
 [16–64], [−0.026\ [−0.649, 0.593]], [−0.006\ [−0.592, 0.579]], [−0.020\ [−0.867, 0.829]],
 [24–64], [0.934\ [−0.605, 2.448]], [−1.082\ [−2.557, 0.378]], [2.016\ [−0.108, 4.104]],
)

The sole confirmatory whole-vector permutation test gives $p=0.6615$ with exact conditional null standard deviation 0.182312. The 16–64 and 24–64 fits remain descriptive; their broad intervals are reported without choosing a new primary cutoff. Original and fresh cohorts are separate. Non-reproduction is compatible with sampling variation but does not establish the original cause or prove independent replicas.

Hamilton checked 24,576 matchings, 49,152 observable values and 55,296 files. A checksummed compact archive was downloaded; full raw tilings remain on Hamilton. Local verification rechecked coverage, seeds, the frozen plan and table hashes, then independently reproduced 72 moments, six coefficients and three contrasts. Maximum moment and coefficient differences were below $9 times 10^(-16)$ and $1.5 times 10^(-14)$. Repeating the fixed 9,999-permutation test gave the same p-value. Original bootstrap intervals were retained rather than independently rerun.

*Next scientific step:* review the completed matched-law evidence and mathematical sampler validation with Sunil. The Gamma contrast remains $-0.046$ with interval $[-0.366,0.282]$; the control confirmation cannot turn it into evidence for different growth laws. A wider matched window would require a new fixed design and a precision/runtime justification. No further campaign has been submitted. The supervisor draft is updated but unsent; Alberto's historical AI-attribution account and independent human review remain outstanding.
#source("docs/validation/control_investigation_20261006/confirmation_verification.json; square_glauber_python/outputs/control_confirmation_20261006/confirmation/analysis/analysis.json")

#pagebreak()
= Correspondence and scope: 7 October 2026
#note[*A completed experiment is not a settled hypothesis.* The newly supplied July–September emails clarify the intended comparison and several differences from the experiment now completed. This chapter summarises that correspondence, not a new numerical analysis.]

Sunil first requested fresh-environment LERWs, then shared-environment pairs, Gamma-disordered Aztec heights and a spatial variance decomposition. On 4 August he explicitly endorsed comparing paired spatial increments between domains. His later concern was whether random-weight spanning-tree sampling represented the desired dimer law; this motivated direct Glauber sampling. On 29 September he preferred the original Aztec weights on ordinary squares and endorsed CFTP. Alberto replied that he would reproduce the original setup and continue CFTP.

*What matches:* the new comparison uses the original Gamma shapes 0.2 and 0.25, scale one, other weights one; two replicas share each fixed environment. Physical directions, endpoint separations and height units are matched. A row of the staggered Aztec height table is physically diagonal, which explains the square probe direction. The Gamma family means remain 0.2 and 0.25, not one.

*What changed:* the historical August Aztec experiment reported orders 128–1300 and four fractions, 1/32, 1/16, 1/8 and 1/4. The fresh matched comparison reaches only 64 and uses the last two fractions; the separate corrected Aztec cohort reaches 128. The matched campaign uses shuffling after CFTP proved impractical, with Float64 reduction qualifications. It must not be described as a CFTP production run or full replication of the historical large-size study. Historical positive fits need reassessment because of the orientation finding; that finding alone does not show that their entire annealed signal was artificial.

*Which theoretical observable?* The July email refers to a colleague's double-dimer height-difference result; the August interpretation concerns spatial disorder covariance. For conditionally independent replicas, $Var(X_1-X_2)/2=T$ and $Cov(X_1,X_2)=D$, but those identities alone impose no growth rate on either component. The correspondence does not supply the precise theorem, coupling, centring or limit regime. Obtain the statement from Sunil before identifying it with our fitted $D$, or asserting that $T$ must be logarithmic for this model.

*AI attribution:* Alberto's 2 August email reports ChatGPT help with precision-weighted analysis. His reply to 29 September reports explanations, debugging, language porting and SLURM help. These are historical self-reports, not a complete authorship audit. Later Codex work includes implementation, experiment design, execution, analysis and drafting. Sunil's August preference against AI-assisted manuscript writing and September acceptance of disclosed AI assistance should be reconciled by the human authors before paper drafting; neither email verifies a journal policy.

The next review should cover the sampler correction, change of algorithm, smaller size window, precise theoretical observable and inconclusive matched contrast. A larger fixed experiment may be justified after feasibility and precision planning. The historical claim of no positive square signal cannot establish absence for the new model. No new sampling is launched by this correspondence review.
#source("docs/correspondence/reconciliation_20261007.md; docs/AI_CONTRIBUTIONS.md. Source: email text supplied by Alberto; referenced attachments not independently inspected.")

#pagebreak()
= Larger campaign planning: 7 October 2026
#note[*A larger physical size window is worth investigating.* Alberto requested consideration of a larger campaign without time pressure. The bounded feasibility gate has completed successfully; the full production design is frozen and submitted, as recorded in the next chapter. More samples reduce random error, but do not remove finite-size bias or guarantee a contrast.]

The present square shuffling implementation supports host orders through 127, hence square sides through 128. Cached geometry and reduction storage grow approximately cubically. Reaching 256 requires a deliberate implementation extension and validation, not merely a longer scheduler time limit.

Array 18960990 tests square sides 96 and 128 with eight Gamma pairs, two deterministic-family-mean pairs and two all-one pairs each: 24 pairs total. Two additional numerical tasks compare all possible creation probabilities with 80- and 160-decimal-digit reductions, using the original tolerances. Four tasks have one CPU, 16 GiB and one hour each, with at most two concurrent. Dependent job 18960991 verifies saved matchings and observables. Pilot observations are excluded from scientific production and no curvature fit is made from them.

Both numerical tasks passed. At side 96, maximum probability discrepancy was $3.32 times 10^(-15)$ and twice the sum of discrepancies was $4.02 times 10^(-11)$. At side 128 the corresponding values were $4.91 times 10^(-15)$ and $1.17 times 10^(-10)$. Audit peak memory was 1.35 and 3.09 GiB respectively. These fixed-environment checks do not constitute universal rounding bounds. All four pilot tasks and dependent verification completed with exit code zero and empty error logs. Downloaded raw data independently passed verification: 24 pairs, 48 matchings and 96 observable values. The side-128 Gamma cell took 115.9 seconds for eight pairs and peaked at about 1.73 GiB. The SSH session subsequently expired.

== Illustrative precision and runtime
Using only the old completed Gamma data, estimate the covariance of centred-product influence vectors for the two sample-covariance estimates. Preserve dependence between probes, then combine independent size/domain blocks with the curvature fit weights. Above size 64, extrapolate the influence magnitude with squared logarithms and separately stress its covariance by a factor four. This is a planning scenario, not a confidence bound or a prediction of the unknown larger-size law.

#table(columns: (1.2fr, 1fr, 1.35fr, 1fr), inset: 4pt, stroke: 0.4pt,
 [Size window], [Pairs/size/domain], [Approx. 95% half-width for $Delta c$], [Square Gamma CPU-hours],
 [16–128, 7 sizes], [4,096], [0.319–0.526], [19.2],
 [16–128, 7 sizes], [8,192], [0.225–0.372], [38.4],
 [16–128, 7 sizes], [16,384], [0.159–0.263], [76.7],
 [32–128, 5 sizes], [8,192], [0.586–0.946], [38.1],
 [32–256, 7 sizes], [8,192], [0.320–0.604], [306.8],
)

CPU-hours extrapolate the old side-64 runtime cubically and exclude controls, validation, raw verification and queue waits. They are aggregate computation, not promised elapsed times. Raising the lower cutoff can sharply reduce precision by narrowing the log-distance range. Pilot timings must inform the final resource budget.

The next chapter records the frozen 258,048-pair design, passed validation and submitted jobs. Report every outcome; do not extend sampling merely to obtain a preferred sign or p-value.
#source("docs/history/LARGER_CAMPAIGN_20261007.md; docs/validation/larger_feasibility_20261007/precision_planning.json; plan.json; hamilton_status.json")


#pagebreak()
= Frozen larger campaign: 8 October 2026
#note[*Validated on Hamilton and submitted.* The new cohort contains 258,048 paired environments. Validation job 18965117 passed in 9 minutes 10 seconds. Production arrays 18965699 and 18965700 are submitted, followed by analysis 18965701. Results remain pending.]

#table(columns: (1fr, 3fr), inset: 4pt, stroke: 0.4pt,
 [Domains and sizes], [Aztec and ordinary square: 16, 24, 32, 48, 64, 96, 128.],
 [Gamma law], [Original independent shape-0.2 and shape-0.25 families, scale one; other weights one.],
 [Pairs per size/domain], [16,384 Gamma; 1,024 family-mean controls; 1,024 all-one controls.],
 [Observable], [Two symmetric physical diagonal increments, fractions 1/8 and 1/4; common height scale one.],
 [Primary comparison], [$Delta c=c_(Aztec)-c_(square)$ over 16–128; fixed sensitivity windows 32–128 and 48–128.],
 [Uncertainty], [10,000 whole-environment bootstrap resamples; report every window.],
 [Control diagnostic], [9,999 whole-vector permutations; maximum standardized absolute curvature contrast across two controls and three windows.],
 [Batching], [504 independent batches of 512 pairs; at most two concurrent sampling tasks.],
)

The new source archive contains 164 files. All 155 files from the previous verified sampling snapshot are unchanged. The additions implement batching, explicit retries, receipt checks and the larger analysis; the Julia wrapper extends allowed sizes and counts without changing the sampling kernel. This remains shuffling with the documented Float64 qualifications, not a CFTP production run.

Each batch has its own fixed master seed. All 753,648 checked role seeds across old, pilot, confirmation and new scientific/validation cohorts are distinct. Interrupted attempts are preserved. An explicit retry regenerates the same fixed streams and checks any completed records for agreement; it never replaces a difficult environment. Completed batches are reused after checking their receipts. Workers reload raw matchings and recompute observables before certification. Final analysis requires every batch and rechecks its certificates and compact data; missing batches block fitting.

Local verification passed all 152 Python tests. A small integration run verified eight pairs through both language implementations, and the new control diagnostic reproduced the old six contrasts and permutation p-value 0.0093. Hamilton validation passed all 42 domain/model/size cells with four pairs each, all 152 Python tests and 21,970 Julia assertions, completed-batch reuse and the validation-only analysis. Exit status was zero and stderr empty. All 336 matchings and 672 observable values passed raw verification. Validation observations are excluded from production.

Sampling requests one CPU and 4 GiB per task. The 360 smaller batches have 30-minute limits; the 144 larger batches have four-hour limits and start only after the smaller array succeeds. The final analysis has one CPU, 8 GiB and eight hours. Pilot-based estimates suggest roughly 105 CPU-hours for square Gamma sampling and 12 for square controls, plus Aztec sampling, startup and verification. With two workers, allow several days plus queue time; these are estimates, not completion guarantees.

The fixed cohort is independent of earlier experiments. It tests a wider finite window with greater precision, not the historical 128–1300 regime or an asymptotic theorem. A control discrepancy will be reported alongside the Gamma result, without selective deletion or withholding of the scientific outcome. Source hashes, configuration hashes, checks and AI assistance are recorded in the repository.
#source("docs/history/LARGER_CAMPAIGN_20261007.md; docs/validation/larger_matched_20261007/; square_glauber_python/configs/larger_matched_20261007.json; docs/AI_CONTRIBUTIONS.md")


#pagebreak()
= Larger results and public data: 10 October 2026
#note[*Both geometries show a positive fitted squared-log contribution over sizes 16–128.* The main scientific results are the coefficient for each model. The larger experiment contains 258,048 paired environments, including controls, and all 504 batches passed verification.]

Use the original Gamma law: independent families of shape 0.2 and 0.25, scale one, remaining weights one. Each fixed environment has two conditionally independent tilings. Spatial increments use two matched diagonal probes, with fractions 1/8 and 1/4, common physical distances and height units. The disorder contribution is $D=Cov(X_1,X_2)$. Fit $D_p(d)=a_p+b_p log d+c(log d)^2$ separately in each geometry, allowing probe-specific intercepts and ordinary-log slopes.

#table(columns: (1fr, 1.5fr, 1.5fr), inset: 5pt, stroke: 0.4pt,
 [Size window], [Aztec $c$ [95% CI]], [Ordinary square $c$ [95% CI]],
 [16–128, primary], [*0.17994* [0.07819, 0.27740]], [*0.15882* [0.06446, 0.25529]],
 [32–128], [0.24780 [−0.01861, 0.51194]], [0.07641 [−0.17496, 0.32768]],
 [48–128], [0.31780 [−0.32242, 0.94459]], [0.29556 [−0.30549, 0.90156]],
)

Both primary intervals exclude zero. This supports positive finite-range squared-log curvature on the ordinary square grid using the Aztec parameters, as well as on the Aztec diamond. Higher-cutoff intervals include zero, so persistence at large scales and the asymptotic physics conjecture remain unresolved. The original domain contrast is retained as an additional diagnostic; it is not the criterion for a positive coefficient in the square model. Its primary estimate is 0.02112 with interval [−0.11813, 0.15902].

The predeclared joint six-contrast control permutation test gives $p=0.2039$. The all-one primary contrast has a pointwise interval narrowly excluding zero: −0.39645 [−0.77829, −0.01349]. Both findings are retained; the joint non-rejection does not prove correctness. No data were selectively excluded or added after inspecting the outcomes.

Sampling arrays 18965699 and 18965700 and analysis 18965701 completed with exit zero and empty error logs. Analysis finished at 16:56:17 on 10 October. The checksummed compact archive was retrieved. Local verification confirmed all 504 batch records and reproduced 252 moments, 18 coefficients and nine contrasts, with maximum discrepancies $1.78 times 10^(-15)$ and $2.44 times 10^(-14)$ for moments and coefficients. Original uncertainty outputs were retained during the first verification. The public-data workflow subsequently reran all 10,000 bootstrap draws and 9,999 permutations: coefficient and interval discrepancies were below $3.83 times 10^(-14)$, moment intervals agreed exactly, and the joint p-value was again 0.2039.

== Simple reproduction and repository organisation
The public data bundle contains all 516,096 measurement rows, full pair identifiers and seeds, frozen plan, recorded intervals, provenance and exact Hamilton source snapshot. Full weights and tiling arrays remain on Hamilton, where workers verified them before issuing certificates. Sampling used corrected shuffling with the recorded numerical qualifications; this is not a CFTP production run.

Run #raw("make analysis") to verify the public checksums and recompute every point estimate and the figure using recorded uncertainty. Run #raw("make bootstrap") to rerun the fixed 10,000 whole-environment resamples and 9,999 control permutations. The README leads with both coefficients, and the analysis guide explains the equations and code. Historical files are preserved locally and in Git history. Alberto has emailed Sunil; human scientific review remains the next step.
#source("docs/RESULTS.md; analysis/README.md; data/README.md; data/verification.json; docs/AI_CONTRIBUTIONS.md")
