"""
Reference Glauber dynamics for weighted domino tilings of an even square grid.

The state is a height matrix on the faces of a `2L × 2L` dimer region, with
the tileable boundary convention preserved in
`aztec/reference/glauber_reference.jl`.
Across a dimer edge the height difference has absolute value three; across an
unoccupied edge it has absolute value one.  A local height update is therefore
exactly a flip of two opposite dimers around one square face.

This module deliberately implements *literal random-face* heat-bath updates.
It does not use the tempting "choose a flippable site" jump chain: that chain
has state-dependent holding rates and must not be used for unweighted
stationary samples without a separate rejection-free proof.
"""
module GlauberSquareGrid

using LinearAlgebra
using Random
using Statistics

export EdgeWeights,
       max_height_configuration,
       min_height_configuration,
       constant_edge_weights,
       random_edge_weights,
       validate_height_configuration,
       local_configurations,
       local_edge_weights,
       heatbath_update!,
       run_updates!,
       AcceleratedChain,
       accelerated_chain,
       run_accelerated_updates!,
       center_height,
       enumerate_height_configurations,
       matching_weight,
       matching_log_weight,
       kasteleyn_matrix,
       height_difference_moments_kasteleyn,
       center_height_moments_kasteleyn,
       tempered_edge_weights,
       replica_exchange_log_ratio,
       parallel_tempering_swap!,
       run_parallel_tempering_updates!,
       sample_center_height_chain_parallel_tempering,
       exact_center_distribution,
       detailed_balance_residual,
       integrated_autocorrelation_time,
       effective_sample_size,
       sample_center_height_chain,
       sample_center_height_chain_accelerated,
       compare_extremal_starts,
       sample_center_height_pair

"Positive weights on all height-gradient edges of a `2L × 2L` dimer region."
struct EdgeWeights
    # `vertical[i,j]` lies between height entries `(i,j)` and `(i+1,j)`.
    vertical::Matrix{Float64}
    # `horizontal[i,j]` lies between height entries `(i,j)` and `(i,j+1)`.
    horizontal::Matrix{Float64}

    function EdgeWeights(vertical::AbstractMatrix{<:Real}, horizontal::AbstractMatrix{<:Real})
        size(vertical, 2) == size(horizontal, 1) || throw(ArgumentError(
            "vertical/horizontal height-edge dimensions are incompatible"))
        size(vertical, 1) == size(horizontal, 2) || throw(ArgumentError(
            "vertical/horizontal height-edge dimensions are incompatible"))
        converted_vertical = Matrix{Float64}(vertical)
        converted_horizontal = Matrix{Float64}(horizontal)
        all(isfinite, converted_vertical) && all(>(0), converted_vertical) || throw(ArgumentError(
            "vertical edge weights must be finite and strictly positive"))
        all(isfinite, converted_horizontal) && all(>(0), converted_horizontal) || throw(ArgumentError(
            "horizontal edge weights must be finite and strictly positive"))
        new(converted_vertical, converted_horizontal)
    end
end

"Boundary helper from the supplied historical code, without distributed globals."
function boundary_corner_value(i::Int, j::Int, n::Int)
    if i == 0 || i == n
        return mod(j, 2)
    elseif j == 0 || j == n
        return -mod(i, 2)
    end
    return 0
end

function extremal_entry(i::Int, j::Int, k::Int, L::Int, sign::Int)
    side = 2L
    if i == k || i == side - k || j == k || j == side - k
        return sign * 2k + boundary_corner_value(i - k, j - k, side - 2k)
    end
    return 0
end

"Maximum tileable height configuration for a `2L × 2L` dimer square."
function max_height_configuration(L::Integer)
    L > 0 || throw(ArgumentError("L must be positive"))
    order = Int(L)
    side = 2order
    return reshape([
        extremal_entry(i, j, min(i, j, side - i, side - j), order, 1)
        for i in 0:side for j in 0:side
    ], side + 1, side + 1)
end

"Minimum tileable height configuration for a `2L × 2L` dimer square."
function min_height_configuration(L::Integer)
    L > 0 || throw(ArgumentError("L must be positive"))
    order = Int(L)
    side = 2order
    return reshape([
        extremal_entry(i, j, min(i, j, side - i, side - j), order, -1)
        for i in 0:side for j in 0:side
    ], side + 1, side + 1)
end

function check_weights_for_L(weights::EdgeWeights, L::Integer)
    side = 2Int(L)
    size(weights.vertical) == (side, side + 1) || throw(ArgumentError(
        "vertical weights do not match L=$L"))
    size(weights.horizontal) == (side + 1, side) || throw(ArgumentError(
        "horizontal weights do not match L=$L"))
    return nothing
end

"All-one edge weights: the no-disorder dimer control."
function constant_edge_weights(L::Integer)
    L > 0 || throw(ArgumentError("L must be positive"))
    side = 2Int(L)
    return EdgeWeights(ones(side, side + 1), ones(side + 1, side))
end

# Marsaglia--Tsang Gamma sampler, kept local so this reference module has no
# dependency on a particular random-weight implementation elsewhere in the repo.
function rand_gamma(rng::AbstractRNG, shape::Float64, scale::Float64)
    if shape < 1
        uniform_draw = rand(rng)
        while uniform_draw == 0
            uniform_draw = rand(rng)
        end
        return rand_gamma(rng, shape + 1, scale) * uniform_draw^(1 / shape)
    end
    d = shape - 1 / 3
    c = inv(sqrt(9d))
    while true
        normal_draw = randn(rng)
        base = 1 + c * normal_draw
        base > 0 || continue
        candidate = base^3
        uniform_draw = rand(rng)
        if uniform_draw < 1 - 0.0331 * normal_draw^4 ||
           log(uniform_draw) < normal_draw^2 / 2 + d * (1 - candidate + log(candidate))
            return scale * d * candidate
        end
    end
end

"""
    random_edge_weights(rng, L; distribution=:gamma, parameter=0.5)

Draw one frozen i.i.d. edge environment.  Every supported law has mean one:
Gamma uses `(shape=parameter, scale=1/parameter)`, lognormal uses log-scale
standard deviation `parameter`, and uniform uses `Uniform(0, 2)`.
"""
function random_edge_weights(
    rng::AbstractRNG,
    L::Integer;
    distribution::Symbol=:gamma,
    parameter::Real=0.5,
)
    L > 0 || throw(ArgumentError("L must be positive"))
    parameter > 0 || throw(ArgumentError("distribution parameter must be positive"))
    distribution in (:gamma, :lognormal, :uniform) || throw(ArgumentError(
        "distribution must be :gamma, :lognormal, or :uniform"))
    draw = if distribution === :gamma
        () -> rand_gamma(rng, Float64(parameter), inv(Float64(parameter)))
    elseif distribution === :lognormal
        () -> exp(Float64(parameter) * randn(rng) - Float64(parameter)^2 / 2)
    else
        () -> 2 * rand(rng)
    end
    side = 2Int(L)
    return EdgeWeights(
        [draw() for _ in 1:side, _ in 1:(side + 1)],
        [draw() for _ in 1:(side + 1), _ in 1:side],
    )
end

@inline is_dimer_difference(value::Integer) = abs(value) == 3
@inline is_legal_difference(value::Integer) = abs(value) == 1 || abs(value) == 3

"Return the four height differences around one primal vertex/height plaquette."
function plaquette_differences(height::AbstractMatrix{<:Integer}, i::Int, j::Int)
    return (
        height[i, j + 1] - height[i, j],
        height[i + 1, j + 1] - height[i + 1, j],
        height[i + 1, j] - height[i, j],
        height[i + 1, j + 1] - height[i, j + 1],
    )
end

"""
    validate_height_configuration(height, L)

Check the fixed boundary and the dimer constraint: every primal vertex has
exactly one incident `|Δh| = 3` edge.  This is an independent matching-level
invariant of the height encoding.
"""
function validate_height_configuration(height::AbstractMatrix{<:Integer}, L::Integer)
    L > 0 || return (valid=false, reason="L must be positive", dimer_edges=0)
    side = 2Int(L)
    size(height) == (side + 1, side + 1) || return (
        valid=false, reason="wrong height matrix dimensions", dimer_edges=0)
    boundary = max_height_configuration(L)
    boundary_ok = height[1, :] == boundary[1, :] &&
                  height[end, :] == boundary[end, :] &&
                  height[:, 1] == boundary[:, 1] &&
                  height[:, end] == boundary[:, end]
    boundary_ok || return (valid=false, reason="boundary condition changed", dimer_edges=0)

    dimer_edge_incidence = 0
    for i in 1:side, j in 1:side
        differences = plaquette_differences(height, i, j)
        all(is_legal_difference, differences) || return (
            valid=false, reason="non-dimer height difference", dimer_edges=0)
        local_dimers = count(is_dimer_difference, differences)
        local_dimers == 1 || return (
            valid=false, reason="a primal vertex does not have one dimer", dimer_edges=0)
        dimer_edge_incidence += local_dimers
    end
    # Every dimer is incident to two primal vertices.
    return (valid=true, reason="ok", dimer_edges=dimer_edge_incidence ÷ 2)
end

"Return local edge weights `(a,b,c,d)` in clockwise top-first order."
function local_edge_weights(weights::EdgeWeights, i::Int, j::Int)
    return (
        weights.vertical[i - 1, j],     # top
        weights.horizontal[i, j],       # right
        weights.vertical[i, j],         # bottom
        weights.horizontal[i, j - 1],   # left
    )
end

"""
    local_configurations(height, i, j)

For one interior face-height entry, return the two legal local alternatives as
`(ac=..., bd=...)`, where `ac` occupies the top/bottom edges and `bd` the
right/left edges.  Return `nothing` when the face is not flippable.
"""
function local_configurations(height::AbstractMatrix{<:Integer}, i::Int, j::Int)
    2 <= i < size(height, 1) && 2 <= j < size(height, 2) || throw(BoundsError(height, (i, j)))
    north, east, south, west = height[i - 1, j], height[i, j + 1],
                               height[i + 1, j], height[i, j - 1]
    ac = nothing
    bd = nothing
    for candidate in unique((north - 3, north - 1, north + 1, north + 3))
        differences = (candidate - north, candidate - east,
                       candidate - south, candidate - west)
        all(is_legal_difference, differences) || continue
        occupied = is_dimer_difference.(differences)
        if occupied == (true, false, true, false)
            ac = candidate
        elseif occupied == (false, true, false, true)
            bd = candidate
        end
    end
    if isnothing(ac) || isnothing(bd) || height[i, j] != ac && height[i, j] != bd
        return nothing
    end
    return (ac=ac::Int, bd=bd::Int)
end

"""
    heatbath_update!(rng, height, weights)

Choose one interior face uniformly.  If it is flippable, resample its two
local dimer pairings with probabilities `ac/(ac+bd)` and `bd/(ac+bd)`.
The returned named tuple records the attempted face and whether it was
flippable/changed; an update that redraws the current local state is valid.
"""
function heatbath_update!(rng::AbstractRNG, height::Matrix{Int}, weights::EdgeWeights)
    side = size(height, 1) - 1
    iseven(side) && side >= 2 || throw(ArgumentError("height side must be positive and even"))
    L = side ÷ 2
    check_weights_for_L(weights, L)
    i = rand(rng, 2:side)
    j = rand(rng, 2:side)
    configurations = local_configurations(height, i, j)
    isnothing(configurations) && return (i=i, j=j, flippable=false, changed=false)
    a, b, c, d = local_edge_weights(weights, i, j)
    probability_ac = (a * c) / (a * c + b * d)
    old_height = height[i, j]
    height[i, j] = rand(rng) < probability_ac ? configurations.ac : configurations.bd
    return (i=i, j=j, flippable=true, changed=height[i, j] != old_height)
end

"Run a prescribed number of literal random-face heat-bath attempts in place."
function run_updates!(rng::AbstractRNG, height::Matrix{Int}, weights::EdgeWeights, attempts::Integer)
    attempts >= 0 || throw(ArgumentError("attempts must be nonnegative"))
    flippable = 0
    changed = 0
    for _ in 1:Int(attempts)
        result = heatbath_update!(rng, height, weights)
        flippable += result.flippable
        changed += result.changed
    end
    return (attempts=Int(attempts), flippable=flippable, changed=changed)
end

"Fenwick-tree update for exact rejection-free Glauber event selection."
function fenwick_add!(tree::Vector{Float64}, index::Int, delta::Float64)
    while index < length(tree)
        tree[index] += delta
        index += index & -index
    end
    return tree
end

"Find the first one-based rate index whose Fenwick prefix exceeds `target`."
function fenwick_find(tree::Vector{Float64}, target::Float64)
    target >= 0 || throw(ArgumentError("Fenwick target must be nonnegative"))
    index = 0
    bit = 1
    while bit << 1 < length(tree)
        bit <<= 1
    end
    while bit != 0
        candidate = index + bit
        if candidate < length(tree) && tree[candidate] <= target
            index = candidate
            target -= tree[candidate]
        end
        bit >>= 1
    end
    return min(index + 1, length(tree) - 1)
end

"""
Mutable cache for an exact rejection-free realization of random-face Glauber
dynamics.  `rates[f]` is the probability that a literal heat-bath attempt at
face `f` changes the state, conditional on choosing that face.
"""
mutable struct AcceleratedChain
    height::Matrix{Int}
    weights::EdgeWeights
    side::Int
    rates::Vector{Float64}
    fenwick::Vector{Float64}
    total_rate::Float64
end

@inline face_count(chain::AcceleratedChain) = (chain.side - 1)^2
@inline face_index(side::Int, i::Int, j::Int) = (i - 2) * (side - 1) + j - 1
@inline function face_coordinates(side::Int, index::Int)
    quotient, remainder = divrem(index - 1, side - 1)
    return quotient + 2, remainder + 2
end

function local_change_rate(height::Matrix{Int}, weights::EdgeWeights, i::Int, j::Int)
    configurations = local_configurations(height, i, j)
    isnothing(configurations) && return 0.0
    a, b, c, d = local_edge_weights(weights, i, j)
    probability_ac = (a * c) / (a * c + b * d)
    return height[i, j] == configurations.ac ? 1 - probability_ac : probability_ac
end

function refresh_face_rate!(chain::AcceleratedChain, i::Int, j::Int)
    2 <= i <= chain.side && 2 <= j <= chain.side || return chain
    index = face_index(chain.side, i, j)
    old_rate = chain.rates[index]
    new_rate = local_change_rate(chain.height, chain.weights, i, j)
    delta = new_rate - old_rate
    delta == 0 && return chain
    chain.rates[index] = new_rate
    fenwick_add!(chain.fenwick, index, delta)
    chain.total_rate += delta
    return chain
end

"Create an accelerated cache from one valid initial height configuration."
function accelerated_chain(height::Matrix{Int}, weights::EdgeWeights)
    side = size(height, 1) - 1
    iseven(side) && side >= 2 || throw(ArgumentError("height side must be positive and even"))
    L = side ÷ 2
    check_weights_for_L(weights, L)
    validate_height_configuration(height, L).valid || throw(ArgumentError(
        "accelerated chain requires a valid height configuration"))
    rates = zeros(Float64, (side - 1)^2)
    fenwick = zeros(Float64, length(rates) + 1)
    chain = AcceleratedChain(height, weights, side, rates, fenwick, 0.0)
    for i in 2:side, j in 2:side
        refresh_face_rate!(chain, i, j)
    end
    return chain
end

"Draw the number of literal random-face attempts through the next state change."
function geometric_wait(rng::AbstractRNG, probability::Float64)
    0 < probability <= 1 || throw(ArgumentError("geometric probability must lie in (0,1]"))
    probability == 1 && return 1
    uniform_draw = rand(rng)
    while uniform_draw == 0
        uniform_draw = rand(rng)
    end
    return floor(Int, log(uniform_draw) / log1p(-probability)) + 1
end

"Perform one guaranteed state-changing event selected at its exact conditional rate."
function accelerated_event!(rng::AbstractRNG, chain::AcceleratedChain)
    chain.total_rate > 0 || return false
    # `rand` is strictly below one, so the final positive-rate face remains selectable.
    index = fenwick_find(chain.fenwick, rand(rng) * chain.total_rate)
    i, j = face_coordinates(chain.side, index)
    configurations = local_configurations(chain.height, i, j)
    isnothing(configurations) && error("stale accelerated face cache at ($i,$j)")
    previous = chain.height[i, j]
    chain.height[i, j] = previous == configurations.ac ? configurations.bd : configurations.ac
    # Only this face and its four neighbours depend on the changed height.
    refresh_face_rate!(chain, i, j)
    refresh_face_rate!(chain, i - 1, j)
    refresh_face_rate!(chain, i + 1, j)
    refresh_face_rate!(chain, i, j - 1)
    refresh_face_rate!(chain, i, j + 1)
    return true
end

"""
    run_accelerated_updates!(rng, chain, attempts)

Advance exactly `attempts` *literal random-face update times* while skipping
only self-loops.  Thus this has the same transition law at fixed attempted
times as `run_updates!`, unlike an unweighted active-site jump chain.
"""
function run_accelerated_updates!(rng::AbstractRNG, chain::AcceleratedChain, attempts::Integer)
    attempts >= 0 || throw(ArgumentError("attempts must be nonnegative"))
    remaining = Int(attempts)
    changed = 0
    while remaining > 0 && chain.total_rate > 0
        probability = chain.total_rate / face_count(chain)
        wait = geometric_wait(rng, probability)
        if wait > remaining
            remaining = 0
            break
        end
        remaining -= wait
        accelerated_event!(rng, chain) || error("nonpositive accelerated event rate")
        changed += 1
    end
    return (attempts=Int(attempts), changed=changed,
            changed_rate=changed / max(Int(attempts), 1),
            total_rate=chain.total_rate)
end

"Height of the central face in the supplied boundary convention."
function center_height(height::AbstractMatrix{<:Integer})
    side = size(height, 1) - 1
    size(height, 2) == side + 1 && iseven(side) || throw(ArgumentError(
        "height must be a square matrix with an even face side"))
    return height[side ÷ 2 + 1, side ÷ 2 + 1]
end

@inline height_key(height::AbstractMatrix{<:Integer}) = Tuple(vec(height))

"""
    enumerate_height_configurations(L)

Enumerate every domino tiling state by breadth-first local flips from the
maximum height.  This is intentionally restricted to tiny boxes and provides
the independent reference distribution used by tests.
"""
function enumerate_height_configurations(L::Integer)
    L > 0 || throw(ArgumentError("L must be positive"))
    initial = max_height_configuration(L)
    states = Dict{Tuple{Vararg{Int}},Matrix{Int}}(height_key(initial) => initial)
    queue = Matrix{Int}[initial]
    side = 2Int(L)
    cursor = 1
    while cursor <= length(queue)
        state = queue[cursor]
        cursor += 1
        for i in 2:side, j in 2:side
            configurations = local_configurations(state, i, j)
            isnothing(configurations) && continue
            for value in (configurations.ac, configurations.bd)
                value == state[i, j] && continue
                neighbour = copy(state)
                neighbour[i, j] = value
                validation = validate_height_configuration(neighbour, L)
                validation.valid || throw(ErrorException(
                    "local flip broke height/dimer invariant: $(validation.reason)"))
                key = height_key(neighbour)
                if !haskey(states, key)
                    states[key] = neighbour
                    push!(queue, neighbour)
                end
            end
        end
    end
    return collect(values(states))
end

"Product of weights of the dimers encoded by a valid height configuration."
function matching_weight(height::AbstractMatrix{<:Integer}, weights::EdgeWeights)
    side = size(height, 1) - 1
    iseven(side) || throw(ArgumentError("height side must be even"))
    L = side ÷ 2
    check_weights_for_L(weights, L)
    validate_height_configuration(height, L).valid || throw(ArgumentError(
        "height is not a valid dimer configuration"))
    product_weight = 1.0
    # A gradient edge is part of the dimer graph iff it touches at least one
    # interior face-height entry.  Boundary-to-boundary matrix segments only
    # encode the prescribed exterior boundary and are not dimer edges.
    for i in 1:side, j in 1:(side + 1)
        (2 <= i || 2 <= i + 1) && 2 <= j <= side || continue
        is_dimer_difference(height[i + 1, j] - height[i, j]) || continue
        product_weight *= weights.vertical[i, j]
    end
    for i in 1:(side + 1), j in 1:side
        2 <= i <= side && (2 <= j || 2 <= j + 1) || continue
        is_dimer_difference(height[i, j + 1] - height[i, j]) || continue
        product_weight *= weights.horizontal[i, j]
    end
    return product_weight
end

"Logarithm of `matching_weight`, avoiding underflow for strong disorder."
function matching_log_weight(height::AbstractMatrix{<:Integer}, weights::EdgeWeights)
    side = size(height, 1) - 1
    iseven(side) || throw(ArgumentError("height side must be even"))
    L = side ÷ 2
    check_weights_for_L(weights, L)
    validate_height_configuration(height, L).valid || throw(ArgumentError(
        "height is not a valid dimer configuration"))
    log_weight = 0.0
    for i in 1:side, j in 1:(side + 1)
        (2 <= i || 2 <= i + 1) && 2 <= j <= side || continue
        is_dimer_difference(height[i + 1, j] - height[i, j]) || continue
        log_weight += log(weights.vertical[i, j])
    end
    for i in 1:(side + 1), j in 1:side
        2 <= i <= side && (2 <= j || 2 <= j + 1) || continue
        is_dimer_difference(height[i, j + 1] - height[i, j]) || continue
        log_weight += log(weights.horizontal[i, j])
    end
    return log_weight
end

"Infer `L` from a complete square-grid edge environment and validate its shape."
function edge_weight_order(weights::EdgeWeights)
    side = size(weights.vertical, 1)
    iseven(side) && side >= 2 || throw(ArgumentError(
        "vertical weights must describe a positive even-sided square"))
    L = side ÷ 2
    check_weights_for_L(weights, L)
    return L
end

"""
    kasteleyn_matrix(weights)

Construct the finite bipartite Kasteleyn matrix for the same `2L × 2L`
weighted square-grid dimer model used by the Glauber sampler. Rows index
checkerboard-even primal vertices and columns index checkerboard-odd vertices.
Horizontal primal edges are real and vertical primal edges carry phase `im`.

The absolute determinant is the weighted dimer partition function. This dense
reference construction is intended for validation and moderate sizes; a sparse
nested-dissection implementation is required before much larger grids.
"""
function kasteleyn_matrix(
    weights::EdgeWeights;
    number_type::Type{T}=Float64,
) where {T<:AbstractFloat}
    L = edge_weight_order(weights)
    side = 2L
    black_vertices = [(row, column) for row in 1:side for column in 1:side
                      if iseven(row + column)]
    white_vertices = [(row, column) for row in 1:side for column in 1:side
                      if isodd(row + column)]
    black_index = Dict(vertex => index for (index, vertex) in enumerate(black_vertices))
    white_index = Dict(vertex => index for (index, vertex) in enumerate(white_vertices))
    matrix = zeros(Complex{T}, length(black_vertices), length(white_vertices))

    for ((row, column), row_index) in black_index
        if column > 1
            matrix[row_index, white_index[(row, column - 1)]] =
                T(weights.vertical[row, column])
        end
        if column < side
            matrix[row_index, white_index[(row, column + 1)]] =
                T(weights.vertical[row, column + 1])
        end
        if row > 1
            matrix[row_index, white_index[(row - 1, column)]] =
                complex(zero(T), T(weights.horizontal[row, column]))
        end
        if row < side
            matrix[row_index, white_index[(row + 1, column)]] =
                complex(zero(T), T(weights.horizontal[row + 1, column]))
        end
    end
    return (
        matrix=matrix,
        black_vertices=black_vertices,
        white_vertices=white_vertices,
        black_index=black_index,
        white_index=white_index,
    )
end

height_point(point::CartesianIndex{2}) = Tuple(point)
function height_point(point::Tuple{<:Integer,<:Integer})
    return (Int(point[1]), Int(point[2]))
end
height_point(point) = throw(ArgumentError(
    "height paths must contain CartesianIndex{2} values or integer pairs"))

"Map one directed dual-height step to its crossed undirected primal edge."
function crossed_primal_edge(first_point, second_point, side::Int)
    first_row, first_column = first_point
    second_row, second_column = second_point
    row_step = second_row - first_row
    column_step = second_column - first_column
    abs(row_step) + abs(column_step) == 1 || throw(ArgumentError(
        "successive height-path points must be nearest neighbours"))

    if row_step != 0
        row = min(first_row, second_row)
        column = first_column
        1 <= row <= side && 2 <= column <= side || throw(ArgumentError(
            "height path crosses an exterior edge outside the dimer graph"))
        return ((row, column - 1), (row, column))
    end

    row = first_row
    column = min(first_column, second_column)
    2 <= row <= side && 1 <= column <= side || throw(ArgumentError(
        "height path crosses an exterior edge outside the dimer graph"))
    return ((row - 1, column), (row, column))
end

"""
    height_difference_moments_kasteleyn(weights, path)

Compute the finite-volume conditional mean and variance of
`height[path[end]] - height[path[1]]` from Kasteleyn edge correlations. Path
coordinates are one-based entries of the `(2L+1) × (2L+1)` height matrix and
successive entries must be nearest neighbours.

The result is deterministic for a frozen environment. Its numerical
diagnostics expose the selected linear-solve residual and imaginary residuals
that should vanish in the real occupation probabilities and covariances.
"""
function height_difference_moments_kasteleyn(
    weights::EdgeWeights,
    path;
    number_type::Type{T}=Float64,
) where {T<:AbstractFloat}
    L = edge_weight_order(weights)
    side = 2L
    points = height_point.(collect(path))
    length(points) >= 2 || throw(ArgumentError("height path needs at least two points"))
    all(point -> 1 <= point[1] <= side + 1 && 1 <= point[2] <= side + 1,
        points) || throw(BoundsError((side + 1, side + 1), points))

    reference = max_height_configuration(L)
    system = kasteleyn_matrix(weights; number_type=number_type)
    coefficients = Dict{Tuple{Int,Int},Int}()
    reference_increment = 0

    for (first_point, second_point) in zip(points[1:(end - 1)], points[2:end])
        first_row, first_column = first_point
        second_row, second_column = second_point
        primal_first, primal_second = crossed_primal_edge(first_point, second_point, side)
        reference_difference = reference[second_row, second_column] -
                               reference[first_row, first_column]
        abs(reference_difference) in (1, 3) || error(
            "reference height has an illegal path increment")
        unoccupied_sign = abs(reference_difference) == 1 ? reference_difference :
                          -reference_difference ÷ 3
        reference_increment += unoccupied_sign

        black = iseven(sum(primal_first)) ? primal_first : primal_second
        white = black == primal_first ? primal_second : primal_first
        key = (system.black_index[black], system.white_index[white])
        coefficients[key] = get(coefficients, key, 0) - 4unoccupied_sign
    end

    filter!(pair -> !iszero(last(pair)), coefficients)
    edge_keys = sort!(collect(keys(coefficients)))
    edge_coefficients = T[coefficients[key] for key in edge_keys]
    matrix_order = size(system.matrix, 1)
    factorization = lu(system.matrix)
    log_partition, partition_phase = logabsdet(factorization)

    if isempty(edge_keys)
        return (
            mean=T(reference_increment),
            variance=zero(T),
            reference_increment=reference_increment,
            crossed_edges=0,
            matrix_order=matrix_order,
            log_partition=log_partition,
            partition_phase=partition_phase,
            edge_probabilities=T[],
            relative_solve_residual=zero(T),
            maximum_probability_imaginary_residual=zero(T),
            maximum_covariance_imaginary_residual=zero(T),
        )
    end

    black_columns = sort!(unique(first(key) for key in edge_keys))
    black_column_index = Dict(index => column for (column, index) in enumerate(black_columns))
    right_hand_side = zeros(Complex{T}, matrix_order, length(black_columns))
    for (column, black) in enumerate(black_columns)
        right_hand_side[black, column] = 1
    end
    selected_inverse = factorization \ right_hand_side
    relative_solve_residual = norm(system.matrix * selected_inverse - right_hand_side) /
                              max(norm(right_hand_side), eps(T))

    edge_count = length(edge_keys)
    kernel = Matrix{Complex{T}}(undef, edge_count, edge_count)
    for first_edge in 1:edge_count, second_edge in 1:edge_count
        black, white = edge_keys[first_edge]
        second_black = first(edge_keys[second_edge])
        kernel[first_edge, second_edge] = system.matrix[black, white] *
            selected_inverse[white, black_column_index[second_black]]
    end

    complex_probabilities = diag(kernel)
    maximum_probability_imaginary_residual = maximum(abs, imag.(complex_probabilities))
    numerical_tolerance = max(T(100) * eps(T) * matrix_order, sqrt(eps(T)))
    maximum_probability_imaginary_residual <= numerical_tolerance || error(
        "Kasteleyn occupation probabilities have a material imaginary residual")
    probabilities = real.(complex_probabilities)
    all(probability -> -numerical_tolerance <= probability <= 1 + numerical_tolerance,
        probabilities) || error(
        "Kasteleyn occupation probability lies outside [0,1]")
    probabilities = clamp.(probabilities, zero(T), one(T))

    covariance = Matrix{T}(undef, edge_count, edge_count)
    maximum_covariance_imaginary_residual = zero(T)
    for first_edge in 1:edge_count, second_edge in 1:edge_count
        complex_value = if first_edge == second_edge
            complex_probabilities[first_edge] * (1 - complex_probabilities[first_edge])
        else
            -kernel[first_edge, second_edge] * kernel[second_edge, first_edge]
        end
        maximum_covariance_imaginary_residual = max(
            maximum_covariance_imaginary_residual, abs(imag(complex_value)))
        covariance[first_edge, second_edge] = real(complex_value)
    end
    maximum_covariance_imaginary_residual <= numerical_tolerance || error(
        "Kasteleyn occupation covariances have a material imaginary residual")

    mean_difference = reference_increment + dot(edge_coefficients, probabilities)
    variance = dot(edge_coefficients, covariance * edge_coefficients)
    variance_tolerance = numerical_tolerance * max(one(T), sum(abs2, edge_coefficients))
    variance >= -variance_tolerance || error(
        "Kasteleyn height variance is materially negative")

    return (
        mean=mean_difference,
        variance=max(variance, zero(T)),
        reference_increment=reference_increment,
        crossed_edges=edge_count,
        matrix_order=matrix_order,
        log_partition=log_partition,
        partition_phase=partition_phase,
        edge_probabilities=probabilities,
        relative_solve_residual=relative_solve_residual,
        maximum_probability_imaginary_residual=maximum_probability_imaginary_residual,
        maximum_covariance_imaginary_residual=maximum_covariance_imaginary_residual,
    )
end

"""
    center_height_moments_kasteleyn(weights)

Compute the conditional mean and variance of the same central-face height used
by the Glauber campaigns, without Markov-chain sampling. The calculation uses
a straight dual path from the fixed top boundary to the central face.
"""
function center_height_moments_kasteleyn(
    weights::EdgeWeights;
    number_type::Type{T}=Float64,
) where {T<:AbstractFloat}
    L = edge_weight_order(weights)
    center_column = L + 1
    path = [(row, center_column) for row in 1:(L + 1)]
    difference = height_difference_moments_kasteleyn(
        weights, path; number_type=number_type)
    boundary_height = max_height_configuration(L)[1, center_column]
    return merge(difference, (
        mean=boundary_height + difference.mean,
        boundary_height=boundary_height,
        path=path,
    ))
end

"Weights for the Gibbs law proportional to `matching_weight(height, weights)^beta`."
function tempered_edge_weights(weights::EdgeWeights, beta::Real)
    isfinite(beta) && beta >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
    exponent = Float64(beta)
    return EdgeWeights(weights.vertical .^ exponent, weights.horizontal .^ exponent)
end

"Log Metropolis ratio for swapping states between two inverse-temperature replicas."
function replica_exchange_log_ratio(
    left_height::AbstractMatrix{<:Integer},
    left_beta::Real,
    right_height::AbstractMatrix{<:Integer},
    right_beta::Real,
    base_weights::EdgeWeights,
)
    isfinite(left_beta) && isfinite(right_beta) || throw(ArgumentError("betas must be finite"))
    return (Float64(left_beta) - Float64(right_beta)) *
           (matching_log_weight(right_height, base_weights) -
            matching_log_weight(left_height, base_weights))
end

"Rebuild cached event rates after an accepted replica exchange."
function rebuild_accelerated_cache!(chain::AcceleratedChain)
    rebuilt = accelerated_chain(chain.height, chain.weights)
    chain.rates .= rebuilt.rates
    chain.fenwick .= rebuilt.fenwick
    chain.total_rate = rebuilt.total_rate
    return chain
end

"""
    parallel_tempering_swap!(rng, left, left_beta, right, right_beta, base_weights)

Attempt an adjacent replica exchange.  Each replica retains its own tempered
Glauber kernel; the Metropolis acceptance rule preserves their joint product
law, hence the `beta=1` replica remains an exact sample from the requested
weighted-dimer Gibbs distribution.
"""
function parallel_tempering_swap!(
    rng::AbstractRNG,
    left::AcceleratedChain,
    left_beta::Real,
    right::AcceleratedChain,
    right_beta::Real,
    base_weights::EdgeWeights,
)
    size(left.height) == size(right.height) || throw(ArgumentError(
        "replica heights must have the same dimensions"))
    log_ratio = replica_exchange_log_ratio(
        left.height, left_beta, right.height, right_beta, base_weights)
    accepted = log(rand(rng)) <= min(0.0, log_ratio)
    if accepted
        temporary = copy(left.height)
        left.height .= right.height
        right.height .= temporary
        rebuild_accelerated_cache!(left)
        rebuild_accelerated_cache!(right)
    end
    return (accepted=accepted, log_acceptance_ratio=log_ratio)
end

function validate_tempering_betas(betas)
    values = Float64.(collect(betas))
    length(values) >= 2 || throw(ArgumentError("parallel tempering needs at least two replicas"))
    issorted(values) || throw(ArgumentError("tempering betas must be sorted"))
    all(beta -> 0 <= beta <= 1, values) || throw(ArgumentError(
        "tempering betas must lie in [0, 1]"))
    values[end] == 1.0 || throw(ArgumentError("the target replica must have beta=1"))
    return values
end

"""
    run_parallel_tempering_updates!(rng, chains, betas, base_weights, attempts;
                                    swap_interval_attempts)

Advance every tempered replica by the same number of literal attempted-update
times, with alternating adjacent replica exchanges.  Pass the returned
`swap_round` and `attempts_since_swap` into the next call when splitting one
chain across burn-in or retained-sample intervals.  The last (`beta=1`)
replica is the only one used for the physical observable.
"""
function run_parallel_tempering_updates!(
    rng::AbstractRNG,
    chains::AbstractVector{<:AcceleratedChain},
    betas,
    base_weights::EdgeWeights,
    attempts::Integer;
    swap_interval_attempts::Integer,
    swap_round_offset::Integer=0,
    attempts_since_swap::Integer=0,
)
    attempts >= 0 || throw(ArgumentError("attempts must be nonnegative"))
    swap_interval_attempts > 0 || throw(ArgumentError(
        "swap_interval_attempts must be positive"))
    swap_round_offset >= 0 || throw(ArgumentError("swap_round_offset must be nonnegative"))
    0 <= attempts_since_swap < swap_interval_attempts || throw(ArgumentError(
        "attempts_since_swap must lie in [0, swap_interval_attempts)"))
    values = validate_tempering_betas(betas)
    length(chains) == length(values) || throw(ArgumentError(
        "one accelerated chain is required per beta"))
    remaining = Int(attempts)
    swap_round = Int(swap_round_offset)
    since_swap = Int(attempts_since_swap)
    accepted_swaps_by_pair = zeros(Int, length(chains) - 1)
    attempted_swaps_by_pair = zeros(Int, length(chains) - 1)
    target_changed = 0
    while remaining > 0
        block = min(remaining, Int(swap_interval_attempts) - since_swap)
        for (index, chain) in enumerate(chains)
            report = run_accelerated_updates!(rng, chain, block)
            index == length(chains) && (target_changed += report.changed)
        end
        since_swap += block
        remaining -= block
        if since_swap == swap_interval_attempts
            first_index = iseven(swap_round) ? 1 : 2
            for index in first_index:2:(length(chains) - 1)
                report = parallel_tempering_swap!(rng, chains[index], values[index],
                                                   chains[index + 1], values[index + 1],
                                                   base_weights)
                attempted_swaps_by_pair[index] += 1
                accepted_swaps_by_pair[index] += report.accepted
            end
            swap_round += 1
            since_swap = 0
        end
    end
    return (attempts=Int(attempts), target_changed=target_changed,
            target_changed_rate=target_changed / max(Int(attempts), 1),
            attempted_swaps=sum(attempted_swaps_by_pair),
            accepted_swaps=sum(accepted_swaps_by_pair),
            attempted_swaps_by_pair=attempted_swaps_by_pair,
            accepted_swaps_by_pair=accepted_swaps_by_pair,
            swap_round=swap_round, attempts_since_swap=since_swap)
end

"""
Sample central heights from the `beta=1` replica of an exact parallel-tempering
chain.  Auxiliary replicas alter mixing only; their tempered laws and the
exchange kernel leave the requested `beta=1` weighted-dimer law invariant.
"""
function sample_center_height_chain_parallel_tempering(
    rng::AbstractRNG,
    L::Integer,
    base_weights::EdgeWeights;
    start::Symbol=:max,
    betas=collect(0.0:0.2:1.0),
    burn_in_attempts::Integer,
    thin_attempts::Integer,
    samples::Integer,
    swap_interval_attempts::Integer=1_000,
)
    L > 0 || throw(ArgumentError("L must be positive"))
    burn_in_attempts >= 0 || throw(ArgumentError("burn_in_attempts must be nonnegative"))
    thin_attempts > 0 || throw(ArgumentError("thin_attempts must be positive"))
    samples > 0 || throw(ArgumentError("samples must be positive"))
    check_weights_for_L(base_weights, L)
    values = validate_tempering_betas(betas)
    initial = start === :max ? max_height_configuration(L) :
              start === :min ? min_height_configuration(L) :
              throw(ArgumentError("start must be :max or :min"))
    chains = [accelerated_chain(copy(initial), tempered_edge_weights(base_weights, beta))
              for beta in values]
    burn_in = run_parallel_tempering_updates!(rng, chains, values, base_weights,
                                               burn_in_attempts;
                                               swap_interval_attempts=swap_interval_attempts)
    observed = Vector{Int}(undef, samples)
    target_changed = 0
    accepted_swaps_by_pair = copy(burn_in.accepted_swaps_by_pair)
    attempted_swaps_by_pair = copy(burn_in.attempted_swaps_by_pair)
    swap_round = burn_in.swap_round
    attempts_since_swap = burn_in.attempts_since_swap
    for index in eachindex(observed)
        report = run_parallel_tempering_updates!(rng, chains, values, base_weights,
                                                  thin_attempts;
                                                  swap_interval_attempts=swap_interval_attempts,
                                                  swap_round_offset=swap_round,
                                                  attempts_since_swap=attempts_since_swap)
        observed[index] = center_height(chains[end].height)
        target_changed += report.target_changed
        accepted_swaps_by_pair .+= report.accepted_swaps_by_pair
        attempted_swaps_by_pair .+= report.attempted_swaps_by_pair
        swap_round = report.swap_round
        attempts_since_swap = report.attempts_since_swap
    end
    accepted_swaps = sum(accepted_swaps_by_pair)
    attempted_swaps = sum(attempted_swaps_by_pair)
    pair_acceptance_rates = [attempted == 0 ? NaN : accepted / attempted
                             for (accepted, attempted) in
                             zip(accepted_swaps_by_pair, attempted_swaps_by_pair)]
    minimum_pair_swap_acceptance_rate = any(isnan, pair_acceptance_rates) ? NaN :
                                        minimum(pair_acceptance_rates)
    return (
        heights=observed,
        mean=mean(observed),
        variance=length(observed) > 1 ? var(observed; corrected=true) : 0.0,
        diagnostics=(
            integrated_autocorrelation_time=integrated_autocorrelation_time(observed),
            effective_sample_size=effective_sample_size(observed),
            target_changed_rate=target_changed / (samples * thin_attempts),
            swap_acceptance_rate=accepted_swaps / max(attempted_swaps, 1),
            attempted_swaps=attempted_swaps,
            attempted_swaps_by_pair=attempted_swaps_by_pair,
            accepted_swaps_by_pair=accepted_swaps_by_pair,
            pair_swap_acceptance_rates=pair_acceptance_rates,
            minimum_pair_swap_acceptance_rate=minimum_pair_swap_acceptance_rate,
            target_exchange_attempts=attempted_swaps_by_pair[end],
            target_exchange_acceptance_rate=pair_acceptance_rates[end],
        ),
        final_height=chains[end].height,
        burn_in=burn_in,
    )
end

"Exact tiny-grid centre-height probabilities under the weighted dimer Gibbs law."
function exact_center_distribution(L::Integer, weights::EdgeWeights)
    states = enumerate_height_configurations(L)
    raw_weights = matching_weight.(states, Ref(weights))
    normalizer = sum(raw_weights)
    probabilities = raw_weights ./ normalizer
    support = sort(unique(center_height.(states)))
    masses = Dict(value => sum(probabilities[index] for index in eachindex(states)
                               if center_height(states[index]) == value) for value in support)
    return (states=states, probabilities=probabilities, masses=masses,
            mean=sum(value * probability for (value, probability) in masses),
            variance=sum((value - sum(key * p for (key, p) in masses))^2 * probability
                         for (value, probability) in masses))
end

"Maximum absolute detailed-balance residual over every adjacent pair of tiny states."
function detailed_balance_residual(L::Integer, weights::EdgeWeights)
    states = enumerate_height_configurations(L)
    lookup = Dict(height_key(state) => index for (index, state) in enumerate(states))
    face_count = (2Int(L) - 1)^2
    residual = 0.0
    for (index, state) in enumerate(states), i in 2:(2Int(L)), j in 2:(2Int(L))
        configurations = local_configurations(state, i, j)
        isnothing(configurations) && continue
        other_value = state[i, j] == configurations.ac ? configurations.bd : configurations.ac
        other = copy(state)
        other[i, j] = other_value
        other_index = lookup[height_key(other)]
        a, b, c, d = local_edge_weights(weights, i, j)
        probability_ac = (a * c) / (a * c + b * d)
        transition = (state[i, j] == configurations.ac ? 1 - probability_ac : probability_ac) /
                     face_count
        reverse_transition = (other_value == configurations.ac ? 1 - probability_ac : probability_ac) /
                             face_count
        residual = max(residual, abs(
            matching_weight(state, weights) * transition -
            matching_weight(states[other_index], weights) * reverse_transition,
        ))
    end
    return residual
end

"""
    integrated_autocorrelation_time(values)

Initial-positive-sequence estimate for a scalar trace.  It is a diagnostic,
not a convergence proof; the pilot compares chains from both extremal starts
in addition to reporting this value.
"""
function integrated_autocorrelation_time(values::AbstractVector{<:Real})
    count = length(values)
    count <= 1 && return 1.0
    centred = Float64.(values) .- mean(values)
    denominator = sum(abs2, centred)
    denominator == 0 && return 1.0
    tau = 1.0
    # The truncation prevents a noisy long-lag tail from manufacturing an ESS.
    for lag in 1:min(count - 1, floor(Int, sqrt(count)))
        correlation = sum(centred[index] * centred[index + lag]
                          for index in 1:(count - lag)) / denominator
        correlation <= 0 && break
        tau += 2 * correlation
    end
    return max(tau, 1.0)
end

"Effective sample-size diagnostic derived from `integrated_autocorrelation_time`."
function effective_sample_size(values::AbstractVector{<:Real})
    return length(values) / integrated_autocorrelation_time(values)
end

"Run one chain and retain centre-height samples after burn-in and thinning."
function sample_center_height_chain(
    rng::AbstractRNG,
    L::Integer,
    weights::EdgeWeights;
    start::Symbol=:max,
    burn_in_attempts::Integer,
    thin_attempts::Integer,
    samples::Integer,
)
    L > 0 || throw(ArgumentError("L must be positive"))
    burn_in_attempts >= 0 || throw(ArgumentError("burn_in_attempts must be nonnegative"))
    thin_attempts > 0 || throw(ArgumentError("thin_attempts must be positive"))
    samples > 0 || throw(ArgumentError("samples must be positive"))
    check_weights_for_L(weights, L)
    height = start === :max ? max_height_configuration(L) :
             start === :min ? min_height_configuration(L) :
             throw(ArgumentError("start must be :max or :min"))
    burn_in = run_updates!(rng, height, weights, burn_in_attempts)
    observed = Vector{Int}(undef, samples)
    flippable = 0
    changed = 0
    for index in eachindex(observed)
        report = run_updates!(rng, height, weights, thin_attempts)
        observed[index] = center_height(height)
        flippable += report.flippable
        changed += report.changed
    end
    return (
        heights=observed,
        mean=mean(observed),
        variance=length(observed) > 1 ? var(observed; corrected=true) : 0.0,
        diagnostics=(
            integrated_autocorrelation_time=integrated_autocorrelation_time(observed),
            effective_sample_size=effective_sample_size(observed),
            flippable_rate=flippable / (samples * thin_attempts),
            changed_rate=changed / (samples * thin_attempts),
        ),
        final_height=height,
        burn_in=burn_in,
        sampling=(attempts=Int(samples * thin_attempts), flippable=flippable, changed=changed),
    )
end

"""
As `sample_center_height_chain`, using the exact rejection-free accelerator.
The output is sampled at the same literal attempted-update schedule, so it can
be compared directly with the reference random-face chain.
"""
function sample_center_height_chain_accelerated(
    rng::AbstractRNG,
    L::Integer,
    weights::EdgeWeights;
    start::Symbol=:max,
    burn_in_attempts::Integer,
    thin_attempts::Integer,
    samples::Integer,
)
    L > 0 || throw(ArgumentError("L must be positive"))
    burn_in_attempts >= 0 || throw(ArgumentError("burn_in_attempts must be nonnegative"))
    thin_attempts > 0 || throw(ArgumentError("thin_attempts must be positive"))
    samples > 0 || throw(ArgumentError("samples must be positive"))
    height = start === :max ? max_height_configuration(L) :
             start === :min ? min_height_configuration(L) :
             throw(ArgumentError("start must be :max or :min"))
    chain = accelerated_chain(height, weights)
    burn_in = run_accelerated_updates!(rng, chain, burn_in_attempts)
    observed = Vector{Int}(undef, samples)
    changed = 0
    for index in eachindex(observed)
        report = run_accelerated_updates!(rng, chain, thin_attempts)
        observed[index] = center_height(chain.height)
        changed += report.changed
    end
    return (
        heights=observed,
        mean=mean(observed),
        variance=length(observed) > 1 ? var(observed; corrected=true) : 0.0,
        diagnostics=(
            integrated_autocorrelation_time=integrated_autocorrelation_time(observed),
            effective_sample_size=effective_sample_size(observed),
            changed_rate=changed / (samples * thin_attempts),
            final_total_change_rate=chain.total_rate / face_count(chain),
        ),
        final_height=chain.height,
        burn_in=burn_in,
        sampling=(attempts=Int(samples * thin_attempts), changed=changed),
    )
end

"""
    compare_extremal_starts(rng_max, rng_min, L, weights; ...)

Run independent chains from the maximum and minimum valid height
configurations in one frozen environment.  Agreement is a useful pilot
diagnostic, but is deliberately not labelled a proof of equilibration.
"""
function compare_extremal_starts(
    rng_max::AbstractRNG,
    rng_min::AbstractRNG,
    L::Integer,
    weights::EdgeWeights;
    burn_in_attempts::Integer,
    thin_attempts::Integer,
    samples::Integer,
)
    maximum = sample_center_height_chain(rng_max, L, weights;
        start=:max, burn_in_attempts=burn_in_attempts,
        thin_attempts=thin_attempts, samples=samples)
    minimum = sample_center_height_chain(rng_min, L, weights;
        start=:min, burn_in_attempts=burn_in_attempts,
        thin_attempts=thin_attempts, samples=samples)
    return (maximum=maximum, minimum=minimum,
            mean_gap=maximum.mean - minimum.mean,
            final_height_gap=center_height(maximum.final_height) - center_height(minimum.final_height))
end

@inline function splitmix64(value::UInt64)
    value += 0x9e3779b97f4a7c15
    value = (value ⊻ (value >> 30)) * 0xbf58476d1ce4e5b9
    value = (value ⊻ (value >> 27)) * 0x94d049bb133111eb
    return value ⊻ (value >> 31)
end

"""
Draw one shared frozen environment and sample two independently seeded chains.
The result preserves the environmental pairing needed for later variance
decomposition, while the present pilot can inspect the marginal centre height.
"""
function sample_center_height_pair(
    sample_seed::Integer,
    L::Integer;
    distribution::Symbol=:gamma,
    parameter::Real=0.5,
    burn_in_attempts::Integer,
    thin_attempts::Integer,
    samples::Integer,
    start::Symbol=:max,
)
    base = UInt64(sample_seed)
    environment_seed = splitmix64(base ⊻ 0x656e7669726f6e6d)
    replica_1_seed = splitmix64(base ⊻ 0x7265706c69636131)
    replica_2_seed = splitmix64(base ⊻ 0x7265706c69636132)
    weights = random_edge_weights(Random.Xoshiro(environment_seed), L;
                                  distribution=distribution, parameter=parameter)
    first = sample_center_height_chain(Random.Xoshiro(replica_1_seed), L, weights;
        start=start, burn_in_attempts=burn_in_attempts,
        thin_attempts=thin_attempts, samples=samples)
    second = sample_center_height_chain(Random.Xoshiro(replica_2_seed), L, weights;
        start=start, burn_in_attempts=burn_in_attempts,
        thin_attempts=thin_attempts, samples=samples)
    return (environment_seed=environment_seed, replica_1_seed=replica_1_seed,
            replica_2_seed=replica_2_seed, weights=weights, replica_1=first, replica_2=second)
end

end
