"""
Reference Glauber dynamics for weighted domino tilings of an even square grid.

The state is a height matrix on the faces of a `2L × 2L` dimer region, with
the fixed boundary convention used in Sunil's supplied `glaubertwo.jl`.
Across a dimer edge the height difference has absolute value three; across an
unoccupied edge it has absolute value one.  A local height update is therefore
exactly a flip of two opposite dimers around one square face.

This module deliberately implements *literal random-face* heat-bath updates.
It does not use the tempting "choose a flippable site" jump chain: that chain
has state-dependent holding rates and must not be used for unweighted
stationary samples without a separate rejection-free proof.
"""
module GlauberSquareGrid

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
       center_height,
       enumerate_height_configurations,
       matching_weight,
       exact_center_distribution,
       detailed_balance_residual,
       integrated_autocorrelation_time,
       effective_sample_size,
       sample_center_height_chain,
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

"Return the local edge weights `(a,b,c,d)` in Sunil's clockwise top-first order."
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
