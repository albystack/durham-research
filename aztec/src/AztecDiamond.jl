"""
Core simulation code for weighted Aztec-diamond domino tilings.

The public API deliberately separates four stages:

1. draw or supply edge weights;
2. reduce those weights and draw creation coins;
3. replay the coins through domino shuffling;
4. evaluate or validate the resulting tiling.

Keeping these stages separate makes it possible to test the mathematical
recurrence independently of the shuffling implementation.  It also makes the
shared-environment double-dimer construction explicit: reduce one environment
once, but draw two independent sets of creation coins.
"""
module AztecDiamond

using Random
using Printf

export random_uniform_weights,
       random_gamma_weights,
       gamma_disordered_weights,
       creation_probabilities,
       gamma_disordered_probabilities,
       gamma_disordered_creation_choices,
       gamma_disordered_creation_choice_pair,
       uniform_creation_choices,
       sample_tiling,
       sample_tiling_from_choices,
       sample_gamma_center_height,
       sample_gamma_center_height_pair,
       sample_gamma_spatial_increment_pair,
       sample_uniform_spatial_increment_pair,
       height_function,
       face_height,
       center_face_index,
       center_height,
       symmetric_face_columns,
       symmetric_height_increment,
       validate_tiling,
       orientation_counts,
       write_table,
       write_svg

# This file is intentionally dependency-light.  It is arranged in the same
# order as the simulation pipeline:
#
#   weights -> reduction/creation coins -> domino shuffling -> height
#
# Validation and rendering helpers are at the end of the file.

# ---------------------------------------------------------------------------
# Random weights
# ---------------------------------------------------------------------------

"""
    random_uniform_weights(rng, order)

Return the `2order × 2order` i.i.d. `Uniform(0,1)` positive weight table
used by the generic random-weight example.
"""
function random_uniform_weights(rng::AbstractRNG, order::Integer)
    order > 0 || throw(ArgumentError("order must be positive"))
    return rand(rng, Float64, 2 * order, 2 * order)
end

function rand_gamma(rng::AbstractRNG, shape::Float64, scale::Float64)
    # Marsaglia--Tsang sampler.  The standard rejection sampler requires
    # shape >= 1.  For smaller shapes, boost by one and apply the exact
    # transformation G(shape) = G(shape+1) * U^(1/shape).  This branch matters
    # here because the production shapes are 0.2 and 0.25.
    if shape < 1
        uniform_draw = rand(rng)
        while uniform_draw == 0
            uniform_draw = rand(rng)
        end
        return rand_gamma(rng, shape + 1, scale) * uniform_draw^(1 / shape)
    end

    # Marsaglia and Tsang's cubic normal proposal for shape >= 1.
    d = shape - 1 / 3
    c = inv(sqrt(9 * d))
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

# ---------------------------------------------------------------------------
# Weight reduction and creation decisions
# ---------------------------------------------------------------------------

"""
    random_gamma_weights(rng, rows, cols; shape, scale)

Return a matrix of independent Gamma variables, with the shape/scale
parameterization `mean = shape * scale`.
"""
function random_gamma_weights(
    rng::AbstractRNG,
    rows::Integer,
    cols::Integer;
    shape::Real,
    scale::Real=1.0,
)
    rows > 0 && cols > 0 || throw(ArgumentError("matrix dimensions must be positive"))
    shape > 0 || throw(ArgumentError("Gamma shape must be positive"))
    scale > 0 || throw(ArgumentError("Gamma scale must be positive"))
    result = Matrix{Float64}(undef, rows, cols)
    float_shape = Float64(shape)
    float_scale = Float64(scale)
    for index in eachindex(result)
        result[index] = rand_gamma(rng, float_shape, float_scale)
    end
    return result
end

"""
    gamma_disordered_weights(rng, order; alpha=0.2, beta=0.25)

Sample the two independent weight families in Definition 1.1 of Duits and
Van Peski, *The Gamma-disordered Aztec diamond*:
`a[i,j] ~ Gamma(alpha, 1)` and `b[i,j] ~ Gamma(beta, 1)`.
The two other edge families are gauge-fixed to weight 1.
"""
function gamma_disordered_weights(
    rng::AbstractRNG,
    order::Integer;
    alpha::Real=0.2,
    beta::Real=0.25,
)
    order > 0 || throw(ArgumentError("order must be positive"))
    a = random_gamma_weights(rng, order, order; shape=alpha, scale=1.0)
    b = random_gamma_weights(rng, order, order; shape=beta, scale=1.0)
    return (a=a, b=b)
end

function check_weights(weights::AbstractMatrix{<:Real})
    rows, cols = size(weights)
    rows == cols || throw(ArgumentError("weight table must be square"))
    iseven(rows) || throw(ArgumentError("weight table side length must be even"))
    rows > 0 || throw(ArgumentError("weight table cannot be empty"))
    all(isfinite, weights) || throw(ArgumentError("weights must be finite"))
    all(>(0), weights) || throw(ArgumentError(
        "this random-weight implementation requires strictly positive weights",
    ))
    return rows ÷ 2
end

"""
    creation_probabilities(weights)

Compute the sequence of creation-probability matrices for domino shuffling.
This is the strictly-positive-weight specialization of `d3pslim` followed by
`probsslim` in the supplied 2022 Julia code. The first matrix is `1 × 1` and
the last is `order × order`.
"""
function creation_probabilities(weights::AbstractMatrix{<:Real})
    order = check_weights(weights)
    reductions = Vector{Matrix{Float64}}(undef, order)
    reductions[1] = Float64.(weights)

    for level in 2:order
        # Each urban-renewal step removes the outermost two rows and columns.
        # The zero-based indices below mirror the parity convention in the
        # original domino-shuffling recurrence; converting them in one place
        # avoids scattering `+1` corrections through the formula.
        previous = reductions[level - 1]
        previous_side = size(previous, 1)
        current_side = previous_side - 2
        current = Matrix{Float64}(undef, current_side, current_side)

        for i0 in 0:(current_side - 1), j0 in 0:(current_side - 1)
            i1 = i0 + 2 * (i0 % 2) + 1
            i2 = i0 + 2
            j1 = j0 + 2 * (j0 % 2) + 1
            j2 = j0 + 2

            denominator =
                previous[i1, j1] * previous[i2, j2] +
                previous[i1, j2] * previous[i2, j1]
            isfinite(denominator) && denominator > 0 ||
                throw(ArgumentError("non-finite reduction at level $level"))
            current[i0 + 1, j0 + 1] = previous[i1, j1] / denominator
        end

        reductions[level] = current
    end

    probabilities = Vector{Matrix{Float64}}(undef, order)
    for k0 in 0:(order - 1)
        reduced = reductions[order - k0]
        probability = Matrix{Float64}(undef, k0 + 1, k0 + 1)
        for i0 in 0:k0, j0 in 0:k0
            diagonal =
                reduced[2 * i0 + 2, 2 * j0 + 2] *
                reduced[2 * i0 + 1, 2 * j0 + 1]
            off_diagonal =
                reduced[2 * i0 + 2, 2 * j0 + 1] *
                reduced[2 * i0 + 1, 2 * j0 + 2]
            p = diagonal / (diagonal + off_diagonal)
            isfinite(p) && 0 <= p <= 1 ||
                throw(ArgumentError("invalid creation probability at level $(k0 + 1)"))
            probability[i0 + 1, j0 + 1] = p
        end
        probabilities[k0 + 1] = probability
    end
    return probabilities
end

"""
    gamma_disordered_probabilities(a, b)

Compute the creation probabilities for the biased Gamma-disordered model
using equations (1.22)–(1.23) of Duits and Van Peski. In this code's matrix
encoding, the returned probability is for the green/yellow pair and is
therefore `b/(a+b)`.
"""
function gamma_disordered_probabilities(
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real},
)
    order = check_gamma_weights(a, b)

    # probabilities[k] belongs to the order-k diamond generated during the
    # forward shuffle.  We start from order L and reduce backwards to order 1.
    probabilities = Vector{Matrix{Float64}}(undef, order)
    current_a = Float64.(a)
    current_b = Float64.(b)
    probabilities[order] = current_b ./ (current_a + current_b)
    for level in order:-1:2
        reduced_a, reduced_b = reduce_gamma_weights(current_a, current_b)
        all(isfinite, reduced_a) && all(isfinite, reduced_b) ||
            throw(ArgumentError("non-finite Gamma reduction at level $level"))
        probability = reduced_b ./ (reduced_a + reduced_b)
        all(p -> isfinite(p) && 0 <= p <= 1, probability) ||
            throw(ArgumentError("invalid Gamma creation probability at level $(level - 1)"))
        probabilities[level - 1] = probability
        current_a = reduced_a
        current_b = reduced_b
    end
    return probabilities
end

"""Validate the two square, positive arrays that define one Gamma environment."""
function check_gamma_weights(
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real},
)
    size(a) == size(b) || throw(ArgumentError("a and b must have the same size"))
    rows, cols = size(a)
    rows == cols || throw(ArgumentError("a and b must be square"))
    rows > 0 || throw(ArgumentError("weight tables cannot be empty"))
    all(isfinite, a) && all(isfinite, b) ||
        throw(ArgumentError("Gamma weights must be finite"))
    all(>(0), a) && all(>(0), b) ||
        throw(ArgumentError("Gamma weights must be strictly positive"))
    return rows
end

"""
Apply one Duits--Van Peski Gamma reduction step.

An `L x L` environment becomes an `(L-1) x (L-1)` environment. Keeping this
step in one function prevents the probability and bit-coin samplers from
silently drifting apart.
"""
function reduce_gamma_weights(
    current_a::AbstractMatrix{<:Real},
    current_b::AbstractMatrix{<:Real},
)
    size(current_a) == size(current_b) ||
        throw(ArgumentError("a and b must have the same size"))
    level, columns = size(current_a)
    level == columns || throw(ArgumentError("a and b must be square"))
    level > 1 || throw(ArgumentError("cannot reduce an order-one environment"))
    reduced_a = Matrix{Float64}(undef, level - 1, level - 1)
    reduced_b = Matrix{Float64}(undef, level - 1, level - 1)

    reduce_gamma_weights!(
        reduced_a,
        reduced_b,
        current_a,
        current_b,
        level,
    )
    return reduced_a, reduced_b
end

"""In-place form of one Gamma reduction step, used to reuse large buffers."""
function reduce_gamma_weights!(
    reduced_a::AbstractMatrix{Float64},
    reduced_b::AbstractMatrix{Float64},
    current_a::AbstractMatrix{<:Real},
    current_b::AbstractMatrix{<:Real},
    level::Int,
)
    size(current_a, 1) >= level && size(current_a, 2) >= level ||
        throw(ArgumentError("current a buffer is too small"))
    size(current_b, 1) >= level && size(current_b, 2) >= level ||
        throw(ArgumentError("current b buffer is too small"))
    size(reduced_a, 1) >= level - 1 && size(reduced_a, 2) >= level - 1 ||
        throw(ArgumentError("reduced a buffer is too small"))
    size(reduced_b, 1) >= level - 1 && size(reduced_b, 2) >= level - 1 ||
        throw(ArgumentError("reduced b buffer is too small"))

    # This is equations (1.22) of Duits--Van Peski written without temporary
    # arrays.  Superscript [level] is represented by current_a/current_b and
    # superscript [level-1] by reduced_a/reduced_b.
    @inbounds for i in 1:(level - 1), j in 1:(level - 1)
        reduced_a[i, j] =
            current_a[i, j] / (current_a[i, j] + current_b[i, j]) *
            (current_a[i + 1, j] + current_b[i + 1, j])
        reduced_b[i, j] =
            current_b[i, j + 1] /
            (current_a[i, j + 1] + current_b[i, j + 1]) *
            (current_a[i + 1, j + 1] + current_b[i + 1, j + 1])
    end
    return nothing
end

function draw_creation_choices!(
    rng::AbstractRNG,
    choices::AbstractMatrix{Bool},
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real},
)
    # `true` selects the green/yellow diagonal in this matrix encoding.  Its
    # probability is b/(a+b), the second probability in equation (1.23).
    for index in eachindex(choices, a, b)
        choices[index] = rand(rng) < b[index] / (a[index] + b[index])
    end
    return choices
end

"""
    gamma_disordered_creation_choices(rng, a, b)

Reduce the Gamma weights while pre-drawing every independent creation coin.
The result stores one bit per potential creation rather than a `Float64`
probability. Coins at sites which are not holes during shuffling are ignored,
so this is distributionally identical to drawing a coin only when needed.
The rolling Gamma-weight state uses `O(L^2)` memory. The complete sequence
of creation decisions uses `O(L^3)` bits, rather than `O(L^3)` Float64
probability values.
"""
function gamma_disordered_creation_choices(
    rng::AbstractRNG,
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real},
)
    return only(gamma_disordered_creation_choices(rng, a, b, 1))
end

"""
    gamma_disordered_creation_choices(rng, a, b, copies)

Draw `copies` conditionally independent creation-choice sequences for one
shared Gamma environment. The expensive Gamma reduction is performed once;
each copy receives its own independent Bernoulli coin at every potential
creation site. This is the efficient primitive needed by the double-dimer
experiment.
"""
function gamma_disordered_creation_choices(
    rng::AbstractRNG,
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real},
    copies::Integer,
)
    order = check_gamma_weights(a, b)
    copies > 0 || throw(ArgumentError("copies must be positive"))

    # A BitMatrix is important here: the number of potential creation sites is
    # sum(k^2, k=1..L) = O(L^3).  Packed decisions are roughly eight times
    # smaller than Matrix{Bool} and 64 times smaller than Float64 probabilities.
    choice_sets = [Vector{BitMatrix}(undef, order) for _ in 1:copies]
    current_a = Float64.(a)
    current_b = Float64.(b)
    next_a = similar(current_a)
    next_b = similar(current_b)
    for choices in choice_sets
        choices[order] = draw_creation_choices!(
            rng,
            falses(order, order),
            current_a,
            current_b,
        )
    end

    for level in order:-1:2
        # Reuse two full-size work buffers. Only the leading `(level-1)^2`
        # block contains the next environment. This avoids two large matrix
        # allocations at every reduction level.
        reduce_gamma_weights!(next_a, next_b, current_a, current_b, level)
        next_level = level - 1
        reduced_a = @view next_a[1:next_level, 1:next_level]
        reduced_b = @view next_b[1:next_level, 1:next_level]
        all(isfinite, reduced_a) && all(isfinite, reduced_b) ||
            throw(ArgumentError("non-finite Gamma reduction at level $level"))
        for choices in choice_sets
            # Every copy receives fresh Bernoulli draws, while all copies see
            # exactly the same reduced a/b environment.  Interleaving the RNG
            # draws does not couple them: Xoshiro supplies one independent
            # pseudorandom stream and no draw is reused between copies.
            choices[next_level] = draw_creation_choices!(
                rng,
                falses(next_level, next_level),
                reduced_a,
                reduced_b,
            )
        end
        current_a, next_a = next_a, current_a
        current_b, next_b = next_b, current_b
    end
    return choice_sets
end

"""
    gamma_disordered_creation_choice_pair(rng, a, b)

Return two conditionally independent creation-choice sequences for the same
Gamma environment.
"""
function gamma_disordered_creation_choice_pair(
    rng::AbstractRNG,
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real},
)
    choices = gamma_disordered_creation_choices(rng, a, b, 2)
    return (first=choices[1], second=choices[2])
end

"""
    uniform_creation_choices(rng, order)

Draw all independent fair creation coins for a uniformly random order-`L`
Aztec-diamond tiling.  This is the no-disorder control: every potential
creation orientation has probability one half.
"""
function uniform_creation_choices(rng::AbstractRNG, order::Integer)
    order > 0 || throw(ArgumentError("order must be positive"))
    choices = Vector{BitMatrix}(undef, order)
    for level in 1:order
        level_choices = falses(level, level)
        @inbounds for index in eachindex(level_choices)
            level_choices[index] = rand(rng, Bool)
        end
        choices[level] = level_choices
    end
    return choices
end

# ---------------------------------------------------------------------------
# Domino shuffling
# ---------------------------------------------------------------------------

function delete_and_slide(tiling::AbstractMatrix{Bool})
    side = size(tiling, 1)
    # Matrix{Bool} uses one byte per entry and is substantially faster for the
    # many scalar updates below than Julia's packed BitMatrix. Creation coins
    # remain packed bits, where memory savings matter much more.
    # Embedding an order-k diamond inside order k+1 adds one matrix cell on
    # every side.  The parity flip caused by this offset is exactly the colour
    # swap in the usual geometric description of domino shuffling.
    enlarged = fill(false, side + 2, side + 2)
    enlarged[2:(side + 1), 2:(side + 1)] .= tiling
    old_order = side ÷ 2

    # Deletion: opposite dominoes that would collide after sliding are removed
    # as a pair.  Each 2x2 block can contain at most one of the two patterns.
    @inbounds for i0 in 0:(old_order - 1), j0 in 0:(old_order - 1)
        if enlarged[2 * i0 + 1, 2 * j0 + 1] &&
           enlarged[2 * i0 + 2, 2 * j0 + 2]
            enlarged[2 * i0 + 1, 2 * j0 + 1] = false
            enlarged[2 * i0 + 2, 2 * j0 + 2] = false
        elseif enlarged[2 * i0 + 1, 2 * j0 + 2] &&
               enlarged[2 * i0 + 2, 2 * j0 + 1]
            enlarged[2 * i0 + 1, 2 * j0 + 2] = false
            enlarged[2 * i0 + 2, 2 * j0 + 1] = false
        end
    end

    # Sliding: every surviving domino moves one lattice unit in the direction
    # determined by its parity/orientation class.
    @inbounds for i0 in 0:old_order, j0 in 0:old_order
        if enlarged[2 * i0 + 2, 2 * j0 + 2]
            enlarged[2 * i0 + 1, 2 * j0 + 1] = true
            enlarged[2 * i0 + 2, 2 * j0 + 2] = false
        elseif enlarged[2 * i0 + 1, 2 * j0 + 1]
            enlarged[2 * i0 + 1, 2 * j0 + 1] = false
            enlarged[2 * i0 + 2, 2 * j0 + 2] = true
        elseif enlarged[2 * i0 + 2, 2 * j0 + 1]
            enlarged[2 * i0 + 1, 2 * j0 + 2] = true
            enlarged[2 * i0 + 2, 2 * j0 + 1] = false
        elseif enlarged[2 * i0 + 1, 2 * j0 + 2]
            enlarged[2 * i0 + 2, 2 * j0 + 1] = true
            enlarged[2 * i0 + 1, 2 * j0 + 2] = false
        end
    end
    return enlarged
end

function create_dominoes!(
    rng::AbstractRNG,
    tiling::AbstractMatrix{Bool},
    probabilities::AbstractMatrix{<:Real},
)
    order = size(tiling, 1) ÷ 2
    size(probabilities) == (order, order) ||
        throw(ArgumentError("creation-probability matrix has the wrong size"))

    # Creation is allowed only in a genuinely empty face.  The neighbouring
    # checks make this routine robust to the staggered matrix representation:
    # an empty 2x2 block alone is not sufficient near an existing domino.
    @inbounds for i0 in 0:(order - 1), j0 in 0:(order - 1)
        block_empty =
            !tiling[2 * i0 + 1, 2 * j0 + 1] &&
            !tiling[2 * i0 + 2, 2 * j0 + 1] &&
            !tiling[2 * i0 + 1, 2 * j0 + 2] &&
            !tiling[2 * i0 + 2, 2 * j0 + 2]
        block_empty || continue

        left_empty =
            j0 == 0 ||
            (!tiling[2 * i0 + 1, 2 * j0] && !tiling[2 * i0 + 2, 2 * j0])
        right_empty =
            j0 == order - 1 ||
            (!tiling[2 * i0 + 1, 2 * j0 + 3] &&
             !tiling[2 * i0 + 2, 2 * j0 + 3])
        top_empty =
            i0 == 0 ||
            (!tiling[2 * i0, 2 * j0 + 1] && !tiling[2 * i0, 2 * j0 + 2])
        bottom_empty =
            i0 == order - 1 ||
            (!tiling[2 * i0 + 3, 2 * j0 + 1] &&
             !tiling[2 * i0 + 3, 2 * j0 + 2])
        left_empty && right_empty && top_empty && bottom_empty || continue

        if rand(rng) < probabilities[i0 + 1, j0 + 1]
            tiling[2 * i0 + 1, 2 * j0 + 1] = true
            tiling[2 * i0 + 2, 2 * j0 + 2] = true
        else
            tiling[2 * i0 + 2, 2 * j0 + 1] = true
            tiling[2 * i0 + 1, 2 * j0 + 2] = true
        end
    end
    return tiling
end

function create_dominoes_from_choices!(
    tiling::AbstractMatrix{Bool},
    choices::AbstractMatrix{Bool},
)
    order = size(tiling, 1) ÷ 2
    size(choices) == (order, order) ||
        throw(ArgumentError("creation-choice matrix has the wrong size"))

    # This is the deterministic replay form of create_dominoes!.  Potential
    # sites that are not holes simply ignore their pre-drawn bit, which is
    # distributionally the same as drawing a coin only after finding a hole.
    @inbounds for i0 in 0:(order - 1), j0 in 0:(order - 1)
        block_empty =
            !tiling[2 * i0 + 1, 2 * j0 + 1] &&
            !tiling[2 * i0 + 2, 2 * j0 + 1] &&
            !tiling[2 * i0 + 1, 2 * j0 + 2] &&
            !tiling[2 * i0 + 2, 2 * j0 + 2]
        block_empty || continue

        left_empty =
            j0 == 0 ||
            (!tiling[2 * i0 + 1, 2 * j0] && !tiling[2 * i0 + 2, 2 * j0])
        right_empty =
            j0 == order - 1 ||
            (!tiling[2 * i0 + 1, 2 * j0 + 3] &&
             !tiling[2 * i0 + 2, 2 * j0 + 3])
        top_empty =
            i0 == 0 ||
            (!tiling[2 * i0, 2 * j0 + 1] && !tiling[2 * i0, 2 * j0 + 2])
        bottom_empty =
            i0 == order - 1 ||
            (!tiling[2 * i0 + 3, 2 * j0 + 1] &&
             !tiling[2 * i0 + 3, 2 * j0 + 2])
        left_empty && right_empty && top_empty && bottom_empty || continue

        if choices[i0 + 1, j0 + 1]
            tiling[2 * i0 + 1, 2 * j0 + 1] = true
            tiling[2 * i0 + 2, 2 * j0 + 2] = true
        else
            tiling[2 * i0 + 2, 2 * j0 + 1] = true
            tiling[2 * i0 + 1, 2 * j0 + 2] = true
        end
    end
    return tiling
end

"""
    sample_tiling(rng, probabilities)

Run the deletion, sliding, and creation steps and return the same binary
`2order × 2order` matrix encoding used by the supplied Julia/Mathematica code.
Each `true` entry represents one domino.
"""
function sample_tiling(
    rng::AbstractRNG,
    probabilities::AbstractVector{<:AbstractMatrix},
)
    order = length(probabilities)
    order > 0 || throw(ArgumentError("at least one probability matrix is required"))

    size(probabilities[1]) == (1, 1) ||
        throw(ArgumentError("probability level 1 must have size (1, 1)"))
    first_probability = probabilities[1][1, 1]
    isfinite(first_probability) && 0 <= first_probability <= 1 ||
        throw(ArgumentError("probability level 1 must lie in [0, 1]"))

    tiling = if rand(rng) < first_probability
        Bool[1 0; 0 1]
    else
        Bool[0 1; 1 0]
    end

    for level in 2:order
        tiling = delete_and_slide(tiling)
        create_dominoes!(rng, tiling, probabilities[level])
    end
    return tiling
end

"""
    sample_tiling_from_choices(choices)

Run domino shuffling using independent creation choices which were drawn in
advance by `gamma_disordered_creation_choices`.
"""
function sample_tiling_from_choices(
    choices::AbstractVector{<:AbstractMatrix{Bool}},
)
    order = length(choices)
    order > 0 || throw(ArgumentError("at least one choice matrix is required"))
    size(choices[1]) == (1, 1) ||
        throw(ArgumentError("choice level 1 must have size (1, 1)"))
    tiling = if choices[1][1, 1]
        Bool[1 0; 0 1]
    else
        Bool[0 1; 1 0]
    end
    for level in 2:order
        tiling = delete_and_slide(tiling)
        create_dominoes_from_choices!(tiling, choices[level])
    end
    return tiling
end


"""
    sample_gamma_center_height(seed, order; alpha=0.2, beta=0.25)

Generate one independent Gamma-disordered Aztec-diamond environment,
sample one dimer configuration conditional on that environment, and
return the height at the central face.
"""
function sample_gamma_center_height(
    seed::Integer,
    order::Integer;
    alpha::Real=0.2,
    beta::Real=0.25,
)
    order > 0 || throw(ArgumentError("order must be positive"))
    alpha > 0 || throw(ArgumentError("alpha must be positive"))
    beta > 0 || throw(ArgumentError("beta must be positive"))

    rng = Xoshiro(UInt64(seed))
    weights = gamma_disordered_weights(
        rng,
        order;
        alpha=alpha,
        beta=beta,
    )
    choices = gamma_disordered_creation_choices(
        rng,
        weights.a,
        weights.b,
    )
    tiling = sample_tiling_from_choices(choices)

    count(tiling) == order * (order + 1) ||
        throw(ErrorException(
            "invalid domino count for order=$order seed=$seed",
        ))

    return center_height(tiling)
end

"""
    sample_gamma_center_height_pair(seed, order; alpha=0.2, beta=0.25)

Generate one Gamma environment and two independent domino tilings conditional
on that same environment. Return both central heights and their difference.
This is one observation from the double-dimer experiment.
"""
function sample_gamma_center_height_pair(
    seed::Integer,
    order::Integer;
    alpha::Real=0.2,
    beta::Real=0.25,
)
    order > 0 || throw(ArgumentError("order must be positive"))
    alpha > 0 || throw(ArgumentError("alpha must be positive"))
    beta > 0 || throw(ArgumentError("beta must be positive"))

    # The RNG stream has three non-overlapping pieces: the shared environment,
    # creation coins for copy 1, and creation coins for copy 2.  Only the first
    # piece is shared.  Thus the tilings are independent conditional on a,b.
    rng = Xoshiro(UInt64(seed))
    weights = gamma_disordered_weights(
        rng,
        order;
        alpha=alpha,
        beta=beta,
    )
    choices = gamma_disordered_creation_choice_pair(
        rng,
        weights.a,
        weights.b,
    )

    first_tiling = sample_tiling_from_choices(choices.first)
    second_tiling = sample_tiling_from_choices(choices.second)
    expected_dominoes = order * (order + 1)
    count(first_tiling) == expected_dominoes ||
        throw(ErrorException("invalid first tiling for order=$order seed=$seed"))
    count(second_tiling) == expected_dominoes ||
        throw(ErrorException("invalid second tiling for order=$order seed=$seed"))

    first_height = center_height(first_tiling)
    second_height = center_height(second_tiling)
    return (
        height_1=first_height,
        height_2=second_height,
        difference=first_height - second_height,
    )
end


"""
    sample_gamma_spatial_increment_pair(seed, order, separations; alpha, beta)

Generate one Gamma environment and two conditionally independent tilings.
For each requested face-column separation, return the height increment between
two faces placed symmetrically around the centre row.  The two replicas make
the connected and disorder contributions directly observable:

`Var(delta_1-delta_2)/2` is conditional tiling noise, while
`Cov(delta_1,delta_2)` is disorder-induced fluctuation.
"""
function sample_gamma_spatial_increment_pair(
    seed::Integer,
    order::Integer,
    separations::AbstractVector{<:Integer};
    alpha::Real=0.2,
    beta::Real=0.25,
)
    order > 0 || throw(ArgumentError("order must be positive"))
    alpha > 0 || throw(ArgumentError("alpha must be positive"))
    beta > 0 || throw(ArgumentError("beta must be positive"))
    all(separation -> 0 < separation < order, separations) ||
        throw(ArgumentError("every separation must lie between 1 and order-1"))

    rng = Xoshiro(UInt64(seed))
    weights = gamma_disordered_weights(rng, order; alpha=alpha, beta=beta)
    choices = gamma_disordered_creation_choice_pair(rng, weights.a, weights.b)
    first_tiling = sample_tiling_from_choices(choices.first)
    second_tiling = sample_tiling_from_choices(choices.second)
    expected_dominoes = order * (order + 1)
    count(first_tiling) == expected_dominoes ||
        throw(ErrorException("invalid first spatial tiling for order=$order seed=$seed"))
    count(second_tiling) == expected_dominoes ||
        throw(ErrorException("invalid second spatial tiling for order=$order seed=$seed"))

    return [
        begin
            columns = symmetric_face_columns(order, separation)
            first_increment = symmetric_height_increment(first_tiling, separation)
            second_increment = symmetric_height_increment(second_tiling, separation)
            (
                separation=Int(separation),
                left_column=columns.left,
                right_column=columns.right,
                increment_1=first_increment,
                increment_2=second_increment,
                difference=first_increment - second_increment,
            )
        end
        for separation in separations
    ]
end

"""
    sample_uniform_spatial_increment_pair(seed, order, separations)

Generate two independent uniform Aztec-diamond tilings and measure the same
spatial increments as the Gamma experiment.  There is no shared random
environment, so their cross-replica covariance should be zero; their marginal
and half-difference variances provide the ordinary-log control.
"""
function sample_uniform_spatial_increment_pair(
    seed::Integer,
    order::Integer,
    separations::AbstractVector{<:Integer},
)
    order > 0 || throw(ArgumentError("order must be positive"))
    all(separation -> 0 < separation < order, separations) ||
        throw(ArgumentError("every separation must lie between 1 and order-1"))
    rng = Xoshiro(UInt64(seed))
    first_tiling = sample_tiling_from_choices(uniform_creation_choices(rng, order))
    second_tiling = sample_tiling_from_choices(uniform_creation_choices(rng, order))
    expected_dominoes = order * (order + 1)
    count(first_tiling) == expected_dominoes ||
        throw(ErrorException("invalid first uniform tiling for order=$order seed=$seed"))
    count(second_tiling) == expected_dominoes ||
        throw(ErrorException("invalid second uniform tiling for order=$order seed=$seed"))

    return [
        begin
            columns = symmetric_face_columns(order, separation)
            first_increment = symmetric_height_increment(first_tiling, separation)
            second_increment = symmetric_height_increment(second_tiling, separation)
            (
                separation=Int(separation),
                left_column=columns.left,
                right_column=columns.right,
                increment_1=first_increment,
                increment_2=second_increment,
                difference=first_increment - second_increment,
            )
        end
        for separation in separations
    ]
end


# ---------------------------------------------------------------------------
# Height observable
# ---------------------------------------------------------------------------

"""
    height_function(tiling)

Compute the staggered face-height table used by this study. For an order-`L` tiling,
the result has size `(2L + 1) × (L + 1)`. Even and odd rows represent
staggered face locations. Crossing an occupied dimer edge changes the height
by `+3`; crossing an unoccupied edge changes it by `-1`.
"""
function height_function(tiling::AbstractMatrix{Bool})
    side1, side2 = size(tiling)
    side1 == side2 || throw(ArgumentError("tiling matrix must be square"))
    iseven(side1) || throw(ArgumentError("tiling matrix side length must be even"))
    order = side1 ÷ 2
    heights = zeros(Int, side1 + 1, order + 1)

    for column in 1:(order + 1)
        heights[1, column] = 2 * column - 2
    end
    for row in 2:(side1 + 1), column in 1:order
        heights[row, column] =
            heights[row - 1, column] +
            (tiling[row - 1, 2 * column - 1] ? 3 : -1)
    end
    return heights
end

"""
    face_height(tiling, row, column)

Compute one entry of `height_function(tiling)` directly.  This costs only one
vertical path and avoids allocating the full staggered height table.
"""
function face_height(
    tiling::AbstractMatrix{Bool},
    row::Integer,
    column::Integer,
)
    side1, side2 = size(tiling)
    side1 == side2 || throw(ArgumentError("tiling matrix must be square"))
    iseven(side1) || throw(ArgumentError("tiling matrix side length must be even"))
    order = side1 ÷ 2
    1 <= row <= side1 + 1 || throw(ArgumentError("face row is out of range"))
    1 <= column <= order || throw(ArgumentError("face column is out of range"))

    height = 2 * column - 2
    edge_column = 2 * column - 1
    @inbounds for edge_row in 1:(row - 1)
        height += tiling[edge_row, edge_column] ? 3 : -1
    end
    return height
end

"""
    center_face_index(order)

Return the row and column of the face closest to the geometric center in the
staggered height table.
"""
function center_face_index(order::Integer)
    order > 0 || throw(ArgumentError("order must be positive"))
    return (row=order + 1, column=fld(order, 2) + 1)
end

"""
    center_height(tiling)

Compute only the central face height, without allocating the full height
table. This is exactly `height_function(tiling)[center_face_index(order)...]`.
"""
function center_height(tiling::AbstractMatrix{Bool})
    side1, side2 = size(tiling)
    side1 == side2 || throw(ArgumentError("tiling matrix must be square"))
    iseven(side1) || throw(ArgumentError("tiling matrix side length must be even"))
    order = side1 ÷ 2
    center = center_face_index(order)
    return face_height(tiling, center.row, center.column)
end

"""
    symmetric_face_columns(order, separation)

Return two valid face columns separated by `separation` and centred as closely
as parity permits around the geometric middle of the Aztec diamond.
"""
function symmetric_face_columns(order::Integer, separation::Integer)
    order > 1 || throw(ArgumentError("order must exceed one"))
    0 < separation < order ||
        throw(ArgumentError("separation must lie between 1 and order-1"))
    left = fld(order - separation, 2) + 1
    right = left + separation
    1 <= left < right <= order || error("internal symmetric-column error")
    return (left=left, right=right)
end

"""
    symmetric_height_increment(tiling, separation)

Height at the right member minus height at the left member of a symmetric
pair of faces on the central row.  Only two vertical paths are evaluated.
"""
function symmetric_height_increment(
    tiling::AbstractMatrix{Bool},
    separation::Integer,
)
    side1, side2 = size(tiling)
    side1 == side2 || throw(ArgumentError("tiling matrix must be square"))
    iseven(side1) || throw(ArgumentError("tiling matrix side length must be even"))
    order = side1 ÷ 2
    columns = symmetric_face_columns(order, separation)
    row = order + 1
    return face_height(tiling, row, columns.right) -
           face_height(tiling, row, columns.left)
end

# ---------------------------------------------------------------------------
# Validation and output helpers
# ---------------------------------------------------------------------------

function covered_cells(tiling::AbstractMatrix{Bool})
    side = size(tiling, 1)
    cells = Set{Tuple{Int,Int}}()
    overlap_count = 0

    # A true matrix entry denotes a domino centre/parity class, not a board
    # square.  Convert it to the two lattice cells covered by that domino so
    # overlaps and holes can be checked geometrically.
    for i in axes(tiling, 1), j in axes(tiling, 2)
        tiling[i, j] || continue
        x = j - i
        y = side + 1 - i - j
        domino_cells = if isodd(i) == isodd(j)
            ((x - 1, y), (x + 1, y))
        else
            ((x, y - 1), (x, y + 1))
        end
        for cell in domino_cells
            cell in cells && (overlap_count += 1)
            push!(cells, cell)
        end
    end
    return cells, overlap_count
end

"""
    validate_tiling(tiling)

Check the matrix size, domino count, absence of overlaps, and exact coverage
of the order-`n` Aztec diamond. Return a named tuple of validation statistics.
"""
function validate_tiling(tiling::AbstractMatrix{Bool})
    side1, side2 = size(tiling)
    side1 == side2 || throw(ArgumentError("tiling matrix must be square"))
    iseven(side1) || throw(ArgumentError("tiling matrix side length must be even"))
    order = side1 ÷ 2
    expected_dominoes = order * (order + 1)
    dominoes = count(tiling)
    cells, overlap_count = covered_cells(tiling)

    # The target cells are the odd-parity lattice points inside the L1 ball
    # representing the order-L Aztec diamond.
    target_cells = Set{Tuple{Int,Int}}()
    for x in (-(2 * order - 1)):2:(2 * order - 1)
        for y in (-(2 * order - 1)):2:(2 * order - 1)
            abs(x) + abs(y) <= 2 * order && push!(target_cells, (x, y))
        end
    end

    missing_cells = length(setdiff(target_cells, cells))
    outside_cells = length(setdiff(cells, target_cells))
    valid =
        dominoes == expected_dominoes &&
        overlap_count == 0 &&
        missing_cells == 0 &&
        outside_cells == 0

    return (
        valid=valid,
        order=order,
        matrix_side=side1,
        dominoes=dominoes,
        expected_dominoes=expected_dominoes,
        covered_cells=length(cells),
        expected_cells=2 * expected_dominoes,
        overlaps=overlap_count,
        missing_cells=missing_cells,
        outside_cells=outside_cells,
    )
end

"""
    orientation_counts(tiling)

Count the four parity/orientation classes used by the legacy Mathematica
renderer: green, blue, yellow, and red.
"""
function orientation_counts(tiling::AbstractMatrix{Bool})
    counts = Dict(:green => 0, :blue => 0, :yellow => 0, :red => 0)
    for i in axes(tiling, 1), j in axes(tiling, 2)
        tiling[i, j] || continue
        key = if isodd(i) && isodd(j)
            :green
        elseif isodd(i) && iseven(j)
            :blue
        elseif iseven(i) && iseven(j)
            :yellow
        else
            :red
        end
        counts[key] += 1
    end
    return (
        green=counts[:green],
        blue=counts[:blue],
        yellow=counts[:yellow],
        red=counts[:red],
    )
end

"""
    write_table(path, matrix; delimiter=" ")

Write a numeric or Boolean matrix in a simple, lossless text format. Floating
point values use 17 significant digits, enough to round-trip a `Float64`.
"""
function write_table(
    path::AbstractString,
    matrix::AbstractMatrix;
    delimiter::AbstractString=" ",
)
    mkpath(dirname(path))
    open(path, "w") do io
        for i in axes(matrix, 1)
            for j in axes(matrix, 2)
                j == first(axes(matrix, 2)) || print(io, delimiter)
                value = matrix[i, j]
                if value isa Bool
                    print(io, value ? 1 : 0)
                elseif value isa AbstractFloat
                    @printf(io, "%.17g", value)
                else
                    print(io, value)
                end
            end
            println(io)
        end
    end
    return path
end

function rotate_point(x::Real, y::Real)
    scale = inv(sqrt(2.0))
    return ((x + y) * scale, (x - y) * scale)
end

"""
    write_svg(path, tiling; image_size=1800)

Render the tiling with the same coordinates, `-45°` rotation, and four colour
classes as the Mathematica `AztecPrinter` supplied with the old Julia code.
"""
function write_svg(
    path::AbstractString,
    tiling::AbstractMatrix{Bool};
    image_size::Integer=1800,
)
    validation = validate_tiling(tiling)
    validation.valid || throw(ArgumentError("refusing to draw an invalid tiling"))
    side = size(tiling, 1)
    margin = 2.0
    extent = (side + 2) / sqrt(2.0) + margin
    min_coordinate = -extent
    view_side = 2 * extent

    colors = Dict(
        :green => "#00ff00",
        :blue => "#0000ff",
        :yellow => "#ffff00",
        :red => "#ff0000",
    )

    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        println(
            io,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" ",
            "viewBox=\"$min_coordinate $min_coordinate $view_side $view_side\" ",
            "width=\"$image_size\" height=\"$image_size\">",
        )
        println(io, "<rect x=\"$min_coordinate\" y=\"$min_coordinate\" ",
                "width=\"$view_side\" height=\"$view_side\" fill=\"white\"/>")

        for i in axes(tiling, 1), j in axes(tiling, 2)
            tiling[i, j] || continue
            center_x = j - i
            center_y = side + 1 - i - j
            horizontal = isodd(i) == isodd(j)
            half_width = horizontal ? 2.0 : 1.0
            half_height = horizontal ? 1.0 : 2.0
            corners = (
                (center_x - half_width, center_y - half_height),
                (center_x + half_width, center_y - half_height),
                (center_x + half_width, center_y + half_height),
                (center_x - half_width, center_y + half_height),
            )
            rotated = map(corners) do (x, y)
                x_rotated, y_rotated = rotate_point(x, y)
                (x_rotated, -y_rotated)
            end
            points = join(("$x,$y" for (x, y) in rotated), " ")
            key = if isodd(i) && isodd(j)
                :green
            elseif isodd(i) && iseven(j)
                :blue
            elseif iseven(i) && iseven(j)
                :yellow
            else
                :red
            end
            println(io, "<polygon points=\"$points\" fill=\"$(colors[key])\"/>")
        end
        println(io, "</svg>")
    end
    return path
end

# Square-grid spanning-tree / Temperley-dimer extension.
include("SquareGrid.jl")
export SquareGrid

end
