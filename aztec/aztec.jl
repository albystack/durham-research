#!/usr/bin/env julia
# Weighted domino tilings of an Aztec diamond of order L.
# Each environment has a ~ Gamma(0.2, 1), b ~ Gamma(0.25, 1), and unit weights
# on the other two edge families. For each fixed environment, draw two tilings
# independently. Measure height differences at offsets L/8 and L/4 along each
# coordinate of the diagonal; the physical distances are sqrt(2) times these.
#
# Run: julia aztec.jl --size 8 --environments 2 --output aztec-example
# Requires Julia 1.10+; no extra packages. The output folder must be new.
# Defaults run two small pairs. The figures use the completed campaign data.
#
# Reading order: weights -> reduction -> shuffling -> heights -> batch runner.
# Source: albystack/durham-research, 27b4bff. Sampling steps and seeds are unchanged.
# OpenAI Codex assisted with implementation, testing, analysis and editing.

module AztecDiamond

using Random
using Printf

export gamma_disordered_weights, sample_tiling_from_choices,
       symmetric_height_increment, validate_tiling, write_table


# 1. Random edge weights

"""Draw one Gamma variable using rejection sampling."""
function rand_gamma(rng::AbstractRNG, shape::Float64, scale::Float64)
    # For shape < 1, use Gamma(shape + 1) * U^(1/shape).
    if shape < 1
        uniform_draw = rand(rng)
        while uniform_draw == 0
            uniform_draw = rand(rng)
        end
        return rand_gamma(rng, shape + 1, scale) * uniform_draw^(1 / shape)
    end

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

"""Draw an array of independent Gamma weights; mean = shape * scale."""
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

"""Draw the two independent edge-weight families of one environment."""
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

"""Require two square arrays of positive, finite weights."""
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

# 2. Reduce an order-L environment to order L-1

"""Allocate the smaller arrays and apply one weight reduction."""
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

"""Apply the Gamma reduction, reusing the supplied arrays."""
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

    # These are the two Gamma weight recurrences, applied in paper coordinates.
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

# 3. Delete, slide, and create dominoes

# The 2L-by-2L array records domino centres, not board squares. Row and column
# parity encode orientation; covered_cells converts these entries to board cells.

"""Enlarge the diamond, remove colliding pairs, and slide the survivors."""
function delete_and_slide(tiling::AbstractMatrix{Bool})
    side = size(tiling, 1)
    # Adding one matrix entry on each side changes the parity of old entries.
    enlarged = fill(false, side + 2, side + 2)
    enlarged[2:(side + 1), 2:(side + 1)] .= tiling
    old_order = side ÷ 2

    # Remove opposite pairs that would collide.
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

    # Each survivor moves to the opposite corner of its active cell.
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

"""Fill each empty cell using its already-drawn coin."""
function create_dominoes_from_choices!(
    tiling::AbstractMatrix{Bool},
    choices::AbstractMatrix{Bool},
)
    order = size(tiling, 1) ÷ 2
    size(choices) == (order, order) ||
        throw(ArgumentError("creation-choice matrix has the wrong size"))

    # Ignore coins at cells that are not empty; adjacent dominoes also block creation.
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

"""Grow a tiling from order 1. Each true entry of the final matrix is a domino."""
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

# 4. Height differences and checks

"""Integrate the staggered height table using +3 at occupied entries and -1 elsewhere."""
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

"""Evaluate one height by summing along a vertical path in the matrix."""
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

"""Choose two face columns centred on the diamond."""
function symmetric_face_columns(order::Integer, separation::Integer)
    order > 1 || throw(ArgumentError("order must exceed one"))
    0 < separation < order ||
        throw(ArgumentError("separation must lie between 1 and order-1"))
    left = fld(order - separation, 2) + 1
    right = left + separation
    1 <= left < right <= order || error("internal symmetric-column error")
    return (left=left, right=right)
end

"""Return right height minus left height on the central row."""
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

"""Convert each occupied matrix entry to the two board cells it covers."""
function covered_cells(tiling::AbstractMatrix{Bool})
    side = size(tiling, 1)
    cells = Set{Tuple{Int,Int}}()
    overlap_count = 0

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

"""Check the domino count, overlaps, and exact coverage of the diamond."""
function validate_tiling(tiling::AbstractMatrix{Bool})
    side1, side2 = size(tiling)
    side1 == side2 || throw(ArgumentError("tiling matrix must be square"))
    iseven(side1) || throw(ArgumentError("tiling matrix side length must be even"))
    order = side1 ÷ 2
    expected_dominoes = order * (order + 1)
    dominoes = count(tiling)
    cells, overlap_count = covered_cells(tiling)

    # Board cells have odd coordinates and satisfy |x| + |y| <= 2L.
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

"""Save a matrix as text; 17 significant digits preserve Float64 weights."""
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

end # module AztecDiamond

using .AztecDiamond
using Random, SHA, Printf

# 5. Creation coins for two tilings in the same environment

const ROBUST_GAMMA_VERSION = "gamma-interval-coins-v1"

"""Bound b/(a+b) from below and above, allowing for floating-point rounding."""
function coin_bounds(a::Float64, b::Float64)
    isfinite(a) && isfinite(b) && a > 0 && b > 0 || error("positive finite weights required")
    s = a + b
    isfinite(s) || return (0.0, 1.0)
    slo, shi = prevfloat(s), nextfloat(s)
    lo = max(0.0, prevfloat(b / shi))
    hi = slo > 0 ? min(1.0, nextfloat(b / slo)) : 1.0
    return (lo, hi)
end

"""Add random bits until the interval containing U lies wholly above or below n/d."""
function rational_coin_prefix(rng, n::BigInt, d::BigInt, prefix::BigInt, bits::Int)
    0 <= n <= d && d > 0 && 0 <= prefix < (big(1) << bits) || error("invalid rational coin")
    while true
        threshold = n << bits
        (prefix + 1) * d <= threshold && return true
        prefix * d >= threshold && return false
        prefix = (prefix << 64) + BigInt(rand(rng, UInt64))
        bits += 64
    end
end

"""Draw a coin with probability b/(a+b); use exact fractions when bounds overlap."""
function robust_gamma_coin(rng, a::Float64, b::Float64, bounds=coin_bounds(a, b))
    # The first 53 random bits locate U inside an interval of width 2^-53.
    prefix = rand(rng, UInt64) >> 11
    lower = ldexp(Float64(prefix), -53)
    upper = ldexp(Float64(prefix + 1), -53)
    lo, hi = bounds
    upper <= lo && return (true, false)
    lower >= hi && return (false, false)
    ar, br = Rational{BigInt}(a), Rational{BigInt}(b)
    p = br / (ar + br)
    return (rational_coin_prefix(rng, numerator(p), denominator(p), BigInt(prefix), 53), true)
end

"""Reduce the shared weights once and draw separate creation coins for both replicas."""
function robust_gamma_choice_pair(rng, a, b)
    order = AztecDiamond.check_gamma_weights(a, b)
    sets = [Vector{BitMatrix}(undef, order) for _ in 1:2]
    # The reduction and tiling matrix use opposite row directions: reflect both ways.
    ca, cb = reverse(Float64.(a), dims=1), reverse(Float64.(b), dims=1)
    na, nb = similar(ca), similar(cb)
    fallbacks = 0
    for level in order:-1:1
        first, second = falses(level, level), falses(level, level)
        for j in 1:level, i in 1:level
            bounds = coin_bounds(ca[i, j], cb[i, j])
            first[i, j], f1 = robust_gamma_coin(rng, ca[i, j], cb[i, j], bounds)
            second[i, j], f2 = robust_gamma_coin(rng, ca[i, j], cb[i, j], bounds)
            fallbacks += f1 + f2
        end
        sets[1][level] = reverse(first, dims=1)
        sets[2][level] = reverse(second, dims=1)
        if level > 1
            AztecDiamond.reduce_gamma_weights!(na, nb, ca, cb, level)
            ca, na = na, ca
            cb, nb = nb, cb
        end
    end
    return (first=sets[1], second=sets[2], fallbacks=fallbacks)
end

# Coins are exact for the stored weights. The preceding Float64 reductions
# still have rounding error; small-grid checks do not prove universal exactness.


# 6. Reproducible batches and saved results

filehash(path) = bytes2hex(sha256(read(path)))

"""Derive a repeatable seed for an environment or its pair of tilings."""
function derive_seed(master, model, n, index, role)
    digest = sha256("matched-shuffle-20261004:$master:aztec:$model:$n:$index:$role")
    foldl((x, b) -> (x << 8) | UInt64(b), digest[1:8]; init=UInt64(0))
end

"""For each environment, draw two tilings and save their height increments."""
function main(args)
    iseven(length(args)) || error("arguments must be --name value pairs")
    options = Dict(args[i] => args[i+1] for i in 1:2:length(args))
    Set(keys(options)) == Set([
        "--output", "--size", "--model", "--environments", "--seed", "--plan-sha256", "--validation",
    ]) || error("wrong options")
    n = parse(Int, options["--size"])
    model = options["--model"]
    count = parse(Int, options["--environments"])
    master = parse(Int, options["--seed"])
    validation = options["--validation"]
    planhash = options["--plan-sha256"]
    n in (8, 16, 24, 32, 48, 64, 96, 128) && 1 <= count <= 512 || error("outside larger-batch limits")
    model in ("aztec_gamma", "aztec_mean", "uniform") || error("unknown model")
    validation in ("true", "false") || error("invalid validation flag")
    output = options["--output"]
    ispath(output) && error("fresh output directory required")
    mkpath(output)
    write(
        joinpath(output, "manifest.json"),
        "{\"domain\":\"aztec\",\"size\":$n,\"model\":\"$model\",\"pairs\":$count,\"master_seed\":$master,\"validation_only\":$validation,\"plan_sha256\":\"$planhash\",\"sampler_version\":\"$ROBUST_GAMMA_VERSION\",\"julia\":\"$VERSION\"}\n",
    )
    started = time()
    write(joinpath(output, "status.json"), "{\"state\":\"STARTED\",\"completed\":0}\n")
    temporary = joinpath(output, "measurements.csv.tmp")

    open(temporary, "w") do table
        println(table, "domain,size,model,law_id,environment_id,probe,axis,r,distance,height_scale,value1,value2,environment_seed,coin_seed")
        for index in 1:count
            envseed = derive_seed(master, model, n, index, "environment")
            coinseed = derive_seed(master, model, n, index, "paired-coins")
            envseed != coinseed || error("seed collision")
            if model == "aztec_gamma"
                w = gamma_disordered_weights(Xoshiro(envseed), n; alpha=.2, beta=.25)
                a, b = w.a, w.b
            else
                a, b = model == "aztec_mean" ? (fill(.2, n, n), fill(.25, n, n)) : (ones(n, n), ones(n, n))
            end

            # Both tilings share a and b, but use separate random creation coins.
            choices = robust_gamma_choice_pair(Xoshiro(coinseed), a, b)
            first, second = sample_tiling_from_choices(choices.first), sample_tiling_from_choices(choices.second)
            validate_tiling(first).valid && validate_tiling(second).valid || error("invalid tiling")
            case = joinpath(output, "cases", @sprintf("env-%05d", index))
            mkpath(case)
            for (name, value) in zip(("a.txt", "b.txt", "sample1.txt", "sample2.txt"), (a, b, first, second))
                write_table(joinpath(case, name), value)
            end

            id = "aztec-$model-n$n-env$index"
            rows = IOBuffer()
            for rho in (1/8, 1/4)
                r = round(Int, n * rho)
                law = "\"$model:alpha=0.2:beta=0.25:scale=1:a-north,b-east\""
                label = @sprintf("diagonal:rho=%.12g", rho)
                println(rows, join((
                    "aztec", n, model, law, id, label, "diagonal", r, sqrt(2) * r, 1,
                    symmetric_height_increment(first, r), symmetric_height_increment(second, r), envseed, coinseed,
                ), ','))
            end
            rowtext = String(take!(rows))
            write(table, rowtext)
            write(joinpath(case, "rows.csv"), rowtext)
            write(
                joinpath(case, "sha256.txt"),
                join(["$(filehash(joinpath(case, f))) $f" for f in ("a.txt", "b.txt", "sample1.txt", "sample2.txt", "rows.csv")], "\n") * "\n",
            )
            if index % 64 == 0
                flush(table)
                write(joinpath(output, "status.json"), "{\"state\":\"STARTED\",\"completed\":$index}\n")
                println("aztec $model n=$n completed=$index/$count")
                flush(stdout)
            end
        end
    end

    # Only completed batches receive the final filename and COMPLETE status.
    mv(temporary, joinpath(output, "measurements.csv"))
    write(
        joinpath(output, "status.json"),
        "{\"state\":\"COMPLETE\",\"completed\":$count,\"seconds\":$(time()-started),\"measurements_sha256\":\"$(filehash(joinpath(output,"measurements.csv")))\",\"arithmetic\":\"Float64 reductions with robust creation coins; rounding audited separately\"}\n",
    )
end

"""Use a small demonstration unless explicit command-line options are supplied."""
function standalone_arguments(args)
    if "--help" in args || "-h" in args
        println("Usage: julia aztec.jl --size 8 --environments 2 --output aztec-example")
        println("Optional: --model aztec_gamma|aztec_mean|uniform --seed INTEGER")
        println("Sizes: 8,16,24,32,48,64,96,128. At most 512 pairs per invocation.")
        println("The output directory must not already exist.")
        return nothing
    end
    iseven(length(args)) || error("arguments must be --name value pairs")
    options = Dict(
        "--size" => "8", "--environments" => "2", "--output" => "aztec-example",
        "--model" => "aztec_gamma", "--seed" => "2026101302",
        "--plan-sha256" => "standalone-demonstration", "--validation" => "true",
    )
    for i in 1:2:length(args)
        haskey(options, args[i]) || error("unknown option: $(args[i])")
        options[args[i]] = args[i+1]
    end
    return reduce(vcat, [[key, value] for (key, value) in sort(collect(options))])
end

if abspath(PROGRAM_FILE) == @__FILE__
    options = standalone_arguments(ARGS)
    options === nothing || main(options)
end
