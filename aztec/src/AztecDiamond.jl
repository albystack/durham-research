module AztecDiamond

using Random
using Printf

export random_uniform_weights,
       random_gamma_weights,
       gamma_disordered_weights,
       creation_probabilities,
       gamma_disordered_probabilities,
       sample_tiling,
       validate_tiling,
       orientation_counts,
       write_table,
       write_svg

"""
    random_uniform_weights(rng, order)

Return the `2order × 2order` i.i.d. `Uniform(0,1)` positive weight table
used by the random-weight example at the end of Sunil Chhita's 2022
`simulatorfinal.jl`.
"""
function random_uniform_weights(rng::AbstractRNG, order::Integer)
    order > 0 || throw(ArgumentError("order must be positive"))
    return rand(rng, Float64, 2 * order, 2 * order)
end

function rand_gamma(rng::AbstractRNG, shape::Float64, scale::Float64)
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
    size(a) == size(b) || throw(ArgumentError("a and b must have the same size"))
    rows, cols = size(a)
    rows == cols || throw(ArgumentError("a and b must be square"))
    rows > 0 || throw(ArgumentError("weight tables cannot be empty"))
    all(isfinite, a) && all(isfinite, b) ||
        throw(ArgumentError("Gamma weights must be finite"))
    all(>(0), a) && all(>(0), b) ||
        throw(ArgumentError("Gamma weights must be strictly positive"))

    order = rows
    a_levels = Vector{Matrix{Float64}}(undef, order)
    b_levels = Vector{Matrix{Float64}}(undef, order)
    a_levels[order] = Float64.(a)
    b_levels[order] = Float64.(b)

    for level in order:-1:2
        current_a = a_levels[level]
        current_b = b_levels[level]
        reduced_a = Matrix{Float64}(undef, level - 1, level - 1)
        reduced_b = Matrix{Float64}(undef, level - 1, level - 1)

        for i in 1:(level - 1), j in 1:(level - 1)
            reduced_a[i, j] =
                current_a[i, j] / (current_a[i, j] + current_b[i, j]) *
                (current_a[i + 1, j] + current_b[i + 1, j])
            reduced_b[i, j] =
                current_b[i, j + 1] /
                (current_a[i, j + 1] + current_b[i, j + 1]) *
                (current_a[i + 1, j + 1] + current_b[i + 1, j + 1])
        end
        all(isfinite, reduced_a) && all(isfinite, reduced_b) ||
            throw(ArgumentError("non-finite Gamma reduction at level $level"))
        a_levels[level - 1] = reduced_a
        b_levels[level - 1] = reduced_b
    end

    probabilities = Vector{Matrix{Float64}}(undef, order)
    for level in 1:order
        denominator = a_levels[level] + b_levels[level]
        probability = b_levels[level] ./ denominator
        all(p -> isfinite(p) && 0 <= p <= 1, probability) ||
            throw(ArgumentError("invalid Gamma creation probability at level $level"))
        probabilities[level] = probability
    end
    return probabilities
end

function delete_and_slide(tiling::AbstractMatrix{Bool})
    side = size(tiling, 1)
    enlarged = falses(side + 2, side + 2)
    enlarged[2:(side + 1), 2:(side + 1)] .= tiling
    old_order = side ÷ 2

    for i0 in 0:(old_order - 1), j0 in 0:(old_order - 1)
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

    for i0 in 0:old_order, j0 in 0:old_order
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

    for i0 in 0:(order - 1), j0 in 0:(order - 1)
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

    tiling = if rand(rng) < probabilities[1][1, 1]
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

Count the four parity/orientation classes used by Sunil's Mathematica
`AztecPrinter`: green, blue, yellow, and red.
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

end
