"""
Square-grid weighted spanning trees and their Temperley dimer matchings.

The finite graph has interior vertices `(x,y)` with `max(abs(x),abs(y)) < L`.
Every step from an interior vertex to the square boundary is wired to one root,
while retaining its geometric direction.  Wilson's algorithm samples a rooted
arborescence using either directed site weights or shared undirected
conductances.  The complementary dual tree then gives the exact generalized
Temperley perfect matching.

The spatial observable is a dimer-height increment across a horizontal dual
cut through the centre.  Along this cut all white nodes have degree four and
the deterministic reference-flow terms cancel, so the returned integer is the
full height increment in the convention fixed below, not merely a proxy.
"""
module SquareGrid

using Random

export AbstractSquareGridEnvironment,
       BaselineEnvironment,
       DirectedSiteIIDEnvironment,
       UndirectedConductanceEnvironment,
       TemporalIIDEnvironment,
       FixedEnvironment,
       TreeSample,
       TemperleyMatching,
       NORTH,
       EAST,
       SOUTH,
       WEST,
       edge_weight,
       sample_full_tree,
       validate_tree,
       build_temperley_matching,
       validate_matching,
       symmetric_probe_vertices,
       spatial_height_increment,
       dimer_height_increment_along_path,
       sample_spatial_increment_pair,
       sample_ordinary_ust_spatial_increment_pair,
       sample_temporal_spatial_increment_pair,
       derive_sample_seeds,
       derive_temporal_sample_seeds,
       horizontal_edge_id,
       vertical_edge_id,
       parent_edge_id,
       face_id,
       primal_edge_count

const NORTH = UInt8(1)
const EAST  = UInt8(2)
const SOUTH = UInt8(3)
const WEST  = UInt8(4)
const DIRECTIONS = (NORTH, EAST, SOUTH, WEST)
const DX = (0, 1, 0, -1)
const DY = (1, 0, -1, 0)

abstract type AbstractSquareGridEnvironment end

"No-disorder environment with all four transition weights equal to one."
struct BaselineEnvironment <: AbstractSquareGridEnvironment end

"Directed site-i.i.d. positive weights; opposite orientations are independent."
struct DirectedSiteIIDEnvironment <: AbstractSquareGridEnvironment
    seed::UInt64
    distribution::Symbol
    parameter::Float64
end

"One positive conductance per unoriented geometric edge."
struct UndirectedConductanceEnvironment <: AbstractSquareGridEnvironment
    seed::UInt64
    distribution::Symbol
    parameter::Float64
end

"""
Independent temporal transition weights for one Wilson-tree sample.

Unlike the site and conductance environments, this is not a quenched spatial
environment: every raw transition draws and discards a new N/E/S/W vector.
The mutable RNG belongs to exactly one tree sample and must never be shared
between replicas when estimating the temporal marginal observable.
"""
mutable struct TemporalIIDEnvironment{R<:AbstractRNG} <: AbstractSquareGridEnvironment
    weight_seed::UInt64
    distribution::Symbol
    parameter::Union{Nothing,Float64}
    weight_rng::R
    sampled_weight_vectors::Int64
end

"Materialized transition CDFs for one fixed environment and square size."
struct MaterializedTransitionEnvironment <: AbstractSquareGridEnvironment
    L::Int
    cdf::Vector{NTuple{3,Float64}}
end

"Deterministic direction weights, primarily for exact tiny-grid tests."
struct FixedEnvironment <: AbstractSquareGridEnvironment
    weights::NTuple{4,Float64}
    function FixedEnvironment(weights::NTuple{4,T}) where {T<:Real}
        converted = ntuple(index -> Float64(weights[index]), 4)
        all(value -> isfinite(value) && value > 0, converted) ||
            throw(ArgumentError("fixed weights must be finite and positive"))
        new(converted)
    end
end

function check_distribution(distribution::Symbol, parameter::Real)
    distribution in (:gamma, :lognormal, :uniform) ||
        throw(ArgumentError("distribution must be :gamma, :lognormal, or :uniform"))
    parameter > 0 || throw(ArgumentError("distribution parameter must be positive"))
    return Float64(parameter)
end

function DirectedSiteIIDEnvironment(seed::Integer; distribution::Symbol=:gamma,
                                    parameter::Real=0.5)
    return DirectedSiteIIDEnvironment(
        UInt64(seed), distribution, check_distribution(distribution, parameter))
end

function UndirectedConductanceEnvironment(seed::Integer; distribution::Symbol=:gamma,
                                          parameter::Real=0.5)
    return UndirectedConductanceEnvironment(
        UInt64(seed), distribution, check_distribution(distribution, parameter))
end

function TemporalIIDEnvironment(seed::Integer; distribution::Symbol=:gamma,
                                parameter::Union{Nothing,Real}=0.5)
    distribution in (:baseline, :gamma, :lognormal, :uniform) || throw(ArgumentError(
        "temporal distribution must be :baseline, :gamma, :lognormal, or :uniform"))
    if distribution === :baseline
        parameter === nothing || Float64(parameter) > 0 ||
            throw(ArgumentError("baseline temporal weights need a positive placeholder parameter"))
        return TemporalIIDEnvironment(UInt64(seed), :baseline, nothing,
                                      Random.Xoshiro(UInt64(seed)), 0)
    end
    isnothing(parameter) && throw(ArgumentError("temporal disorder needs a parameter"))
    checked = Float64(parameter)
    checked > 0 || throw(ArgumentError("temporal distribution parameter must be positive"))
    # This is deliberately the historical temporal convention, not the
    # fixed-square-grid Uniform(0, upper) convention.
    distribution === :uniform && checked >= 1 && throw(ArgumentError(
        "temporal uniform parameter must satisfy 0 < a < 1 for Uniform(1-a, 1+a)"))
    return TemporalIIDEnvironment(UInt64(seed), distribution, checked,
                                  Random.Xoshiro(UInt64(seed)), 0)
end

@inline function splitmix64(value::UInt64)::UInt64
    value += 0x9e3779b97f4a7c15
    value = (value ⊻ (value >> 30)) * 0xbf58476d1ce4e5b9
    value = (value ⊻ (value >> 27)) * 0x94d049bb133111eb
    return value ⊻ (value >> 31)
end

@inline signed_word(value::Integer) = UInt64(reinterpret(UInt32, Int32(value)))

@inline function mix_words(seed::UInt64, words::UInt64...)
    value = seed
    for word in words
        value = splitmix64(value ⊻ splitmix64(word))
    end
    return value
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

@inline function sample_scalar_weight(seed::UInt64, distribution::Symbol,
                                      parameter::Float64)::Float64
    rng = Random.Xoshiro(seed)
    if distribution === :gamma
        # Mean-one Gamma(shape=k, scale=1/k).
        return rand_gamma(rng, parameter, inv(parameter))
    elseif distribution === :lognormal
        # Mean-one lognormal with standard deviation parameter sigma in log-space.
        return exp(parameter * randn(rng) - parameter^2 / 2)
    elseif distribution === :uniform
        # The global upper scale cancels from the tree law. Parameter 2 gives
        # a convenient mean-one Uniform(0, 2) representative.
        return max(parameter * rand(rng), floatmin(Float64))
    end
    throw(ArgumentError("unsupported distribution: $distribution"))
end

@inline function direction_delta(direction::UInt8)
    1 <= direction <= 4 || throw(ArgumentError("invalid direction $direction"))
    return DX[Int(direction)], DY[Int(direction)]
end

@inline reverse_direction(direction::UInt8) = UInt8(mod(Int(direction) + 1, 4) + 1)

@inline edge_weight(::BaselineEnvironment, x::Integer, y::Integer,
                    direction::UInt8, L::Integer) = 1.0

@inline function edge_weight(environment::FixedEnvironment, x::Integer, y::Integer,
                             direction::UInt8, L::Integer)
    return environment.weights[Int(direction)]
end

@inline function edge_weight(environment::DirectedSiteIIDEnvironment,
                             x::Integer, y::Integer, direction::UInt8, L::Integer)
    key = mix_words(
        environment.seed,
        0x6469726563746564,
        signed_word(x), signed_word(y), UInt64(direction),
    )
    return sample_scalar_weight(key, environment.distribution, environment.parameter)
end

@inline function edge_weight(environment::UndirectedConductanceEnvironment,
                             x::Integer, y::Integer, direction::UInt8, L::Integer)
    dx, dy = direction_delta(direction)
    neighbour_x = Int(x) + dx
    neighbour_y = Int(y) + dy
    ax, ay = Int(x), Int(y)
    bx, by = neighbour_x, neighbour_y
    if (bx, by) < (ax, ay)
        ax, ay, bx, by = bx, by, ax, ay
    end
    key = mix_words(
        environment.seed,
        0x756e646972656374,
        signed_word(ax), signed_word(ay), signed_word(bx), signed_word(by),
    )
    return sample_scalar_weight(key, environment.distribution, environment.parameter)
end

@inline function temporal_scalar_weight(
    rng::AbstractRNG,
    distribution::Symbol,
    parameter::Union{Nothing,Float64},
)::Float64
    distribution === :baseline && return 1.0
    distribution === :gamma && return rand_gamma(rng, parameter::Float64, inv(parameter::Float64))
    distribution === :lognormal && return exp(
        (parameter::Float64) * randn(rng) - (parameter::Float64)^2 / 2)
    distribution === :uniform && return 1.0 - (parameter::Float64) +
                                      2.0 * (parameter::Float64) * rand(rng)
    throw(ArgumentError("unsupported temporal distribution: $distribution"))
end

"Draw and discard one fresh directional vector for a single raw transition."
@inline function temporal_weights!(environment::TemporalIIDEnvironment)
    environment.sampled_weight_vectors += 1
    return ntuple(
        _ -> temporal_scalar_weight(
            environment.weight_rng, environment.distribution, environment.parameter),
        4,
    )
end

@inline grid_side(L::Integer) = 2 * Int(L) - 1
@inline vertex_count(L::Integer) = grid_side(L)^2

@inline function check_L(L::Integer)
    L >= 1 || throw(ArgumentError("L must be at least one"))
    L <= (typemax(Int32) - 1) ÷ 2 || throw(ArgumentError("L is too large"))
    return Int(L)
end

@inline function vertex_index(L::Integer, x::Integer, y::Integer)::Int
    size = grid_side(L)
    (-L < x < L && -L < y < L) ||
        throw(ArgumentError("($x,$y) is not an interior vertex for L=$L"))
    return (Int(y) + Int(L) - 1) * size + (Int(x) + Int(L))
end

@inline function vertex_coordinates(L::Integer, index::Integer)
    size = grid_side(L)
    1 <= index <= size^2 || throw(BoundsError(1:size^2, index))
    zero_based = Int(index) - 1
    column = mod(zero_based, size) + 1
    row = fld(zero_based, size) + 1
    return (x=column - Int(L), y=row - Int(L))
end

@inline function neighbour_index(L::Integer, x::Integer, y::Integer,
                                 direction::UInt8)::Int
    dx, dy = direction_delta(direction)
    next_x = Int(x) + dx
    next_y = Int(y) + dy
    return abs(next_x) == L || abs(next_y) == L ? 0 :
           vertex_index(L, next_x, next_y)
end

@inline function choose_direction(environment::AbstractSquareGridEnvironment,
                                  x::Integer, y::Integer, L::Integer,
                                  rng::AbstractRNG)::UInt8
    weights = ntuple(
        index -> edge_weight(environment, x, y, UInt8(index), L), 4)
    total = weights[1] + weights[2] + weights[3] + weights[4]
    isfinite(total) && total > 0 ||
        throw(ErrorException("transition weights must have positive finite sum"))
    threshold = rand(rng) * total
    threshold <= weights[1] && return NORTH
    threshold <= weights[1] + weights[2] && return EAST
    threshold <= weights[1] + weights[2] + weights[3] && return SOUTH
    return WEST
end

@inline function choose_direction(environment::MaterializedTransitionEnvironment,
                                  x::Integer, y::Integer, L::Integer,
                                  rng::AbstractRNG)::UInt8
    Int(L) == environment.L || throw(ArgumentError("materialized environment has wrong L"))
    north, east, south = environment.cdf[vertex_index(L, x, y)]
    draw = rand(rng)
    draw <= north && return NORTH
    draw <= east && return EAST
    draw <= south && return SOUTH
    return WEST
end

@inline function choose_direction(environment::TemporalIIDEnvironment,
                                  x::Integer, y::Integer, L::Integer,
                                  rng::AbstractRNG)::UInt8
    # x, y, and L intentionally do not enter the weight draw: revisits receive
    # fresh temporal disorder rather than a cached site value.
    weights = temporal_weights!(environment)
    total = weights[1] + weights[2] + weights[3] + weights[4]
    isfinite(total) && total > 0 || throw(ErrorException(
        "temporal transition weights must have positive finite sum"))
    threshold = rand(rng) * total
    threshold <= weights[1] && return NORTH
    threshold <= weights[1] + weights[2] && return EAST
    threshold <= weights[1] + weights[2] + weights[3] && return SOUTH
    return WEST
end

materialize_environment(environment::Union{BaselineEnvironment,FixedEnvironment},
                        L::Integer) = environment

"Precompute each row-normalised transition law once for reuse by both replicas."
function materialize_environment(
    environment::Union{DirectedSiteIIDEnvironment,UndirectedConductanceEnvironment},
    L::Integer,
)
    checked_L = check_L(L)
    cdf = Vector{NTuple{3,Float64}}(undef, vertex_count(checked_L))
    for index in eachindex(cdf)
        point = vertex_coordinates(checked_L, index)
        weights = ntuple(
            direction -> edge_weight(
                environment, point.x, point.y, UInt8(direction), checked_L),
            4,
        )
        total = sum(weights)
        isfinite(total) && total > 0 || throw(ErrorException(
            "transition weights must have positive finite sum"))
        cdf[index] = (
            weights[1] / total,
            (weights[1] + weights[2]) / total,
            (weights[1] + weights[2] + weights[3]) / total,
        )
    end
    return MaterializedTransitionEnvironment(checked_L, cdf)
end

"A rooted tree sampled on the wired square. Parent directions point toward the root."
struct TreeSample
    L::Int
    parent_direction::Vector{UInt8}
    raw_steps::Int64
    branches::Int
end

@inline function parent_direction(tree::TreeSample, x::Integer, y::Integer)::UInt8
    return tree.parent_direction[vertex_index(tree.L, x, y)]
end

"""
    sample_full_tree(environment, L, rng; max_steps_per_walk=nothing)

Sample a complete weighted rooted spanning tree with Wilson's algorithm.
The order of starting vertices is deterministic, making fixed-seed results
independent of Julia thread scheduling at the campaign level.
"""
function sample_full_tree(environment::AbstractSquareGridEnvironment,
                          L::Integer, rng::AbstractRNG;
                          max_steps_per_walk::Union{Nothing,Integer}=nothing)
    checked_L = check_L(L)
    count = vertex_count(checked_L)
    in_tree = falses(count)
    parents = zeros(UInt8, count)
    marks = zeros(Int32, count)
    positions = zeros(Int32, count)
    walk_number = Int32(0)
    raw_steps_total = Int64(0)
    branch_count = 0
    cap = isnothing(max_steps_per_walk) ?
        max(100_000, 2_000 * checked_L^2) : Int(max_steps_per_walk)
    cap > 0 || throw(ArgumentError("max_steps_per_walk must be positive"))

    for start in 1:count
        in_tree[start] && continue
        walk_number == typemax(Int32) && error("too many Wilson branches")
        walk_number += Int32(1)
        branch_count += 1
        path = Int[start]
        directions = UInt8[]
        marks[start] = walk_number
        positions[start] = Int32(1)
        current = start
        branch_steps = 0

        while true
            branch_steps += 1
            branch_steps <= cap || throw(ErrorException(
                "Wilson walk from vertex $start exceeded $cap raw steps"))
            raw_steps_total += 1
            point = vertex_coordinates(checked_L, current)
            direction = choose_direction(
                environment, point.x, point.y, checked_L, rng)
            next_index = neighbour_index(
                checked_L, point.x, point.y, direction)

            if next_index == 0 || in_tree[next_index]
                push!(directions, direction)
                break
            elseif marks[next_index] == walk_number
                keep = Int(positions[next_index])
                @inbounds for removed in path[(keep + 1):end]
                    marks[removed] = Int32(0)
                    positions[removed] = Int32(0)
                end
                resize!(path, keep)
                resize!(directions, keep - 1)
                current = next_index
            else
                push!(directions, direction)
                push!(path, next_index)
                marks[next_index] = walk_number
                positions[next_index] = Int32(length(path))
                current = next_index
            end
        end

        length(directions) == length(path) ||
            error("internal Wilson path/direction mismatch")
        @inbounds for index in eachindex(path)
            vertex = path[index]
            parents[vertex] = directions[index]
            in_tree[vertex] = true
        end
    end

    all(!=(UInt8(0)), parents) || error("Wilson sampler left an unset parent")
    return TreeSample(checked_L, parents, raw_steps_total, branch_count)
end

"Check that every vertex has one parent and reaches the wired root without a cycle."
function validate_tree(tree::TreeSample)
    count = length(tree.parent_direction)
    count == vertex_count(tree.L) ||
        return (valid=false, reason="wrong parent vector length", vertices=count)
    state = zeros(UInt8, count) # 0 unseen, 1 active, 2 complete
    for start in 1:count
        state[start] == 2 && continue
        current = start
        stack = Int[]
        while current != 0 && state[current] == 0
            direction = tree.parent_direction[current]
            direction in DIRECTIONS ||
                return (valid=false, reason="invalid parent direction", vertices=count)
            state[current] = 1
            push!(stack, current)
            point = vertex_coordinates(tree.L, current)
            current = neighbour_index(tree.L, point.x, point.y, direction)
        end
        if current != 0 && state[current] == 1
            return (valid=false, reason="directed cycle", vertices=count)
        end
        for vertex in stack
            state[vertex] = 2
        end
    end
    return (valid=true, reason="ok", vertices=count)
end

@inline function horizontal_edge_id(L::Integer, x::Integer, y::Integer)::Int
    n = grid_side(L)
    (-L < x < L - 1 && -L < y < L) ||
        throw(ArgumentError("invalid horizontal edge ($x,$y) for L=$L"))
    return (Int(y) + Int(L) - 1) * (n - 1) + (Int(x) + Int(L))
end

@inline function vertical_edge_id(L::Integer, x::Integer, y::Integer)::Int
    n = grid_side(L)
    (-L < x < L && -L < y < L - 1) ||
        throw(ArgumentError("invalid vertical edge ($x,$y) for L=$L"))
    horizontal_count = n * (n - 1)
    return horizontal_count +
           (Int(y) + Int(L) - 1) * n + (Int(x) + Int(L))
end

@inline function north_boundary_edge_id(L::Integer, x::Integer)::Int
    n = grid_side(L)
    -L < x < L || throw(ArgumentError("invalid north boundary edge"))
    return 2 * n * (n - 1) + (Int(x) + Int(L))
end

@inline function south_boundary_edge_id(L::Integer, x::Integer)::Int
    n = grid_side(L)
    -L < x < L || throw(ArgumentError("invalid south boundary edge"))
    return 2 * n * (n - 1) + n + (Int(x) + Int(L))
end

@inline function east_boundary_edge_id(L::Integer, y::Integer)::Int
    n = grid_side(L)
    -L < y < L || throw(ArgumentError("invalid east boundary edge"))
    return 2 * n * (n - 1) + 2n + (Int(y) + Int(L))
end

@inline function west_boundary_edge_id(L::Integer, y::Integer)::Int
    n = grid_side(L)
    -L < y < L || throw(ArgumentError("invalid west boundary edge"))
    return 2 * n * (n - 1) + 3n + (Int(y) + Int(L))
end

@inline primal_edge_count(L::Integer) = 4 * Int(L) * (2 * Int(L) - 1)

@inline function parent_edge_id(L::Integer, x::Integer, y::Integer,
                                direction::UInt8)::Int
    if direction == NORTH
        return y == L - 1 ? north_boundary_edge_id(L, x) :
               vertical_edge_id(L, x, y)
    elseif direction == EAST
        return x == L - 1 ? east_boundary_edge_id(L, y) :
               horizontal_edge_id(L, x, y)
    elseif direction == SOUTH
        return y == -L + 1 ? south_boundary_edge_id(L, x) :
               vertical_edge_id(L, x, y - 1)
    elseif direction == WEST
        return x == -L + 1 ? west_boundary_edge_id(L, y) :
               horizontal_edge_id(L, x - 1, y)
    end
    throw(ArgumentError("invalid direction $direction"))
end

@inline function face_id(L::Integer, lower_left_x::Integer,
                         lower_left_y::Integer)::Int
    m = 2 * Int(L)
    (-L <= lower_left_x < L && -L <= lower_left_y < L) ||
        throw(ArgumentError("invalid face ($lower_left_x,$lower_left_y) for L=$L"))
    return (Int(lower_left_y) + Int(L)) * m +
           (Int(lower_left_x) + Int(L)) + 1
end

"Exact Temperley matching represented by primal-tree and dual-tree ownership."
struct TemperleyMatching
    tree::TreeSample
    dual_parent_edge::Vector{Int}
    edge_owner::Vector{UInt8} # 1 vertex node, 2 face node
    outer_face::Int
end

@inline function tree_uses_horizontal(tree::TreeSample, x::Int, y::Int)
    return parent_direction(tree, x, y) == EAST ||
           parent_direction(tree, x + 1, y) == WEST
end

@inline function tree_uses_vertical(tree::TreeSample, x::Int, y::Int)
    return parent_direction(tree, x, y) == NORTH ||
           parent_direction(tree, x, y + 1) == SOUTH
end

"Build the complementary dual tree and the exact generalized Temperley matching."
function build_temperley_matching(tree::TreeSample)
    validation = validate_tree(tree)
    validation.valid || throw(ArgumentError("invalid primal tree: $(validation.reason)"))
    L = tree.L
    n = grid_side(L)
    faces_per_side = 2L
    number_of_faces = faces_per_side^2
    outer_face = face_id(L, -L, -L)
    adjacency = [Tuple{Int,Int}[] for _ in 1:number_of_faces]
    complement_count = 0

    function add_complement(edge_id::Int, first_face::Int, second_face::Int,
                            selected::Bool)
        if !selected
            push!(adjacency[first_face], (second_face, edge_id))
            push!(adjacency[second_face], (first_face, edge_id))
            complement_count += 1
        end
        return nothing
    end

    # Interior horizontal edges, then interior vertical edges. This order must
    # agree with the public edge-ID formulas above.
    for y in (-L + 1):(L - 1), x in (-L + 1):(L - 2)
        edge = horizontal_edge_id(L, x, y)
        add_complement(
            edge,
            face_id(L, x, y),
            face_id(L, x, y - 1),
            tree_uses_horizontal(tree, x, y),
        )
    end
    for y in (-L + 1):(L - 2), x in (-L + 1):(L - 1)
        edge = vertical_edge_id(L, x, y)
        add_complement(
            edge,
            face_id(L, x - 1, y),
            face_id(L, x, y),
            tree_uses_vertical(tree, x, y),
        )
    end

    for x in (-L + 1):(L - 1)
        add_complement(
            north_boundary_edge_id(L, x),
            face_id(L, x - 1, L - 1),
            face_id(L, x, L - 1),
            parent_direction(tree, x, L - 1) == NORTH,
        )
    end
    for x in (-L + 1):(L - 1)
        add_complement(
            south_boundary_edge_id(L, x),
            face_id(L, x - 1, -L),
            face_id(L, x, -L),
            parent_direction(tree, x, -L + 1) == SOUTH,
        )
    end
    for y in (-L + 1):(L - 1)
        add_complement(
            east_boundary_edge_id(L, y),
            face_id(L, L - 1, y),
            face_id(L, L - 1, y - 1),
            parent_direction(tree, L - 1, y) == EAST,
        )
    end
    for y in (-L + 1):(L - 1)
        add_complement(
            west_boundary_edge_id(L, y),
            face_id(L, -L, y),
            face_id(L, -L, y - 1),
            parent_direction(tree, -L + 1, y) == WEST,
        )
    end

    complement_count == number_of_faces - 1 || throw(ErrorException(
        "tree complement has $complement_count edges; expected $(number_of_faces - 1)"))

    dual_parent_edge = zeros(Int, number_of_faces)
    visited = falses(number_of_faces)
    queue = Vector{Int}(undef, number_of_faces)
    head = 1
    tail = 1
    queue[1] = outer_face
    visited[outer_face] = true
    while head <= tail
        face = queue[head]
        head += 1
        for (neighbour, edge) in adjacency[face]
            visited[neighbour] && continue
            tail += 1
            queue[tail] = neighbour
            visited[neighbour] = true
            dual_parent_edge[neighbour] = edge
        end
    end
    all(visited) || throw(ErrorException("complementary dual graph is disconnected"))

    owner = zeros(UInt8, primal_edge_count(L))
    for index in eachindex(tree.parent_direction)
        point = vertex_coordinates(L, index)
        edge = parent_edge_id(L, point.x, point.y, tree.parent_direction[index])
        owner[edge] == 0 || throw(ErrorException("primal edge matched twice"))
        owner[edge] = 1
    end
    for face in 1:number_of_faces
        face == outer_face && continue
        edge = dual_parent_edge[face]
        edge != 0 || throw(ErrorException("bounded face has no dual parent"))
        owner[edge] == 0 || throw(ErrorException("edge matched by vertex and face"))
        owner[edge] = 2
    end
    all(!=(UInt8(0)), owner) || throw(ErrorException("Temperley matching is not perfect"))
    return TemperleyMatching(tree, dual_parent_edge, owner, outer_face)
end

function validate_matching(matching::TemperleyMatching)
    L = matching.tree.L
    expected_edges = primal_edge_count(L)
    length(matching.edge_owner) == expected_edges ||
        return (valid=false, reason="wrong edge-owner length", edges=length(matching.edge_owner))
    all(owner -> owner == 1 || owner == 2, matching.edge_owner) ||
        return (valid=false, reason="unmatched or multiply matched edge", edges=expected_edges)
    vertex_owned = count(==(UInt8(1)), matching.edge_owner)
    face_owned = count(==(UInt8(2)), matching.edge_owner)
    vertex_owned == vertex_count(L) ||
        return (valid=false, reason="wrong vertex matching count", edges=expected_edges)
    face_owned == (2L)^2 - 1 ||
        return (valid=false, reason="wrong face matching count", edges=expected_edges)

    # Reconstruct every matched incidence from the primal and dual parents.
    # This is stronger than checking owner totals: it detects a duplicate or an
    # owner label that is not incident to the white node claiming it.
    claimed = falses(expected_edges)
    for index in eachindex(matching.tree.parent_direction)
        point = vertex_coordinates(L, index)
        edge = parent_edge_id(
            L, point.x, point.y, matching.tree.parent_direction[index])
        matching.edge_owner[edge] == 1 ||
            return (valid=false, reason="vertex parent has wrong owner", edges=expected_edges)
        claimed[edge] &&
            return (valid=false, reason="matching edge claimed twice", edges=expected_edges)
        claimed[edge] = true
    end
    number_of_faces = (2L)^2
    length(matching.dual_parent_edge) == number_of_faces ||
        return (valid=false, reason="wrong dual-parent length", edges=expected_edges)
    matching.dual_parent_edge[matching.outer_face] == 0 ||
        return (valid=false, reason="outer face has a dual parent", edges=expected_edges)
    for face in 1:number_of_faces
        face == matching.outer_face && continue
        edge = matching.dual_parent_edge[face]
        1 <= edge <= expected_edges ||
            return (valid=false, reason="invalid face parent edge", edges=expected_edges)
        matching.edge_owner[edge] == 2 ||
            return (valid=false, reason="face parent has wrong owner", edges=expected_edges)
        claimed[edge] &&
            return (valid=false, reason="matching edge claimed twice", edges=expected_edges)
        claimed[edge] = true
    end
    all(claimed) ||
        return (valid=false, reason="unclaimed matching edge", edges=expected_edges)
    return (
        valid=true,
        reason="ok",
        edges=expected_edges,
        vertex_matches=vertex_owned,
        face_matches=face_owned,
    )
end

@inline function interior_black_edge_id(L::Int, scaled_x::Int, scaled_y::Int)
    if isodd(scaled_x) && iseven(scaled_y)
        return horizontal_edge_id(L, fld(scaled_x, 2), fld(scaled_y, 2))
    elseif iseven(scaled_x) && isodd(scaled_y)
        return vertical_edge_id(L, fld(scaled_x, 2), fld(scaled_y, 2))
    end
    throw(ArgumentError("($scaled_x,$scaled_y) is not an interior black edge-node"))
end

"Return `(matched, white_to_black_dx, white_to_black_dy)` for one graph link."
function matching_link_data(
    matching::TemperleyMatching,
    first::Tuple{Int,Int},
    second::Tuple{Int,Int},
)
    abs(first[1] - second[1]) + abs(first[2] - second[2]) == 1 ||
        throw(ArgumentError("Temperley graph links must join nearest neighbours"))
    first_is_white = iseven(first[1]) == iseven(first[2])
    second_is_white = iseven(second[1]) == iseven(second[2])
    first_is_white != second_is_white ||
        throw(ArgumentError("link does not join a white node to a black node"))
    white = first_is_white ? first : second
    black = first_is_white ? second : first
    edge = interior_black_edge_id(matching.tree.L, black...)
    matched = if iseven(white[1])
        x, y = fld(white[1], 2), fld(white[2], 2)
        parent_edge_id(matching.tree.L, x, y, parent_direction(matching.tree, x, y)) == edge
    else
        face = face_id(
            matching.tree.L, fld(white[1] - 1, 2), fld(white[2] - 1, 2))
        face != matching.outer_face && matching.dual_parent_edge[face] == edge
    end
    return (
        matched=matched,
        dx=black[1] - white[1],
        dy=black[2] - white[2],
    )
end

"""
    dimer_height_increment_along_path(matching, path)

Independently integrate matching flow minus the degree-four reference flow
along a nearest-neighbour path of dual cells in the scaled Temperley lattice.
Cell `(i,j)` is the unit square with lower-left graph coordinate `(i,j)`.
The result is exact in quarter-height units and is path-independent away from
the removed root/outer-face singularity.
"""
function dimer_height_increment_along_path(
    matching::TemperleyMatching,
    path::AbstractVector{<:Tuple{Int,Int}},
)
    length(path) >= 2 || throw(ArgumentError("height path needs at least two cells"))
    quarter_units = 0
    for index in 1:(length(path) - 1)
        first, second = path[index], path[index + 1]
        step_x, step_y = second[1] - first[1], second[2] - first[2]
        abs(step_x) + abs(step_y) == 1 ||
            throw(ArgumentError("height path must use nearest-neighbour dual steps"))
        if step_x != 0
            crossing_x = max(first[1], second[1])
            lower_y = first[2]
            link_first = (crossing_x, lower_y)
            link_second = (crossing_x, lower_y + 1)
        else
            lower_x = first[1]
            crossing_y = max(first[2], second[2])
            link_first = (lower_x, crossing_y)
            link_second = (lower_x + 1, crossing_y)
        end
        link = matching_link_data(matching, link_first, link_second)
        flow_quarters = link.matched ? 3 : -1
        intersection_sign = step_x * link.dy - step_y * link.dx
        quarter_units += intersection_sign * flow_quarters
    end
    return quarter_units // 4
end

"Probe vertices separated by `separation` and centred as closely as parity allows."
function symmetric_probe_vertices(L::Integer, separation::Integer)
    separation > 0 || throw(ArgumentError("separation must be positive"))
    left_x = -fld(Int(separation), 2)
    right_x = left_x + Int(separation)
    (-L < left_x < L && -L < right_x < L) ||
        throw(ArgumentError("separation $separation leaves the interior domain"))
    return (left=(x=left_x, y=0), right=(x=right_x, y=0))
end

"""
    spatial_height_increment(matching, separation)

Return the exact dimer-height increment along the horizontal dual segment at
scaled height `1/2`, from the left central probe to the right central probe.
The segment alternately crosses vertex-edge and face-edge links of the
Temperley graph. With white-to-black flow orientation, upward matched links
contribute `+1` and downward matched links contribute `-1`.

Every crossed white node has degree four, with equally many upward and
downward crossings. The deterministic `1/4` reference-flow contribution
therefore cancels exactly, leaving this integer height increment.
"""
function spatial_height_increment(matching::TemperleyMatching,
                                  separation::Integer)::Int
    L = matching.tree.L
    probes = symmetric_probe_vertices(L, separation)
    total = 0
    for scaled_x in (2 * probes.left.x + 1):(2 * probes.right.x)
        if iseven(scaled_x)
            x = div(scaled_x, 2)
            parent_direction(matching.tree, x, 0) == NORTH && (total += 1)
        else
            x = div(scaled_x - 1, 2)
            face = face_id(L, x, 0)
            edge = horizontal_edge_id(L, x, 0)
            matching.dual_parent_edge[face] == edge && (total -= 1)
        end
    end
    return total
end

function derive_sample_seeds(sample_seed::Integer)
    seed = UInt64(sample_seed)
    return (
        environment=splitmix64(seed ⊻ 0x656e7669726f6e6d),
        replica_1=splitmix64(seed ⊻ 0x7265706c69636131),
        replica_2=splitmix64(seed ⊻ 0x7265706c69636132),
    )
end

"Separate temporal weight and direction streams for two independent tree replicas."
function derive_temporal_sample_seeds(sample_seed::Integer)
    seed = UInt64(sample_seed)
    return (
        replica_1_weight=splitmix64(seed ⊻ 0x74656d705f726577),
        replica_1_direction=splitmix64(seed ⊻ 0x74656d705f726564),
        replica_2_weight=splitmix64(seed ⊻ 0x74656d705f726631),
        replica_2_direction=splitmix64(seed ⊻ 0x74656d705f726632),
    )
end

function make_environment(model::Symbol, environment_seed::UInt64,
                          distribution::Symbol, parameter::Real)
    model === :baseline && return BaselineEnvironment()
    model === :directed_site_iid && return DirectedSiteIIDEnvironment(
        environment_seed; distribution=distribution, parameter=parameter)
    model === :undirected_conductance && return UndirectedConductanceEnvironment(
        environment_seed; distribution=distribution, parameter=parameter)
    throw(ArgumentError(
        "environment model must be :baseline, :directed_site_iid, or :undirected_conductance"))
end

"""
    sample_spatial_increment_pair(sample_seed, L, separations; ...)

Draw one environment and two conditionally independent complete trees, convert
both to exact Temperley matchings, and evaluate all separations on those same
two replicas. Returns `(rows, diagnostics)`.
"""
function sample_spatial_increment_pair(
    sample_seed::Integer,
    L::Integer,
    separations::AbstractVector{<:Integer};
    environment_model::Symbol=:directed_site_iid,
    distribution::Symbol=:gamma,
    parameter::Real=0.5,
    max_steps_per_walk::Union{Nothing,Integer}=nothing,
)
    isempty(separations) && throw(ArgumentError("at least one separation is required"))
    length(unique(separations)) == length(separations) ||
        throw(ArgumentError("separations must be unique"))
    seeds = derive_sample_seeds(sample_seed)
    environment = make_environment(
        environment_model, seeds.environment, distribution, parameter)
    environment = materialize_environment(environment, L)
    first_tree = sample_full_tree(
        environment, L, Random.Xoshiro(seeds.replica_1); max_steps_per_walk=max_steps_per_walk)
    second_tree = sample_full_tree(
        environment, L, Random.Xoshiro(seeds.replica_2); max_steps_per_walk=max_steps_per_walk)
    first_matching = build_temperley_matching(first_tree)
    second_matching = build_temperley_matching(second_tree)
    rows = [
        begin
            probes = symmetric_probe_vertices(L, separation)
            first = spatial_height_increment(first_matching, separation)
            second = spatial_height_increment(second_matching, separation)
            (
                separation=Int(separation),
                left_x=probes.left.x,
                right_x=probes.right.x,
                increment_1=first,
                increment_2=second,
                difference=first - second,
            )
        end
        for separation in separations
    ]
    diagnostics = (
        environment_seed=seeds.environment,
        replica_1_seed=seeds.replica_1,
        replica_2_seed=seeds.replica_2,
        raw_steps_1=first_tree.raw_steps,
        raw_steps_2=second_tree.raw_steps,
        branches_1=first_tree.branches,
        branches_2=second_tree.branches,
    )
    return (rows=rows, diagnostics=diagnostics)
end

"""
    sample_ordinary_ust_spatial_increment_pair(sample_seed, L, separations; ...)

Draw two independent ordinary simple-random-walk UST samples, then send both
through the same `TreeSample -> TemperleyMatching -> spatial_height_increment`
path as the temporal refreshed-weight experiment.  The recorded environment
seed is an unused deterministic namespace marker retained for the common
campaign diagnostic schema; no spatial environment is constructed or shared.
"""
function sample_ordinary_ust_spatial_increment_pair(
    sample_seed::Integer,
    L::Integer,
    separations::AbstractVector{<:Integer};
    max_steps_per_walk::Union{Nothing,Integer}=nothing,
)
    isempty(separations) && throw(ArgumentError("at least one separation is required"))
    length(unique(separations)) == length(separations) ||
        throw(ArgumentError("separations must be unique"))
    seeds = derive_sample_seeds(sample_seed)
    first_tree = sample_full_tree(
        BaselineEnvironment(), L, Random.Xoshiro(seeds.replica_1);
        max_steps_per_walk=max_steps_per_walk)
    second_tree = sample_full_tree(
        BaselineEnvironment(), L, Random.Xoshiro(seeds.replica_2);
        max_steps_per_walk=max_steps_per_walk)
    first_matching = build_temperley_matching(first_tree)
    second_matching = build_temperley_matching(second_tree)
    rows = [
        begin
            probes = symmetric_probe_vertices(L, separation)
            first = spatial_height_increment(first_matching, separation)
            second = spatial_height_increment(second_matching, separation)
            (
                separation=Int(separation),
                left_x=probes.left.x,
                right_x=probes.right.x,
                increment_1=first,
                increment_2=second,
                difference=first - second,
            )
        end
        for separation in separations
    ]
    diagnostics = (
        environment_seed=seeds.environment,
        replica_1_seed=seeds.replica_1,
        replica_2_seed=seeds.replica_2,
        raw_steps_1=first_tree.raw_steps,
        raw_steps_2=second_tree.raw_steps,
        branches_1=first_tree.branches,
        branches_2=second_tree.branches,
    )
    return (rows=rows, diagnostics=diagnostics)
end

"""
    sample_temporal_spatial_increment_pair(sample_seed, L, separations; ...)

Draw two fully independent Wilson trees whose directional weights are refreshed
at every raw transition.  The pair is an independence control for the temporal
marginal variance: it has no shared quenched environment and must not be
interpreted as a disorder-covariance estimator.
"""
function sample_temporal_spatial_increment_pair(
    sample_seed::Integer,
    L::Integer,
    separations::AbstractVector{<:Integer};
    distribution::Symbol=:gamma,
    parameter::Union{Nothing,Real}=0.5,
    max_steps_per_walk::Union{Nothing,Integer}=nothing,
)
    isempty(separations) && throw(ArgumentError("at least one separation is required"))
    length(unique(separations)) == length(separations) ||
        throw(ArgumentError("separations must be unique"))
    seeds = derive_temporal_sample_seeds(sample_seed)
    first_environment = TemporalIIDEnvironment(
        seeds.replica_1_weight; distribution=distribution, parameter=parameter)
    second_environment = TemporalIIDEnvironment(
        seeds.replica_2_weight; distribution=distribution, parameter=parameter)
    first_tree = sample_full_tree(
        first_environment, L, Random.Xoshiro(seeds.replica_1_direction);
        max_steps_per_walk=max_steps_per_walk)
    second_tree = sample_full_tree(
        second_environment, L, Random.Xoshiro(seeds.replica_2_direction);
        max_steps_per_walk=max_steps_per_walk)
    first_matching = build_temperley_matching(first_tree)
    second_matching = build_temperley_matching(second_tree)
    rows = [
        begin
            probes = symmetric_probe_vertices(L, separation)
            first = spatial_height_increment(first_matching, separation)
            second = spatial_height_increment(second_matching, separation)
            (
                separation=Int(separation),
                left_x=probes.left.x,
                right_x=probes.right.x,
                increment_1=first,
                increment_2=second,
                difference=first - second,
            )
        end
        for separation in separations
    ]
    diagnostics = (
        replica_1_weight_seed=seeds.replica_1_weight,
        replica_1_direction_seed=seeds.replica_1_direction,
        replica_2_weight_seed=seeds.replica_2_weight,
        replica_2_direction_seed=seeds.replica_2_direction,
        sampled_weight_vectors_1=first_environment.sampled_weight_vectors,
        sampled_weight_vectors_2=second_environment.sampled_weight_vectors,
        raw_steps_1=first_tree.raw_steps,
        raw_steps_2=second_tree.raw_steps,
        branches_1=first_tree.branches,
        branches_2=second_tree.branches,
    )
    return (rows=rows, diagnostics=diagnostics)
end

end # module SquareGrid
