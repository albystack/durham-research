const PointKey = UInt64
const MAX_SITE_CACHE_CAPACITY = 1 << 23

@inline function pack_point(x::Int32, y::Int32)::PointKey
    return (UInt64(reinterpret(UInt32, x)) << 32) | UInt64(reinterpret(UInt32, y))
end

@inline function unpack_point(point::PointKey)::Tuple{Int32,Int32}
    x = reinterpret(Int32, UInt32(point >> 32))
    y = reinterpret(Int32, UInt32(point & 0xffffffff))
    return x, y
end

@inline function is_boundary(x::Integer, y::Integer, L::Integer)::Bool
    return abs(x) == L || abs(y) == L
end

"A bounded direct-mapped cache; evicted sites are regenerated deterministically."
mutable struct SiteWeightCache
    keys::Vector{PointKey}
    values::Vector{NTuple{4,Float64}}
    occupied::BitVector
    count::Int
end

function SiteWeightCache(capacity::Integer)
    capacity >= 1 && ispow2(capacity) ||
        throw(ArgumentError("site cache capacity must be a positive power of two"))
    return SiteWeightCache(Vector{PointKey}(undef, capacity),
                           Vector{NTuple{4,Float64}}(undef, capacity),
                           falses(capacity), 0)
end

Base.length(cache::SiteWeightCache) = cache.count

@inline function splitmix64(value::UInt64)::UInt64
    value += 0x9e3779b97f4a7c15
    value = (value ⊻ (value >> 30)) * 0xbf58476d1ce4e5b9
    value = (value ⊻ (value >> 27)) * 0x94d049bb133111eb
    return value ⊻ (value >> 31)
end

function site_cache_capacity(L::Integer)::Int
    target = clamp(cld(Int128(L) * Int128(L), 4), Int128(1 << 10),
                   Int128(MAX_SITE_CACHE_CAPACITY))
    return Int(nextpow(2, target))
end

mutable struct SiteIIDEnvironment
    environment_seed::UInt64
    model::Symbol
    parameter::Union{Nothing,Float64}
    weights::SiteWeightCache
    weight_seed::UInt64
end

function SiteIIDEnvironment(seed::Integer, model::Symbol,
                            parameter::Union{Nothing,Float64}; cache_capacity::Integer=1 << 16)
    environment_seed = UInt64(seed)
    weight_seed = stable_seed("site_iid_weights", environment_seed, model, parameter)
    return SiteIIDEnvironment(environment_seed, model, parameter,
                              SiteWeightCache(cache_capacity), weight_seed)
end

function SiteIIDEnvironment(seed::Integer, distribution::AbstractString, params::AbstractDict;
                            cache_capacity::Integer=1 << 16)
    model, parameter = distribution_spec(distribution, params)
    return SiteIIDEnvironment(seed, model, parameter; cache_capacity)
end

@inline function site_weights!(environment::SiteIIDEnvironment, point::PointKey)
    environment.model === :baseline && return (1.0, 1.0, 1.0, 1.0)
    cache = environment.weights
    index = Int(splitmix64(point) & UInt64(length(cache.keys) - 1)) + 1
    if @inbounds cache.occupied[index] && cache.keys[index] == point
        return @inbounds cache.values[index]
    end
    rng = Random.Xoshiro(splitmix64(environment.weight_seed ⊻ point))
    weights = ntuple(_ -> sample_weight(rng, environment.model, environment.parameter), 4)
    @inbounds begin
        if !cache.occupied[index]
            cache.occupied[index] = true
            cache.count += 1
        end
        cache.keys[index] = point
        cache.values[index] = weights
    end
    return weights
end

@inline function next_point(current::PointKey, environment::SiteIIDEnvironment,
                            walk_rng::AbstractRNG)::PointKey
    weights = site_weights!(environment, current)
    threshold = rand(walk_rng) * (weights[1] + weights[2] + weights[3] + weights[4])
    direction = if threshold <= weights[1]
        1 # north
    elseif threshold <= weights[1] + weights[2]
        2 # east
    elseif threshold <= weights[1] + weights[2] + weights[3]
        3 # south
    else
        4 # west
    end
    x, y = unpack_point(current)
    direction == 1 && return pack_point(x, y + Int32(1))
    direction == 2 && return pack_point(x + Int32(1), y)
    direction == 3 && return pack_point(x, y - Int32(1))
    return pack_point(x - Int32(1), y)
end

"Run a random walk to the square boundary while maintaining loop erasure online."
function loop_erased_walk(L::Integer, environment::SiteIIDEnvironment,
                          walk_rng::AbstractRNG; max_steps::Union{Nothing,Integer}=nothing)
    L >= 1 || throw(ArgumentError("L must be at least 1"))
    L <= typemax(Int32) || throw(ArgumentError("L exceeds Int32 coordinate range"))
    cap = max_steps === nothing ? max(100_000, 2_000 * Int(L)^2) : Int(max_steps)

    origin = pack_point(Int32(0), Int32(0))
    path = PointKey[origin]
    positions = Dict{PointKey,Int32}(origin => Int32(1))
    current = origin

    for raw_steps in 0:cap
        x, y = unpack_point(current)
        is_boundary(x, y, L) && return path, raw_steps
        raw_steps == cap && break

        current = next_point(current, environment, walk_rng)
        previous_index = get(positions, current, Int32(0))
        if previous_index != 0
            keep = Int(previous_index)
            @inbounds for index in (keep + 1):length(path)
                delete!(positions, path[index])
            end
            resize!(path, keep)
        else
            push!(path, current)
            positions[current] = Int32(length(path))
        end
    end
    throw(ErrorException("walk did not hit the boundary within $cap steps"))
end

@inline function direction_code(a::PointKey, b::PointKey)::Int8
    ax, ay = unpack_point(a)
    bx, by = unpack_point(b)
    bx == ax + 1 && by == ay && return Int8(0) # east
    bx == ax && by == ay + 1 && return Int8(1) # north
    bx == ax - 1 && by == ay && return Int8(2) # west
    bx == ax && by == ay - 1 && return Int8(3) # south
    throw(ArgumentError("path contains a non-nearest-neighbour step"))
end

"Return left quarter-turns minus right quarter-turns along a path."
function winding(path::AbstractVector{PointKey})::Int
    length(path) < 3 && return 0
    turns = 0
    before = direction_code(path[1], path[2])
    @inbounds for index in 2:(length(path) - 1)
        after = direction_code(path[index], path[index + 1])
        turn = mod(Int(after) - Int(before), 4)
        turn == 1 && (turns += 1)
        turn == 3 && (turns -= 1)
        turn == 2 && throw(ArgumentError("loop-erased path contains an immediate U-turn"))
        before = after
    end
    return turns
end
