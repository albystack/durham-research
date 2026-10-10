using Test
using AztecDiamond, Random
include("../scripts/comparison_geometry.jl")

@testset "Physical matched-probe geometry" begin
    for n in (4, 8, 16), r in 1:(n-1)
        columns = symmetric_face_columns(n, r)
        left = comparison_face_xy(n, n+1, columns.left)
        right = comparison_face_xy(n, n+1, columns.right)
        @test (right.x-left.x, right.y-left.y) == (r, -r)
    end
end

@testset "Order-two shuffle distribution equals red/yellow Gibbs law" begin
    # Exhaust every creation-coin sequence, aggregate paths reaching the same
    # tiling, and compare with independently computed occupied-edge products.
    a = [0.13 1.7; 2.1 0.42]
    b = [0.9 0.21; 1.3 0.73]
    probabilities = comparison_gamma_probabilities(a,b)
    distribution = Dict{String,Float64}()
    tilings = Dict{String,Matrix{Bool}}()
    for mask in 0:31
        choices = [falses(1,1), falses(2,2)]
        probability = 1.0
        bit = 0
        for level in 1:2, index in eachindex(choices[level])
            yes = (mask >> bit) & 1 == 1
            choices[level][index] = yes
            p = probabilities[level][index]
            probability *= yes ? p : 1-p
            bit += 1
        end
        tiling = sample_tiling_from_choices(choices)
        key = join(Int.(vec(tiling)))
        distribution[key] = get(distribution,key,0.0) + probability
        tilings[key] = tiling
    end
    @test length(distribution) == 8
    weights = comparison_edge_weights(a,b)
    products = Dict(key => prod(weights[tiling]) for (key,tiling) in tilings)
    normalizer = sum(values(products))
    for key in keys(distribution)
        @test distribution[key] ≈ products[key]/normalizer atol=1e-13
        @test validate_tiling(tilings[key]).valid
    end
end

@testset "Corrected Gamma recurrence agrees with generic weighted reduction" begin
    for n in 1:8
        w=gamma_disordered_weights(Xoshiro(731+n),n)
        expected=creation_probabilities(comparison_edge_weights(w.a,w.b))
        actual=comparison_gamma_probabilities(w.a,w.b)
        for k in 1:n
            @test actual[k] ≈ expected[k] atol=1e-13
        end
        pair=comparison_gamma_choice_pair(Xoshiro(71),w.a,w.b)
        @test validate_tiling(sample_tiling_from_choices(pair.first)).valid
        @test validate_tiling(sample_tiling_from_choices(pair.second)).valid
    end
end

@testset "Independent order-two/order-three weighted distribution and height audit" begin
    include("../scripts/audit_comparison_mapping.jl")
    mktempdir() do directory
        path=joinpath(directory,"audit.json")
        audit_mapping(path) # Raises on any incorrect distribution or height path.
        @test isfile(path)
    end
end
