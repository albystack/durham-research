#!/usr/bin/env julia
# Exhaustive finite tests of the *distribution*, not just recurrence agreement.
using AztecDiamond, Random, Printf
if !isdefined(@__MODULE__, :comparison_gamma_probabilities)
    include("comparison_geometry.jl")
end

function exact_shuffle_distribution(probabilities)
    n=length(probabilities)
    bits=sum(k*k for k in 1:n)
    bits <= 14 || error("audit is restricted to orders <=3")
    distribution=Dict{Tuple,Float64}()
    tilings=Dict{Tuple,Matrix{Bool}}()
    for mask in 0:(2^bits-1)
        choices=[falses(k,k) for k in 1:n]
        probability=1.0
        bit=0
        for k in 1:n, index in eachindex(choices[k])
            yes=!iszero(mask & (1<<bit))
            choices[k][index]=yes
            p=probabilities[k][index]
            probability *= yes ? p : 1-p
            bit+=1
        end
        tiling=sample_tiling_from_choices(choices)
        key=Tuple(vec(tiling))
        distribution[key]=get(distribution,key,0.0)+probability
        tilings[key]=tiling
    end
    return distribution,tilings
end

function audit_mapping(output)
    rows=String[]
    for n in (2,3)
        # Moderate non-symmetric weights avoid near-degenerate probabilities.
        rng=Xoshiro(929+n)
        a=0.1 .+ 2rand(rng,n,n)
        b=0.1 .+ 2rand(rng,n,n)
        weights=comparison_edge_weights(a,b)
        generic,tilings=exact_shuffle_distribution(creation_probabilities(weights))
        legacy,_=exact_shuffle_distribution(gamma_disordered_probabilities(a,b))
        corrected,_=exact_shuffle_distribution(comparison_gamma_probabilities(a,b))
        products=Dict(key=>prod(weights[tiling]) for (key,tiling) in tilings)
        total=sum(values(products))
        gibbs=Dict(key=>value/total for (key,value) in products)
        residual(d)=maximum(abs(d[key]-gibbs[key]) for key in keys(gibbs))
        tv(d)=sum(abs(d[key]-gibbs[key]) for key in keys(gibbs))/2
        residual(generic)<1e-12 && residual(corrected)<1e-12 || error("corrected distribution validation failed")
        length(tilings)==2^(n*(n+1)÷2) || error("incomplete tiling enumeration")
        # Independent local dual path: east then south for each diagonal step.
        # At fixed central row its fluctuating height increment has unit scale.
        offsets=Set{Int}()
        for tiling in values(tilings)
            columns=symmetric_face_columns(n,1)
            face=comparison_face_xy(n,n+1,columns.left)
            occupied=Set((j-i,2n+1-i-j) for i in 1:2n for j in 1:2n if tiling[i,j])
            x,y=face.x,face.y
            # Start (2x,2y), cross east at (2x+1,2y), then south at (2x+2,2y-1).
            east=((2x+1,2y) in occupied) ? 1 : 0
            south=((2x+2,2y-1) in occupied) ? 1 : 0
            # Both crossings border the same primal vertex. Their signs agree;
            # Julia global sign depends on order parity; magnitude is one.
            geometric=(iseven(x+y) ? 1 : -1)*(4south+4east-2)
            push!(offsets,symmetric_height_increment(tiling,1)-(iseven(n) ? -1 : 1)*geometric)
        end
        length(offsets)==1 || error("height coordinate/normalisation mismatch: $offsets")
        push!(rows,"{\"order\":$n,\"tilings\":$(length(tilings)),\"legacy_max_probability_error\":$(residual(legacy)),\"legacy_total_variation\":$(tv(legacy)),\"generic_max_probability_error\":$(residual(generic)),\"corrected_max_probability_error\":$(residual(corrected)),\"height_path_deterministic_offset\":$(only(offsets))}")
    end
    mkpath(dirname(output))
    write(output,"{\"description\":\"Exact creation-path enumeration versus occupied red/yellow edge Gibbs products; historical routines unchanged\",\"cases\":["*join(rows,",")*"]}\n")
    println(read(output,String))
end

if abspath(PROGRAM_FILE)==@__FILE__
    audit_mapping(isempty(ARGS) ? "docs/validation/aztec_mapping_audit.json" : only(ARGS))
end
