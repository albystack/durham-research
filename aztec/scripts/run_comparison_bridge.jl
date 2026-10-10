#!/usr/bin/env julia
# Small Aztec companion to the square feasibility pilot. No scheduler calls.
# julia --project=aztec aztec/scripts/run_comparison_bridge.jl --output DIR
using AztecDiamond, Random, Printf, SHA
include("comparison_geometry.jl")

function main(args)
    iseven(length(args)) || error("arguments are --name value pairs")
    options = Dict(args[i] => args[i+1] for i in 1:2:length(args))
    allowed = Set(["--output", "--sizes", "--environments", "--seed"])
    all(key in allowed for key in keys(options)) || error("unknown argument")
    output = get(options,"--output","aztec/output/comparison_bridge")
    sizes = parse.(Int,split(get(options,"--sizes","8,16,24,32,48"),','))
    count = parse(Int,get(options,"--environments","4"))
    seed = parse(Int,get(options,"--seed","2026092905"))
    !isempty(sizes) && length(unique(sizes))==length(sizes) || error("distinct sizes required")
    all(n -> iseven(n) && 4 <= n <= 128,sizes) || error("bridge restricted to even orders 4..128")
    1 <= count <= 64 && length(sizes)*count <= 1024 || error("small bridge only")
    seed >= 0 || error("seed must be nonnegative")
    mkpath(output)
    destination=joinpath(output,"measurements.csv")
    isfile(destination) && error("output exists; choose a new directory to preserve provenance")
    temporary=destination*".tmp"
    # Distinct per-environment seeds are stored explicitly. Julia RNG streams
    # are not claimed bitwise compatible with the Python PCG64 streams.
    rng=Xoshiro(seed)
    seeds=rand(rng,UInt64,length(sizes)*count)
    length(unique(seeds))==length(seeds) || error("seed collision")
    index=0
    open(temporary,"w") do io
        println(io,"domain,size,model,law_id,environment_id,probe,axis,r,distance,height_scale,start_x,start_y,end_x,end_y,value1,value2,environment_seed")
        for n in sizes, environment in 1:count
            index+=1
            envseed=seeds[index]
            erng=Xoshiro(envseed)
            weights=gamma_disordered_weights(erng,n;alpha=0.2,beta=0.25)
            choices=comparison_gamma_choice_pair(erng,weights.a,weights.b)
            first=sample_tiling_from_choices(choices.first)
            second=sample_tiling_from_choices(choices.second)
            validate_tiling(first).valid && validate_tiling(second).valid || error("invalid tiling")
            case=joinpath(output,"cases","n$(n)-env$(environment)")
            mkpath(case)
            write_table(joinpath(case,"a.txt"),weights.a)
            write_table(joinpath(case,"b.txt"),weights.b)
            write_table(joinpath(case,"sample1.txt"),first)
            write_table(joinpath(case,"sample2.txt"),second)
            seen=Set{Int}()
            for rho in (1/32,1/16,1/8,1/4)
                r=floor(Int,n*rho+0.5)
                (0<r<n && !(r in seen)) || continue
                push!(seen,r)
                columns=symmetric_face_columns(n,r)
                left=comparison_face_xy(n,n+1,columns.left)
                right=comparison_face_xy(n,n+1,columns.right)
                x1=symmetric_height_increment(first,r)
                x2=symmetric_height_increment(second,r)
                label=@sprintf("diagonal:rho=%.12g",rho)
                # Quote the law field because it contains a comma.
                law="\"aztec_gamma:alpha=0.2:beta=0.25:scale=1:a-north,b-east\""
                println(io,join(("aztec",n,"aztec_gamma",law,"n$(n)-env$(environment)",label,
                    "diagonal",r,sqrt(2)*r,1,left.x,left.y,right.x,right.y,x1,x2,envseed),','))
            end
        end
    end
    mv(temporary,destination)
    open(joinpath(output,"summary.json"),"w") do io
        println(io,"{\"all_pairs_certified\":true,\"certification_method\":\"exact weighted shuffling, not CFTP\",\"purpose\":\"small-order bridge only\",\"master_seed\":$seed,\"environments\":$index,\"julia_version\":\"$VERSION\",\"measurements_sha256\":\"$(bytes2hex(sha256(read(destination))))\"}")
    end
    println("Saved $index paired Aztec environments to $output; no scaling claim.")
end

main(ARGS)
