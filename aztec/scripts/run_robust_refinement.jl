#!/usr/bin/env julia
# One size per invocation, suitable for a six-task bounded Slurm array.
# Benchmark locally with --environments 4; the 2048/size run belongs on Hamilton.
using AztecDiamond, Random, SHA, Printf
include("comparison_geometry.jl")
include("robust_gamma_coins.jl")

filehash(path)=bytes2hex(sha256(read(path)))
function derive_seed(master,n,index,role)
    digest=sha256("$ROBUST_GAMMA_VERSION:$master:$n:$index:$role")
    foldl((x,b)->(x<<8)|UInt64(b),digest[1:8];init=UInt64(0))
end

function main(args)
    iseven(length(args)) || error("arguments are --name value pairs")
    options=Dict(args[i]=>args[i+1] for i in 1:2:length(args))
    all(k in ("--output","--size","--environments","--seed") for k in keys(options)) || error("unknown option")
    output=get(options,"--output","");isempty(output) && error("explicit output required")
    n=parse(Int,get(options,"--size","16"))
    count=parse(Int,get(options,"--environments","4"))
    master=parse(Int,get(options,"--seed","2026100102"))
    n in (16,32,48,64,96,128) && 1<=count<=2048 || error("outside fixed refinement limits")
    mkpath(output)
    # Lock is an atomically created directory, protecting concurrent writers.
    lock=joinpath(output,"writer.lock");mkdir(lock)
    try
        source=join([filehash(p) for p in (@__FILE__,joinpath(@__DIR__,"robust_gamma_coins.jl"),joinpath(@__DIR__,"comparison_geometry.jl"),joinpath(@__DIR__,"../src/AztecDiamond.jl"))],":")
        manifest="{\"version\":\"$ROBUST_GAMMA_VERSION\",\"order\":$n,\"environments\":$count,\"master_seed\":$master,\"source_hashes\":\"$source\",\"julia\":\"$VERSION\",\"probes\":[0.125,0.25]}\n"
        mp=joinpath(output,"manifest.json")
        isfile(mp) && read(mp,String)!=manifest && error("manifest changed; use new output")
        write(mp,manifest)
        header="domain,size,model,law_id,environment_id,probe,axis,r,distance,height_scale,value1,value2,environment_seed,coin_seed"
        temporary=joinpath(output,"measurements.csv.tmp")
        started=time();newpairs=0
        open(temporary,"w") do table
            println(table,header)
            for index in 1:count
                envseed=derive_seed(master,n,index,"environment")
                coinseed=derive_seed(master,n,index,"paired-coins")
                envseed!=coinseed || error("seed collision")
                id="robust-n$(n)-env$(index)"
                case=joinpath(output,"cases",id);mkpath(case)
                files=["a.txt","b.txt","sample1.txt","sample2.txt","rows.csv","metadata.json"]
                certificate=joinpath(case,"sha256.txt")
                if isfile(certificate)
                    actual=join(["$(filehash(joinpath(case,f))) $f" for f in files],"\n")*"\n"
                    read(certificate,String)==actual || error("case checksum mismatch")
                else
                    w=gamma_disordered_weights(Xoshiro(envseed),n;alpha=.2,beta=.25)
                    choices=robust_gamma_choice_pair(Xoshiro(coinseed),w.a,w.b)
                    first=sample_tiling_from_choices(choices.first)
                    second=sample_tiling_from_choices(choices.second)
                    validate_tiling(first).valid && validate_tiling(second).valid || error("invalid tiling")
                    for (file,value) in zip(files[1:4],(w.a,w.b,first,second))
                        write_table(joinpath(case,file),value)
                    end
                    open(joinpath(case,"rows.csv"),"w") do io
                        for rho in (1/8,1/4)
                            r=round(Int,n*rho)
                            law="\"aztec_gamma:alpha=0.2:beta=0.25:scale=1:a-north,b-east\""
                            label=@sprintf("diagonal:rho=%.12g",rho)
                            println(io,join(("aztec",n,"aztec_gamma",law,id,label,"diagonal",r,sqrt(2)*r,1,
                                symmetric_height_increment(first,r),symmetric_height_increment(second,r),envseed,coinseed),','))
                        end
                    end
                    write(joinpath(case,"metadata.json"),"{\"environment_seed\":\"$envseed\",\"coin_seed\":\"$coinseed\",\"fallbacks\":$(choices.fallbacks),\"version\":\"$ROBUST_GAMMA_VERSION\"}\n")
                    hashes=join(["$(filehash(joinpath(case,f))) $f" for f in files],"\n")*"\n"
                    write(certificate*".tmp",hashes);mv(certificate*".tmp",certificate;force=true)
                    newpairs+=1
                end
                write(table,read(joinpath(case,"rows.csv"),String))
                index%64==0 && (println("n=$n completed=$index/$count");flush(stdout))
            end
        end
        target=joinpath(output,"measurements.csv");mv(temporary,target;force=true)
        summary="{\"all_pairs_certified\":true,\"certification_method\":\"weighted shuffling with rational coins for Float64 reduced weights; reduction rounding audited separately\",\"sampler_version\":\"$ROBUST_GAMMA_VERSION\",\"order\":$n,\"environments\":$count,\"new_pairs_this_invocation\":$newpairs,\"seconds\":$(time()-started),\"measurements_sha256\":\"$(filehash(target))\"}\n"
        write(joinpath(output,"summary.json.tmp"),summary)
        mv(joinpath(output,"summary.json.tmp"),joinpath(output,"summary.json");force=true)
        println(summary)
    finally
        rm(lock)
    end
end
main(ARGS)
