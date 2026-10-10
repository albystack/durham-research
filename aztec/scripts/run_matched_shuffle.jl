#!/usr/bin/env julia
# Independent Aztec arm of the frozen October 4 comparison. The Python
# dispatcher reads the JSON plan and passes explicit, recorded arguments.
using AztecDiamond, Random, SHA, Printf
include("comparison_geometry.jl")
include("robust_gamma_coins.jl")

filehash(path) = bytes2hex(sha256(read(path)))
function derive_seed(master, model, n, index, role)
    digest = sha256("matched-shuffle-20261004:$master:aztec:$model:$n:$index:$role")
    foldl((x,b)->(x<<8)|UInt64(b), digest[1:8]; init=UInt64(0))
end

function main(args)
    iseven(length(args)) || error("arguments must be --name value pairs")
    options = Dict(args[i]=>args[i+1] for i in 1:2:length(args))
    Set(keys(options)) == Set(["--output","--size","--model","--environments","--seed","--plan-sha256","--validation"]) || error("wrong options")
    n = parse(Int, options["--size"])
    model = options["--model"]
    count = parse(Int, options["--environments"])
    master = parse(Int, options["--seed"])
    validation = options["--validation"]
    planhash = options["--plan-sha256"]
    n in (8,16,24,32,48,64) && 1<=count<=2048 || error("outside frozen limits")
    model in ("aztec_gamma","aztec_mean","uniform") || error("unknown model")
    validation in ("true","false") || error("invalid validation flag")
    output = options["--output"]
    ispath(output) && error("fresh output directory required")
    mkpath(output)
    write(joinpath(output,"manifest.json"), "{\"domain\":\"aztec\",\"size\":$n,\"model\":\"$model\",\"pairs\":$count,\"master_seed\":$master,\"validation_only\":$validation,\"plan_sha256\":\"$planhash\",\"sampler_version\":\"$ROBUST_GAMMA_VERSION\",\"julia\":\"$VERSION\"}\n")
    started = time()
    # An absent COMPLETE status always blocks analysis, including scheduler kills.
    write(joinpath(output,"status.json"), "{\"state\":\"STARTED\",\"completed\":0}\n")
    temporary = joinpath(output,"measurements.csv.tmp")
    open(temporary,"w") do table
        println(table,"domain,size,model,law_id,environment_id,probe,axis,r,distance,height_scale,value1,value2,environment_seed,coin_seed")
        for index in 1:count
            envseed = derive_seed(master,model,n,index,"environment")
            coinseed = derive_seed(master,model,n,index,"paired-coins")
            envseed != coinseed || error("seed collision")
            if model == "aztec_gamma"
                w = gamma_disordered_weights(Xoshiro(envseed),n;alpha=.2,beta=.25)
                a,b = w.a,w.b
            else
                a,b = model=="aztec_mean" ? (fill(.2,n,n),fill(.25,n,n)) : (ones(n,n),ones(n,n))
            end
            choices = robust_gamma_choice_pair(Xoshiro(coinseed),a,b)
            first,second = sample_tiling_from_choices(choices.first),sample_tiling_from_choices(choices.second)
            validate_tiling(first).valid && validate_tiling(second).valid || error("invalid tiling")
            case = joinpath(output,"cases",@sprintf("env-%05d",index));mkpath(case)
            for (name,value) in zip(("a.txt","b.txt","sample1.txt","sample2.txt"),(a,b,first,second))
                write_table(joinpath(case,name),value)
            end
            id = "aztec-$model-n$n-env$index"
            rows = IOBuffer()
            for rho in (1/8,1/4)
                r = round(Int,n*rho)
                law = "\"$model:alpha=0.2:beta=0.25:scale=1:a-north,b-east\""
                label = @sprintf("diagonal:rho=%.12g",rho)
                println(rows,join(("aztec",n,model,law,id,label,"diagonal",r,sqrt(2)*r,1,
                    symmetric_height_increment(first,r),symmetric_height_increment(second,r),envseed,coinseed),','))
            end
            rowtext = String(take!(rows));write(table,rowtext);write(joinpath(case,"rows.csv"),rowtext)
            write(joinpath(case,"sha256.txt"),join(["$(filehash(joinpath(case,f))) $f" for f in ("a.txt","b.txt","sample1.txt","sample2.txt","rows.csv")],"\n")*"\n")
            if index%64==0
                flush(table)
                write(joinpath(output,"status.json"),"{\"state\":\"STARTED\",\"completed\":$index}\n")
                println("aztec $model n=$n completed=$index/$count");flush(stdout)
            end
        end
    end
    mv(temporary,joinpath(output,"measurements.csv"))
    write(joinpath(output,"status.json"),"{\"state\":\"COMPLETE\",\"completed\":$count,\"seconds\":$(time()-started),\"measurements_sha256\":\"$(filehash(joinpath(output,"measurements.csv")))\",\"arithmetic\":\"Float64 reductions with robust creation coins; rounding audited separately\"}\n")
end
main(ARGS)
