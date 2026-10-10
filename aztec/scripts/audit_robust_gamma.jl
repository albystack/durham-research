#!/usr/bin/env julia
# Compare Float64 reduction with two high-precision reference calculations.
# The result measures discrepancies; it is not interval arithmetic for the entire recurrence.
using AztecDiamond, Random, Printf
include("robust_gamma_coins.jl")

function big_reduce(a,b,bits)
    setprecision(BigFloat,bits) do
        k=size(a,1)-1
        x,y=Matrix{BigFloat}(undef,k,k),Matrix{BigFloat}(undef,k,k)
        for j in 1:k,i in 1:k
            x[i,j]=a[i,j]/(a[i,j]+b[i,j])*(a[i+1,j]+b[i+1,j])
            y[i,j]=b[i,j+1]/(a[i,j+1]+b[i,j+1])*(a[i+1,j+1]+b[i+1,j+1])
        end
        return x,y
    end
end

function audit()
    rows=String[]
    for n in (8,32,128)
        w=gamma_disordered_weights(Xoshiro(2026100101+n),n;alpha=.2,beta=.25)
        a,b=reverse(w.a,dims=1),reverse(w.b,dims=1)
        a256,b256=setprecision(BigFloat,256) do;(BigFloat.(a),BigFloat.(b));end
        a512,b512=setprecision(BigFloat,512) do;(BigFloat.(a),BigFloat.(b));end
        maximum_error,total_error,reference_error=0.,0.,0.
        rounded=0;coins=0
        start=time()
        for level in n:-1:1
            setprecision(BigFloat,512) do
                for j in 1:level,i in 1:level
                    p=BigFloat(b[i,j])/(BigFloat(a[i,j])+BigFloat(b[i,j]))
                    reference=b512[i,j]/(a512[i,j]+b512[i,j])
                    other=b256[i,j]/(a256[i,j]+b256[i,j])
                    error=Float64(abs(p-reference))
                    maximum_error=max(maximum_error,error);total_error+=error
                    reference_error=max(reference_error,Float64(abs(reference-other)))
                    naive=b[i,j]/(a[i,j]+b[i,j])
                    rounded+=naive==0 || naive==1
                    coins+=1
                end
            end
            if level>1
                a,b=AztecDiamond.reduce_gamma_weights(a,b)
                a256,b256=big_reduce(a256,b256,256)
                a512,b512=big_reduce(a512,b512,512)
            end
        end
        maximum_error<1e-10 && 2total_error<1e-6 && reference_error<1e-50 || error("numeric audit tolerance exceeded")
        push!(rows,"{\"order\":$n,\"coins_per_replica\":$coins,\"old_rounded_endpoints\":$rounded,\"max_probability_error_vs_512bit\":$maximum_error,\"two_replica_sum_abs_error_proxy\":$(2total_error),\"max_256_vs_512_discrepancy\":$reference_error,\"seconds\":$(time()-start)}")
        println(last(rows));flush(stdout)
    end
    path=isempty(ARGS) ? "docs/validation/numerical_20261001/aztec_numeric_audit.json" : only(ARGS)
    mkpath(dirname(path))
    write(path,"{\"sampler_version\":\"$ROBUST_GAMMA_VERSION\",\"note\":\"Three fixed environments only; high precision comparison, not a global proof\",\"cases\":["*join(rows,",")*"]}\n")
end

audit()
