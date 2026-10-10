#!/usr/bin/env julia
# Re-read saved matchings and independently recompute the recorded increments.
using AztecDiamond
length(ARGS)==3 || error("usage: ROOT ORDER COUNT")
root,n,count = ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3])
for index in 1:count
    case=joinpath(root,"cases","env-"*lpad(index,5,'0'))
    rows=readlines(joinpath(case,"rows.csv"))
    length(rows)==2 || error("wrong probe count")
    tilings=BitMatrix[]
    for replica in 1:2
        lines=readlines(joinpath(case,"sample$replica.txt"))
        values=[parse.(Int,split(line)) for line in lines]
        length(values)==2n && all(length(v)==2n for v in values) || error("wrong tiling shape")
        all(all(x in (0,1) for x in v) for v in values) || error("nonbinary tiling")
        tiling=BitMatrix(reduce(vcat,permutedims.(values)))
        validate_tiling(tiling).valid || error("invalid saved matching")
        push!(tilings,tiling)
    end
    for (line,rho) in zip(rows,(1/8,1/4))
        # The four final CSV columns are numeric (the quoted law is earlier).
        fields=split(line,',')
        observed=(parse(Int,fields[end-3]),parse(Int,fields[end-2]))
        expected=Tuple(symmetric_height_increment(t,round(Int,n*rho)) for t in tilings)
        observed==expected || error("saved observable mismatch at environment $index")
    end
end
println("Verified $count saved Aztec pairs at n=$n")
