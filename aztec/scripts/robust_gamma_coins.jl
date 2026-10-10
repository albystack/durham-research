# Numerically safe coins for the represented reduced weights. Historical code is unchanged.
# The Float64 reduction itself remains approximate and is audited separately.
using Random, AztecDiamond
const ROBUST_GAMMA_VERSION = "gamma-interval-coins-v1"

"""Outward bounds for b/(a+b), using directed neighbouring Float64 values.

Each correctly rounded sum/division lies between its adjacent floats. If a sum
is infinite or the lower denominator bound is zero, the interval widens to
[0,1]. There is never a decision based solely on a rounded probability of one.
"""
function coin_bounds(a::Float64,b::Float64)
    isfinite(a) && isfinite(b) && a>0 && b>0 || error("positive finite weights required")
    s=a+b
    isfinite(s) || return (0.0,1.0)
    slo,shi=prevfloat(s),nextfloat(s)
    lo=max(0.0,prevfloat(b/shi))
    hi=slo>0 ? min(1.0,nextfloat(b/slo)) : 1.0
    return (lo,hi)
end

"""Resolve U<p with a lazily extended uniform binary prefix and exact integers."""
function rational_coin_prefix(rng,n::BigInt,d::BigInt,prefix::BigInt,bits::Int)
    0<=n<=d && d>0 && 0<=prefix<(big(1)<<bits) || error("invalid rational coin")
    while true
        threshold=n<<bits
        (prefix+1)*d<=threshold && return true
        prefix*d>=threshold && return false
        prefix=(prefix<<64)+BigInt(rand(rng,UInt64))
        bits+=64
    end
end

"""Return (coin, exact_fallback_used); correct for stored Float64 a,b."""
function robust_gamma_coin(rng,a::Float64,b::Float64,bounds=coin_bounds(a,b))
    prefix=rand(rng,UInt64)>>11
    lower=ldexp(Float64(prefix),-53)
    upper=ldexp(Float64(prefix+1),-53)
    lo,hi=bounds
    upper<=lo && return (true,false)
    lower>=hi && return (false,false)
    ar,br=Rational{BigInt}(a),Rational{BigInt}(b)
    p=br/(ar+br)
    return (rational_coin_prefix(rng,numerator(p),denominator(p),BigInt(prefix),53),true)
end

"""Corrected orientation, shared reduction, two independent robust coin streams.

New RNG consumption/version: never substitutes silently into historical runs.
"""
function robust_gamma_choice_pair(rng,a,b)
    order=AztecDiamond.check_gamma_weights(a,b)
    sets=[Vector{BitMatrix}(undef,order) for _ in 1:2]
    ca,cb=reverse(Float64.(a),dims=1),reverse(Float64.(b),dims=1)
    na,nb=similar(ca),similar(cb)
    fallbacks=0
    for level in order:-1:1
        first,second=falses(level,level),falses(level,level)
        for j in 1:level,i in 1:level
            bounds=coin_bounds(ca[i,j],cb[i,j])
            first[i,j],f1=robust_gamma_coin(rng,ca[i,j],cb[i,j],bounds)
            second[i,j],f2=robust_gamma_coin(rng,ca[i,j],cb[i,j],bounds)
            fallbacks+=f1+f2
        end
        sets[1][level]=reverse(first,dims=1)
        sets[2][level]=reverse(second,dims=1)
        if level>1
            AztecDiamond.reduce_gamma_weights!(na,nb,ca,cb,level)
            ca,na=na,ca;cb,nb=nb,cb
        end
    end
    return (first=sets[1],second=sets[2],fallbacks=fallbacks)
end
