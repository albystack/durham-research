using Test, Random, AztecDiamond
include("../scripts/robust_gamma_coins.jl")
mutable struct ScriptedWords <: AbstractRNG
    words::Vector{UInt64}
end
Random.rand(r::ScriptedWords,::Type{UInt64})=popfirst!(r.words)

@testset "Rational rare coins and enclosing intervals" begin
    # Both rare tails remain accessible beyond 53 bits; deterministic adversarial streams.
    onesword=typemax(UInt64)
    @test robust_gamma_coin(ScriptedWords([onesword,onesword]),2.0^-100,1.0)==(false,true)
    @test robust_gamma_coin(ScriptedWords([0,0]),1.0,2.0^-100)==(true,true)
    @test robust_gamma_coin(ScriptedWords([0]),2.0^-100,1.0)==(true,false)
    @test robust_gamma_coin(ScriptedWords([onesword]),1.0,2.0^-100)==(false,false)
    for a in (nextfloat(0.),2.0^-100,0.2,1.,1e200,floatmax(Float64)),
        b in (nextfloat(0.),2.0^-100,0.25,1.,1e200,floatmax(Float64))
        lo,hi=coin_bounds(a,b)
        p=Rational{BigInt}(b)/(Rational{BigInt}(a)+Rational{BigInt}(b))
        @test Rational{BigInt}(lo)<=p<=Rational{BigInt}(hi)
    end
    # Enumerate a complete finite uniform-bit partition for a dyadic coin.
    outcomes=[rational_coin_prefix(ScriptedWords(UInt64[]),big(3),big(8),big(k),3) for k in 0:7]
    @test sum(outcomes)==3
end

@testset "New robust paired choices are reproducible and valid" begin
    w=gamma_disordered_weights(Xoshiro(2),8;alpha=0.2,beta=0.25)
    a=robust_gamma_choice_pair(Xoshiro(3),w.a,w.b)
    b=robust_gamma_choice_pair(Xoshiro(3),w.a,w.b)
    @test a==b
    @test validate_tiling(sample_tiling_from_choices(a.first)).valid
    @test validate_tiling(sample_tiling_from_choices(a.second)).valid
end
