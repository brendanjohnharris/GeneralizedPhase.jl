using SafeTestsets

@safetestset "GeneralizedPhase" begin
using GeneralizedPhase
using Test
using CairoMakie


@testset "Muller example" begin
    include("./example.jl")
end;




end
