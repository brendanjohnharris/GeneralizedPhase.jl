using SafeTestsets

@safetestset "GeneralizedPhase" begin
using GeneralizedPhase
using Test

@testset "Analytic signal" begin
    x = LinRange(0, 20, 10000)
    _y = sin.(5.0.*x) + sin.(6.0.*x)
    y = copy(_y) .|> Complex
    analytic_signal!(y)

    lines(_y)
    lines(angle.(y))
    .....move hilber transform stuff to here......
    ..................compare to hilbert transform........................
end

@testset "Muller example" begin
    include("./example.jl")
end;




end
