using SafeTestsets

@safetestset "GeneralizedPhase" begin
using GeneralizedPhase
using Test
using AxisKeys
using FFTW
using AxisKeys
using Statistics
using LinearAlgebra

println("Testing sfft")
@testset "keyed sfft" begin
    X = randn(100, 100)
    X = KeyedArray(X, axes(X))
    f, g = GeneralizedPhase.plan_sfft(X)
    F = f(X)
    _F = fftshift(fft(X)) |> Array
    # The NFFT uses some opposite signs, so just check the abs
    @test abs.(real(F)) ≈ abs.(real(_F))
    @test abs.(imag(F)) ≈ abs.(imag(_F))
end;


end
