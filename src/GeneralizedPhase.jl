module GeneralizedPhase

export generalized_phase, _generalized_phase

using Requires
function __init__()
    @require AxisKeys="94b1ba4f-4ee9-5380-92f1-94cde586c3c5" include("AxisKeys.jl")
end

using Statistics
using Dierckx
import DSP: hilbert, Bandpass, Butterworth, digitalfilter, filtfilt, Unwrap.unwrap!

rewrap(xp) = @. xp - 2*π*floor((xp-π)/(2*π)) - 2*π

interp1(x, y; k=3, bc="extrapolate", kw...) = Spline1D(x, y; k, bc, kw...) # 3rd order spline instead of hermite
interp1(y; kw...) = interp1(axes(x)..., y; kw...)

naninterp(y; kw...) = (idxs = .!isnan.(y); (interp1(findall(idxs), y[idxs]; kw...), idxs))
naninterp!(y; kw...) = ((f, idxs) = naninterp(y; kw...); y[.!idxs] .= f.(findall(.!idxs)); ())

# The generalized-phase (matlab) package manually implements the Hilbert transform using the Marple method, since the `hilbert` matlab function must be locked inside the signal processing toolbox. The `hilbert` method from DSP does this already.
analytic_signal(x::AbstractArray{<:Real}, dim::Int=ndims(x)) = mapslices(hilbert, x; dims=dim)

# Instantaneous freq. by taking the angle of the product of current and previous 𝑠
ifreq!(𝑠::AbstractVector, 𝛥𝑡) = (𝑠[1:end-1] .= angle.(𝑠[2:end] .* conj.(𝑠[1:end-1]))./(2*π*𝛥𝑡))
ifreq(x, args...) = (y = deepcopy(x); ifreq!(y, args...); y[1:end-1] |> real)

function bwlabel(x::AbstractVector)
    x = x .|> Bool
    L = Vector{Int}(undef, length(x))
    l = 0
    L[1] = x[1] == 1 ? l : 0
    for i in 2:lastindex(x)
        L[i] = x[i] == 1 ? (l += !x[i-1]; l) : 0
    end
    return L
end

nanunwrap!(𝜑) = (idxs = .!isnan.(𝜑); 𝜑[idxs] .= unwrap!(𝜑[idxs]))
rewrap!(𝜑) = (𝜑 .= mod.(𝜑 .+ π, 2*π) .- π)

function _generalized_phase(x::AbstractVector, fs, lp=0.0)
    all(isnan.(x)) && return x
    nwin = 3; # Sets the 'buffer' to interpolate over after neg. freq. periods, in terms of the length of the neg. freq. window.
    𝛥𝑡 = 1/fs
    𝑠 = hilbert(x)
    𝜑 = angle.(𝑠)
    𝜔 = ifreq(𝑠, 𝛥𝑡)
    dir = filter(!isnan, 𝜔) |> mean |> sign
    idx = dir.*(𝜔) .< lp
    idx[1] = false
    L = bwlabel(idx)
    for k in 1:maximum(L)
        idxs = findall(L .== k)
        idx[idxs[1]:(idxs[1] + (idxs[end] - idxs[1])*nwin)] .= true
    end
    nanunwrap!(𝜑)
    𝜑[idx] .= NaN
    naninterp!(𝜑)
    rewrap!(𝜑)
    return 𝜑
end

_phasefilter(x, fs; band=[5, 40]) = filtfilt(digitalfilter(Bandpass(band...; fs), Butterworth(4)), x)
generalized_phase(x, fs, args...; kw...) = _generalized_phase(_phasefilter(x, fs; kw...), fs, args...)

end # module
