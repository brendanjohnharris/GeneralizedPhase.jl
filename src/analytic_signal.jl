using FFTW

function marple!(x::AbstractArray{<:Complex}; dim=ndims(x))
    fft!(x, dim)
    N = size(x, dim)
    cond(x, m) = (m = m-1;
                    if      m == 0;                 0.0
                    elseif  1 ≤ m ≤ N/2 - 1;        2*x
                    elseif  m == N/2;               x
                    elseif  N/2 + 1 ≤ m ≤ N-1;      0.0
                    end
                 )
    xba
    map!(cond, x, x, getindex.(CartesianIndices(x), dim))
    ifft!(x, dim)
end

analytic_signal!(x::AbstractArray{<:Complex}; alg=marple!, kw...) = alg(x; kw...)
analytic_signal(x::AbstractArray{<:Complex}; kw...) = (y = copy(x); analytic_signal!(y; kw...); y)
analytic_signal(x::AbstractArray; kw...) = analytic_signal(x.|>Complex)
