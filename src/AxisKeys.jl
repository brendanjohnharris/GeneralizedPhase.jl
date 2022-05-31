using AxisKeys

function generalized_phase(X::KeyedArray, dim=ndims(X), lp=0.0)
    # We assume the axis keys of the `dim` dimension are equally spaced (`Range`) time points.
    fs = 1/step(axiskeys(X, dim))
    mapslices(x->generalized_phase(x, fs, lp), X; dims=dim)
end
