using AxisKeys

function _generalized_phase(X::KeyedArray, dim=ndims(X), lp=0.0)
    # We assume the axis keys of the `dim` dimension are equally spaced (`Range`) time points.
    @assert axiskeys(X, dim) isa AbstractRange
    axs = axiskeys(X)
    fs = 1/step(axiskeys(X, dim))
    KeyedArray(mapslices(x->_generalized_phase(x, fs, lp), X; dims=dim), axs)
end

function generalized_phase(X::KeyedArray, dim=ndims(X), lp=0.0)
    @assert axiskeys(X, dim) isa AbstractRange
    axs = axiskeys(X)
    fs = 1/step(axiskeys(X, dim))
    KeyedArray(mapslices(x->generalized_phase(x, fs, lp), X; dims=dim), axs)
end
