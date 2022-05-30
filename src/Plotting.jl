using Makie

@recipe(PhaseLines, t, x) do scene
    Theme()
end
plottype(::PhaseLines) = Lines
argument_names(::Type{<: PhaseLines}) = (:t, :x)
function Makie.plot!(p::PhaseLines)
    lines!(p, p[:t], p[:x])
    p
end
