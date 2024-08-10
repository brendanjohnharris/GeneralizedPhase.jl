using GeneralizedPhase
using DSP
using MAT
using HTTP
using CairoMakie

# ? ------------------------------- Downloads ------------------------------ ? #
fp, _ = mktemp()
# HTTP.download("https://github.com/mullerlab/generalized-phase/blob/main/data/exampleData.mat?raw=true", fp)
# X = matread(fp)["exampleData"]
HTTP.download("https://github.com/mullerlab/generalized-phase/blob/main/data/exampleChannel.mat?raw=true", fp)

# ? ------------------------------ Input data ------------------------------ ? #
x = matread(fp)["x"] |> vec
design = Butterworth(4)
fs = 1000
lp = 0
dt = 1/fs
T = length(x)/fs
times = dt:dt:T

# ? -------------------- Preprocessing - lowpass + notch ------------------- ? #
flt(x) = x  |> x -> filtfilt(digitalfilter(Bandpass(5, 200; fs), design), x) #=
       =#   |> x -> filtfilt(digitalfilter(Bandstop(58, 62; fs), design), x) #=
       =#   |> x -> filtfilt(digitalfilter(Bandstop(115, 125; fs), design), x)

y = flt(x)

# ? --------------------------- Generalized phase -------------------------- ? #
ỹ = GeneralizedPhase._phasefilter(y, fs)
𝜑 = generalized_phase(y, fs, lp);

# ? --------------------------------- Plots -------------------------------- ? #
colormap = cyclic #:cyclic_mygbm_30_95_c78_n256_s25 # :cyclic_wrwbw_40_90_c42_n256_s25
f = Figure(size=(1500, 500), backgroundcolor=:transparent)
ax = Axis(f[1, 1], backgroundcolor=:transparent)
lines!(ax, times, y, color=:gray50, linewidth=7, label="Raw signal")
lines!(ax, times, ỹ; color=𝜑, linewidth=7, colormap, label="Filtered signal w/ generalized phase")
xlims!(0.06, 0.65)
hidespines!(ax); hidedecorations!(ax)

ins = Axis(f[1, 1],
            width=Relative(0.3),
            height=Relative(0.3),
            aspect=DataAspect(),
            valign=:top,
            halign=:left,
            backgroundcolor=:transparent)
text!(ins, "GP", align=(:center, :center), fontsize=30, color=:gray50, position=CairoMakie.Point2f0(0, 0))
xlims!(ins, -1.5, 1.5); ylims!(ins, -1.5, 1.5)
hidespines!(ins); hidedecorations!(ins)
t = -pi:0.01:π; _x = cos.(t); _y = sin.(t)
lines!(ins, _x, _y; color=t, linewidth=25, colormap)
file = try; normpath(@__DIR__, "../")*"example.pdf"; catch; mktemp()*".pdf"; end
file = try; normpath(@__DIR__, "../")*"example.png"; catch; mktemp()*".png"; end
save(file, f; px_per_unit=5)
print("Saved figure to "*file)


# ? ------------ Compare generalized phase to temporal filtering ----------- ? #
f = Figure(size=(600, 300), backgroundcolor=:transparent)
ax = Axis(f[2, 1], backgroundcolor=:transparent)
lines!(ax, times, angle.(hilbert(ỹ)); color=(:cornflowerblue, 0.8), linewidth=4, label="Analytic phase after filtering")
lines!(ax, times, generalized_phase(x, fs, lp), color=(:crimson, 0.8), linewidth=4, label="Generalized phase")
xlims!(0.06, 0.65)
hidespines!(ax); hidedecorations!(ax)
L = Legend(f[1, 1], ax,
            width=Relative(0.3),
            height=Relative(0.1),
            backgroundcolor=:transparent,
            framevisible=false,
            labelcolor=:gray50,
            tellwidth=true,
            tellheight=true,
            padding=0.0,
            orientation = :horizontal)
rowsize!(f.layout, 1, Relative(0.05))
file = try; normpath(@__DIR__, "../")*"comparison.png"; catch; mktemp()*".png"; end
save(file, f; px_per_unit=5)
print("Saved figure to "*file)

begin
    using Random
    Random.seed!(2)
ts = 0:0.001:10
x = sin.(0.5π*ts) #.+ (cos.(2π*ts.+0.5).^2)
x = sin.(0.5π*ts) .+ 0.3.*(randn(length(x)))
x = GeneralizedPhase._phasefilter(x, 1000; band=[0.5, 2])
ϕ = angle.(hilbert(x))
# ϕ =  _generalized_phase(x, 1000);
 y = x[6900:7800]
 p = ϕ[6900:7800]
 lines(y)
 lines!(p.*5)
current_figure()

colormap = cyclic #:cyclic_mygbm_30_95_c78_n256_s25 # :cyclic_wrwbw_40_90_c42_n256_s25
f = Figure(size=(1500, 500), backgroundcolor=:transparent)
ax = Axis(f[1, 1], backgroundcolor=:transparent)
lines!(ax, y; color=p.+pi, linewidth=20, colormap, label="Filtered signal w/ generalized phase", linecap=:round)
hidedecorations!(ax)
hidespines!(ax)


ins = Axis(f[1, 1],
            width=Relative(0.3),
            height=Relative(0.3),
            aspect=DataAspect(),
            valign=:top,
            halign=:left,
            backgroundcolor=:transparent)
text!(ins, "GP", align=(:center, :center), fontsize=30, color=:gray50, position=CairoMakie.Point2f0(0, 0))
xlims!(ins, -1.5, 1.5); ylims!(ins, -1.5, 1.5)
hidespines!(ins); hidedecorations!(ins)
t = -pi:0.01:pi; _x = cos.(t); _y = sin.(t)
lines!(ins, _x, _y; color=t, linewidth=25, colormap)
file = try; normpath(@__DIR__, "../")*"example_2.pdf"; catch; mktemp()*".pdf"; end
save(file, f; px_per_unit=5)
f
end
