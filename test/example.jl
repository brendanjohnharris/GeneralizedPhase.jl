using GeneralizedPhase
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
colormap = :cyclic_mygbm_30_95_c78_n256_s25 # :cyclic_wrwbw_40_90_c42_n256_s25
f = Figure(resolution=(1100, 500))
ax = Axis(f[1, 1])
lines!(ax, times, y, color=:black, linewidth=7, label="Raw signal")
lines!(ax, times, ỹ; color=𝜑, linewidth=7, colormap, label="Filtered signal w/ generalized phase")
xlims!(0.06, 0.65)
hidespines!(ax); hidedecorations!(ax)

ins = Axis(f[1, 1],
            width=Relative(0.3),
            height=Relative(0.3),
            aspect=DataAspect(),
            valign=:top,
            halign=:left)
text!(ins, "GP", align=(:center, :center), textsize=30, color=:black, position=CairoMakie.Point2f0(0, 0))
xlims!(ins, -1.5, 1.5)
ylims!(ins, -1.5, 1.5)
hidespines!(ins); hidedecorations!(ins)
t = -0.01:0.01:2π
_x = cos.(t)
_y = sin.(t)
lines!(ins, _x, _y; color=t, linewidth=25, colormap)
f
file = mktemp()[1]*".pdf"
save(file, f)
print("Saved figure to "*file)
