
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex
import NaNStatistics


path = raw"C:\Users\HemmannF\switchdrive\presentations\material\\"

fontsize=18

Plots.gr()
Plots.default(grid=false, 
legend = true, 
dpi=250,
xtickfontsize=fontsize,
ytickfontsize=fontsize,
xguidefontsize=fontsize,
yguidefontsize=fontsize,
legendfontsize=fontsize,
bottom_margin = 3Plots.mm,
linewidth=3, 
thickness_scaling = 1)

# functions to have pi ticks
function pitick(start, stop, denom; mode=:text)
    a = Int(cld(start, 2*π/denom))
    b = Int(fld(stop, 2*π/denom))
    tick = range(a*2*π/denom, b*2*π/denom; step=2*π/denom)
    ticklabel = piticklabel.( 2 .* (a:b) .// denom, Val(mode))
    tick, ticklabel
end

function piticklabel(x::Rational, ::Val{:text})
    iszero(x) && return "0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return S * N * "π"
    S * N * "π/" * repr(d)
end

function piticklabel(x::Rational, ::Val{:latex})
    iszero(x) && return Latex.L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return Latex.L"%$S%$N\pi"
    Latex.L"%$S\frac{%$N\pi}{%$d}"
end


dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\without_ring_size_limitation\\"

function heaviside(t)
    0.5 * (sign(t) + 1)
 end
 
temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

temperature = temperatures[1]

mc_step_vec = collect(0:0.1:20)
temperature_vec = (heaviside.(- (mc_step_vec .- 2)) .* temperature 
.+ (0.3 .- 0.1 .* mc_step_vec ).* heaviside.(- (mc_step_vec .- 3)) .* heaviside.( (mc_step_vec .- 2)))


Plots.plot(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}="*string(temperature), ls= :dash)

temperature = temperatures[3]
temperature_vec = (heaviside.(- (mc_step_vec .- 2)) .* temperature 
.+ (0.35 .- 0.1 .* mc_step_vec) .* heaviside.(- (mc_step_vec .- 3.5)) .* heaviside.( (mc_step_vec .- 2)))

Plots.plot!(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}="*string(temperature), ls= :dash )

Plots.plot!(xlabel="Monte Carlo step", ylabel=Latex.L"kT", xlims=(0, 20), right_margin = 3Plots.mm)

Plots.savefig(path*"temperature_T_0.1_0.15_cool_0.1_per_mc_quenched.png")
