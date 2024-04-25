
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex
import Measurements


path = raw"..\..\presentations\material\\"

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
thickness_scaling = 1,
framestyle = :box)

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



dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

function heaviside(t)
    0.5 * (sign(t) + 1)
 end

temperature = 1.0

mc_step_vec = collect(0:0.01:20)
temperature_vec = heaviside.(- (mc_step_vec .- 0.5)) .* temperature

Plots.plot(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}=1.0", ls= :dash, linecolor=Plots.palette(:tab10)[2])

temperature = 0.1

temperature_vec = heaviside.(- (mc_step_vec .- 0.5)) .* temperature

Plots.plot!(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}=0.1", ls= :dash, linecolor=Plots.palette(:tab10)[1] )

Plots.plot!(xlabel="Monte Carlo step", ylabel=Latex.L"kT", xlims=(0, 20), right_margin = 3Plots.mm)

Plots.savefig(path*"temperature_T_0.1_1.0_quenched.png")
