
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex
import Measurements
import Polynomials

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

path = raw"..\..\presentations\material\\"


function heat_cool_temperature_vec(x, max_temp, heat_cool_rate)

    if x < max_temp/heat_cool_rate 
        return heat_cool_rate * x
    elseif x < 2*max_temp/heat_cool_rate 
        return 2* max_temp - heat_cool_rate*x
    else
        return 0
    end
end

x_vec = collect(0:0.01:10)
y_vec = heat_cool_temperature_vec.(x_vec, 0.1, 0.1)

min_temp = 0.1
max_temp = 0.5
temperatures = [0.1, 0.125, 0.15,  0.2, 0.25, 0.3, 0.4, 0.5]
normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

Plots.plot(x_vec, heat_cool_temperature_vec.(x_vec, 0.4, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[7])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.25, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[5])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.2, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[4])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.15, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[3])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.1, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[1])
Plots.plot!(legend = false)

Plots.savefig(path*"heat_cool_0.1_temperature_profile.png")