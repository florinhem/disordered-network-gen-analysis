
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex
import Measurements
import Polynomials


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


dict_path = raw"..\analysis_data\random_networks\\"
filename = "216_vertices_T_0.2_heated_for_0.5_steps_quenched"

structure_factor_angle_averaged_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_angle_averaged.h5")

hyperuniformity_metric, polynomial_fit = NA.get_hyperuniformity_metric(structure_factor_angle_averaged_dict)

Plots.plot(structure_factor_angle_averaged_dict["wavenumber_vec"], 
                    Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]) , 
                    ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict["structure_factor_vec"]))
Plots.plot!(Polynomials.Polynomial(Measurements.value.( collect(polynomial_fit) )))


Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "structure factor", xtick=pitick(0, 56, 1; mode=:latex), legend = false, xlims=(0, 53), ylims=(0, 2))

Plots.savefig(raw"..\plots\random_networks\\"*filename*".png")