
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



dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\216_vertices_T_"

temperatures = [0.125, 0.25, 0.5, 0.0625]

my_plot = Plots.plot()


bond_angle_std_vec = zeros(length(temperatures))

for i in eachindex(temperatures)

    structure_factor_dict = GU.load_h5_dict(dict_path*string(temperatures[i])*"_quenched_structure_factor_bartlett_isotrope.h5")

    my_plot = Plots.plot!(structure_factor_dict["wavenumber_vec"], structure_factor_dict["structure_factor_vec"], label = Latex.L"kT="*string(temperatures[i]) )
end

my_plot = Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,27.5), ylims=(0,5), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_bartlett_216_vertices_T_0.125_0.5.png")



temperatures = [1, 2, 4, 0.0625]

my_plot = Plots.plot()


bond_angle_std_vec = zeros(length(temperatures))

for i in eachindex(temperatures)

    structure_factor_dict = GU.load_h5_dict(dict_path*string(temperatures[i])*"_quenched_structure_factor_bartlett_isotrope.h5")

    my_plot = Plots.plot!(structure_factor_dict["wavenumber_vec"], structure_factor_dict["structure_factor_vec"], label = Latex.L"kT="*string(temperatures[i]) )
end

my_plot = Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,27.5), ylims=(0,5), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_bartlett_216_vertices_T_1_4.png")



temperatures = [4, 6, 8, 0.0625]

my_plot = Plots.plot()


bond_angle_std_vec = zeros(length(temperatures))

for i in eachindex(temperatures)

    structure_factor_dict = GU.load_h5_dict(dict_path*string(temperatures[i])*"_quenched_structure_factor_bartlett_isotrope.h5")

    my_plot = Plots.plot!(structure_factor_dict["wavenumber_vec"], structure_factor_dict["structure_factor_vec"], label = Latex.L"kT="*string(temperatures[i]) )
end

my_plot = Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,27.5), ylims=(0,5), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_bartlett_216_vertices_T_4_8.png")