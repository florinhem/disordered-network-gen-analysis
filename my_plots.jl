
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


graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\run_1\\"

all_filenames = readdir(graph_dict_path)
filenames = filter(filename -> endswith(filename, "_evolution.h5"), all_filenames)
final_energy_vec = Float64[]

for filename in filenames
    println(filename)

    evolution_dict = GU.load_h5_dict(graph_dict_path*filename)

    push!(final_energy_vec, evolution_dict["total_energy_vec"][end])
end

filenames_sorted = filenames[sortperm(final_energy_vec)   ]
sort!(final_energy_vec)


evolution_dict_high = GU.load_h5_dict(graph_dict_path*filenames_sorted[end])
evolution_dict_middle = GU.load_h5_dict(graph_dict_path*filenames_sorted[13])
evolution_dict_low = GU.load_h5_dict(graph_dict_path*filenames_sorted[1])

Plots.plot(collect(1:length(evolution_dict_high["total_energy_vec"]))./(216*18), evolution_dict_high["total_energy_vec"] ./ 216, ls = :dot, label = "high final energy")
Plots.plot!(collect(1:length(evolution_dict_middle["total_energy_vec"]))./(216*18), evolution_dict_middle["total_energy_vec"] ./ 216, ls = :dot, label = "middle final energy")
Plots.plot!(collect(1:length(evolution_dict_low["total_energy_vec"]))./(216*18), evolution_dict_low["total_energy_vec"] ./ 216, ls = :dot, label = "low final energy", yaxis=:log10)

Plots.plot!(xlabel = "Monte Carlo step", ylabel = "energy per vertex", size=(600, 400), legend=:bottomright)
Plots.savefig(raw"..\plots\random_networks\\216_vertices_anneal_quench_total_energy.png")