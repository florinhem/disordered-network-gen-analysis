
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



l_max = 12

some_vertex = 10

y_vec = collect(0:l_max)

dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_low_t = NG.load_graph_from_h5_and_MGformat(dict_path*
    "1000_vertices_T_0.1_heated_for_0.5_steps_quenched")

graph_dict_high_t = NG.load_graph_from_h5_and_MGformat(dict_path*
    "1000_vertices_T_1.0_heated_for_0.5_steps_quenched")

graph_dict_diamond = NG.load_graph_from_h5_and_MGformat(dict_path*
    "216_vertices_T_0.1_heated_for_0.01_steps_quenched")

q_l_total_network_mean_dict = NA.get_q_l_total_network_mean_dict(graph_dict_high_t, l_max)
q_l_total_network_mean_vec = NA.convert_q_l_dict_to_vec(q_l_total_network_mean_dict, l_max)

Plots.scatter(y_vec, Measurements.value.(q_l_total_network_mean_vec),
        yerr=Measurements.uncertainty.(q_l_total_network_mean_vec), linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=1.0", markerstrokecolor=Plots.palette(:tab10)[2] , markercolor=Plots.palette(:tab10)[2], markersize=6, markerstrokewidth =2)


q_l_total_network_mean_dict = NA.get_q_l_total_network_mean_dict(graph_dict_low_t, l_max)
q_l_total_network_mean_vec = NA.convert_q_l_dict_to_vec(q_l_total_network_mean_dict, l_max)

Plots.scatter!(y_vec, Measurements.value.(q_l_total_network_mean_vec),
    yerr=Measurements.uncertainty.(q_l_total_network_mean_vec), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0.1", markerstrokecolor=Plots.palette(:tab10)[1], markercolor=Plots.palette(:tab10)[1], markersize=6, markerstrokewidth =2)


q_l_total_network_mean_dict = NA.get_q_l_averaged_single_vertex_dict(graph_dict_diamond, some_vertex, l_max)
q_l_total_network_mean_vec = NA.convert_q_l_dict_to_vec(q_l_total_network_mean_dict, l_max)

Plots.scatter!(y_vec, Measurements.value.(q_l_total_network_mean_vec), linecolor=Plots.palette(:tab10)[3], label = "diamond", markerstrokecolor=Plots.palette(:tab10)[3], markercolor=Plots.palette(:tab10)[3], markersize=6, markerstrokewidth =2)
    
Plots.plot!(xlabel=Latex.L"l", 
    ylabel=Latex.L"\overline{q}_l", ylims=(0,1.15), xticks=0:2:12)


Plots.savefig(path*"q_l_total_network_mean_1000_vertices_T_0.1_1.0_diamond.png")