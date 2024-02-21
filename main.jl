
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import Plots
import LaTeXStrings as Latex # to display latex symbols in plot labels
import NaNStatistics

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105


# path where structures are stored
load_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\without_ring_size_limitation\\"


# load graph with 216 vertices
#graph_dict_2 = NG.load_graph_from_h5_and_MGformat(load_path*"216_vertices_T_1_quenched")
#
## determine structure factor
#structure_factor_dict_2 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
#    graph_dict_2;
#    sampling_distance_step_length = 0.025,
#    maximal_sampling_distance = 4*graph_dict_2["supercell_edge_length"],
#    nr_wavevector_samples = 10000,
#    save_result = true,
#    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\216_vertices_T_1_quenched",
#    label = "216_vertices_T_1_quenched")
#
#
## load graph with 512 vertices
#graph_dict_5 = NG.load_graph_from_h5_and_MGformat(load_path*"512_vertices_T_1_quenched")

# determine structure factor
#structure_factor_dict_5 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
#    graph_dict_5;
#    sampling_distance_step_length = 0.025,
#    maximal_sampling_distance = 4*graph_dict_5["supercell_edge_length"],
#    nr_wavevector_samples = 5000,
#    save_result = true,
#    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\512_vertices_T_1_quenched",
#    label = "512_vertices_T_1_quenched")


## load graph with 1000 vertices
#graph_dict_5 = NG.load_graph_from_h5_and_MGformat(load_path*"1000_vertices_T_1_quenched")
#
## determine structure factor
#structure_factor_dict_1 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
#    graph_dict_1;
#    sampling_distance_step_length = 0.025,
#    maximal_sampling_distance = 4*graph_dict_5["supercell_edge_length"],
#    nr_wavevector_samples = 5000,
#    save_result = true,
#    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\1000_vertices_T_1_quenched",
#    label = "1000_vertices_T_1_quenched")
#
## load graph with 1000 vertices for T=4
graph_dict_1_4 = NG.load_graph_from_h5_and_MGformat(load_path*"1000_vertices_T_4_quenched")

# determine structure factor
structure_factor_dict_1_4 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_1_4;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_5["supercell_edge_length"],
    nr_wavevector_samples = 5000,
    save_result = true,
    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\1000_vertices_T_4_quenched",
    label = "1000_vertices_T_4_quenched")


# plot structure factors 

window_size = length(structure_factor_dict_1["wavenumber_vec"])/50.0

Plots.plot(
    structure_factor_dict_1["wavenumber_vec"],
    NaNStatistics.movmean(structure_factor_dict_1_4["structure_factor_vec"], window_size),
    label = Latex.L"kT=4"
)
Plots.plot!(
    structure_factor_dict_1["wavenumber_vec"],
    NaNStatistics.movmean(structure_factor_dict_1["structure_factor_vec"], window_size),
    label = Latex.L"kT=1"
)
Plots.plot!(
    xlabel = "wavenumber",
    ylabel = "structure factor",
    legend = true,
    xlims=(0,32.5), ylims=(0,10)
)
