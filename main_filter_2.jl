
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


import Plots
import Graphs
import MetaGraphsNext
import LinearAlgebra
import Polynomials
import Format
import Statistics
import Random

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_networks\dia\run_2\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\neural_network_networks\dia\\"

structure_factor_dict = GU.load_h5_dict(analysis_data_path*"dia_beta_0.1369_t_max_0.1993_t_gradient_0.9613_structure_factor_bonds_array.h5")

cmax = 2.0

NA.plot_structure_factor_heatmap(
    structure_factor_dict,
    save_path*"dia_beta_0.1369_t_max_0.1993_t_gradient_0.9613_structure_factor_bonds_array_x",
    save_plot = true,
    clims = (0, cmax),
    x_y_lims = nothing,
    wavevector_component_to_fix = 1,
    wavevector_value_fixed = 0)

NA.plot_structure_factor_heatmap(
    structure_factor_dict,
    save_path*"dia_beta_0.1369_t_max_0.1993_t_gradient_0.9613_structure_factor_bonds_array_y",
    save_plot = true,
    clims = (0, cmax),
    x_y_lims = nothing,
    wavevector_component_to_fix = 2,
    wavevector_value_fixed = 0)

NA.plot_structure_factor_heatmap(
    structure_factor_dict,
    save_path*"dia_beta_0.1369_t_max_0.1993_t_gradient_0.9613_structure_factor_bonds_array_z",
    save_plot = true,
    clims = (0, cmax),
    x_y_lims = nothing,
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)