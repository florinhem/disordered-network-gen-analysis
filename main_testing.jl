
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import Format
import MetaGraphsNext
import SphericalHarmonics


network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\random\ctn\run_3\\"

filename = "ctn_beta_3.8161_t_max_4.5843_t_gradient_5.4108"

spatial_network = NG.load_spatial_network_from_gml(
    network_path*filename*".gml")

bond_length_std, bond_length_vec = NA.get_bond_length_std(
    spatial_network)

spatial_network_1 = NG.randomly_displace_all_vertices!(
    deepcopy(spatial_network);
    sigma = 0.1,
    update_total_energy = true)

spatial_network_2 = NG.randomly_displace_all_vertices!(
    deepcopy(spatial_network);
    sigma = 0.2,
    update_total_energy = true)

spatial_network_3 = NG.randomly_displace_all_vertices!(
    deepcopy(spatial_network);
    sigma = 0.3,
    update_total_energy = true)

spatial_network_4 = NG.randomly_displace_all_vertices!(
    deepcopy(spatial_network);
    sigma = 0.4,
    update_total_energy = true)

bond_length_std, bond_length_vec = NA.get_bond_length_std(
    spatial_network)
bond_length_std_1, bond_length_vec_1 = NA.get_bond_length_std(
    spatial_network_1)
bond_length_std_2, bond_length_vec_2 = NA.get_bond_length_std(
    spatial_network_2)
bond_length_std_3, bond_length_vec_3 = NA.get_bond_length_std(
    spatial_network_3)
bond_length_std_4, bond_length_vec_4 = NA.get_bond_length_std(
    spatial_network_4)

bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(
    spatial_network)
bond_angle_std_1, bond_angle_vec_1 = NA.get_bond_angle_std(
    spatial_network_1)
bond_angle_std_2, bond_angle_vec_2 = NA.get_bond_angle_std(
    spatial_network_2)
bond_angle_std_3, bond_angle_vec_3 = NA.get_bond_angle_std(
    spatial_network_3)
bond_angle_std_4, bond_angle_vec_4 = NA.get_bond_angle_std(
    spatial_network_4)

# print all standard deviations
println("bond length std: "*string(bond_length_std))
println("bond length std 1: "*string(bond_length_std_1))
println("bond length std 2: "*string(bond_length_std_2))
println("bond length std 3: "*string(bond_length_std_3))
println("bond length std 4: "*string(bond_length_std_4))
println("--------------------------------------------------")
println("bond angle std: "*string(bond_angle_std))
println("bond angle std 1: "*string(bond_angle_std_1))
println("bond angle std 2: "*string(bond_angle_std_2))
println("bond angle std 3: "*string(bond_angle_std_3))
println("bond angle std 4: "*string(bond_angle_std_4))

NG.plot_spatial_network(spatial_network)