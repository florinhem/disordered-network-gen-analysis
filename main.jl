
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .BinaryDataAnalysis as BDA
import .GeneralUtilities as GU

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007

# julia --threads 20


dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_4\\"
filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_dict = BDA.get_binary_data_from_spatial_network(graph_dict;
    bond_radius = 0.35,
    filename = filename,
    save_result=true)

save_path = raw"..\analysis_data\random_networks\\"*filename

# load dict
complete_autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_complete.h5")

spectral_density_dict = BDA.get_spectral_density_by_wavevector_array_fft(structure_dict;
    save_complete_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = save_path,
    complete_autocovariance_fct_direction_dict = complete_autocovariance_fct_direction_dict)
#     
# complete_autocovariance_fct_direction_dict = BDA.get_complete_autocovariance_fct_by_sampling_vec_array(structure_dict;
#     nr_measurements_per_direction = 1000,
#     save_result = true,
#     save_path = raw"..\analysis_data\sample_name\random_networks\\"*filename)