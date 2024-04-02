
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import MetaGraphsNext
import SphericalHarmonics
import Formatting as Fmt

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007

dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

filename = "216_vertices_T_0.1_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)

NG.save_mesh_from_network(graph_dict, filename*"_thick_bonds"; save_path = dict_path, bond_radius = 0.3131)

filename = "216_vertices_T_0.4_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)

NG.save_mesh_from_network(graph_dict, filename*"_thick_bonds"; save_path = dict_path, bond_radius = 0.3131)