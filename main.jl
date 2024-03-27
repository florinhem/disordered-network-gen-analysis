
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import MetaGraphsNext
import SphericalHarmonics

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007

evolution_dict = NA.get_evolution_dict(nr_vertices = 64, network_type = "diamond")
graph_dict_diamond = NG.get_periodic_network(evolution_dict)

evolution_dict = NA.get_evolution_dict(nr_vertices = 64, network_type = "simple cubic")
graph_dict_sc = NG.get_periodic_network(evolution_dict)

central_vertex = 10

l_max = 8

single_vertex_q_l_sc = NA.get_q_l_averaged_single_vertex_dict(graph_dict_sc,
central_vertex,
l_max)

single_vertex_q_l_diamond = NA.get_q_l_averaged_single_vertex_dict(graph_dict_diamond,
central_vertex,
l_max)
