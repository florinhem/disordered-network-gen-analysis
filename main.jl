
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .BinaryDataAnalysis as BDA
import .GeneralUtilities as GU

import Plots
import LinearAlgebra

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

# julia --threads 23

data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

filename = "216_vertices_T_0.1_heat_cool_0.05_per_mc_quenched"

data_dict = GU.load_h5_dict(data_path*filename*"_autocovariance_fct_direction.h5")

supercell_edge_length = LinearAlgebra.norm(autocovariance_fct_direction_dict["sampling_distance_array"][1,1,1,:] .- 
        autocovariance_fct_direction_dict["sampling_distance_array"][1,1,end,:])