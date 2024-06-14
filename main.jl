
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
# which is the cube root of the number of vertices times 2/sqrt(3)

# julia --threads 23

# filter out all filenames that contain "heated_for_0.5_steps" and
# filter out the same entries from all other vectors

for i in [1,3,4,5]

    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_"*string(i)*"\\"

    order_metrics_dict = NA.get_small_length_scale_order_metrics_all_files(analysis_data_path;
        l_max_steinhardt_q_l = 12,
        save_result = true,)

    println("run $(i) done")

end