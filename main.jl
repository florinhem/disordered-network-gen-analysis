
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

# julia --threads 23

function get_pore_size_distributions()
    file_count = 0

    structure_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_"

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285\run_"

    for i in 1:5

        current_structure_path = structure_path * string(i) * "\\"
        current_save_path = save_path * string(i) * "\\"

        # read all files in the current structure path
        structure_files = readdir(current_structure_path)

        for file in structure_files
            # load structure
            structure_dict = GU.load_h5_dict(current_structure_path * file)

            filename = file[1:end-13]

            # get pore size distribution
            pore_size_distribution_dict = NA.get_pore_size_distribution(structure_dict;
            save_result = true,
            save_path = current_save_path * filename)

            # get pore size distribution second moment
            pore_size_distribution_second_moment = NA.get_pore_size_distribution_second_moment(pore_size_distribution_dict)

            # load small scale order metrics dict 
            small_scale_order_metrics_dict = GU.load_h5_dict(current_save_path * filename * "_small_scale_order_metrics.h5")

            # save small scale order metrics dict with pore size distribution second moment
            small_scale_order_metrics_dict["pore_size_distribution_second_moment"] = pore_size_distribution_second_moment

            GU.save_dict_to_h5(small_scale_order_metrics_dict,
            current_save_path * filename*"_small_scale_order_metrics.h5")

            file_count += 1
            println("File ", file_count, " done.")
            
        end
        
    end

    return

end

get_pore_size_distributions()