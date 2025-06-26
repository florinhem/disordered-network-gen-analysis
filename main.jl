
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import GLMakie
#import Statistics
#import Measurements
import Random

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

order_metric_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_networks\srs\\"
spatial_networks_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\neural_network_networks\srs\run_1\\"
networks_for_simulation_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\neural_network_networks\for_simulation\srs\run_1\\"

order_metrics_dict = GU.load_h5_dict(joinpath(order_metric_path, "all_order_metrics.h5"))

critical_pore_radius_vec = order_metrics_dict["critical_pore_radius_vec"]

# get all indices where the critical_pore_radius is between 0.386 and 0.45
indices = findall(x -> x >= 0.386 && x <= 0.45, critical_pore_radius_vec)

# go through all keys of the dict and get the values for the indices
order_metrics_dict_filtered = Dict{String, Any}()
for (key, value) in order_metrics_dict
    #println("Processing key: ", key)
    if length(size(order_metrics_dict[key])) == 1
        order_metrics_dict_filtered[key] = value[indices]
    else
        order_metrics_dict_filtered[key] = value[:, indices]
    end
end

println(length(order_metrics_dict_filtered["critical_pore_radius_vec"]))


# cut the filenames in the filenames_vec to only contain the last part after the last "\\"
# and remove the "_order_metrics.h5" extension
filenames_vec = order_metrics_dict_filtered["filenames_vec"]
filenames_vec = [replace(basename(filename), "_order_metrics.h5" => "") for filename in filenames_vec]

# shuffle the vector
shuffled_filenames_vec = Random.shuffle(filenames_vec)

# save the shuffled filenames to a txt file
open(networks_for_simulation_path*"shuffled_filenames_vec.txt", "w") do io
    for line in shuffled_filenames_vec
        println(io, line)
    end
end

# for all filenames, load the spatial network and save the network for simulation

for filename in shuffled_filenames_vec
    # load the spatial network
    spatial_network = NG.load_spatial_network_from_gml(joinpath(spatial_networks_path, filename * ".gml"))
    # save the network for simulation
    spatial_network = NG.get_spatial_network_for_simulation!(
    spatial_network;
    save_result = true,
    filename = filename,
    save_path= networks_for_simulation_path,)
end
