
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU



# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

analysis_data_path = raw"../analysis_data/neural_network_networks/srs/run_1//"
#analysis_data_path = raw"../analysis_data/neural_network_networks/ctn/run_1//"

all_filenames = []
all_roots = []
for (root, dirs, fs) in walkdir(analysis_data_path)
    for file in fs
        push!(all_filenames, file)
        push!(all_roots, root*"//")
    end
end

order_metrics_filenames = all_filenames[(occursin.("order_metrics", all_filenames)
        .&& .! occursin.("all_order_metrics", all_filenames))]

roots = all_roots[(occursin.("order_metrics", all_filenames)
        .&& .! occursin.("all_order_metrics", all_filenames))]

for i in eachindex(roots)
    order_metric_dict = GU.load_h5_dict(roots[i]*order_metrics_filenames[i])

    structure_factor_filename = replace(order_metrics_filenames[i], "order_metrics" => "structure_factor_array")
    structure_factor_dict = GU.load_h5_dict(roots[i]*structure_factor_filename)

    filename = replace(order_metrics_filenames[i], "_order_metrics.h5" => "")

    println("Processing file: ", filename)

    structure_factor_angle_averaged_dict = (
        NA.get_structure_factor_angle_averaged(
            structure_factor_dict;
            consider_bonds = true,
            gaussian_filter = true,
            gaussian_filter_sigma_x = 2*pi/25, 
            gaussian_filter_filtered_data_x_step_length = 2*pi/25,
            save_result = true,
            save_path = roots[i]*filename,
            label = nothing))

    spatial_network_path = replace(roots[i], "analysis_data" => "structures")

    order_metrics_dict = NA.get_order_metrics(
        filename,
    spatial_network_path,
    roots[i];
    l_max_steinhardt_q_l = 12,
    save_result = true,
    )
end



#network = NG.load_spatial_network_from_gml(
#    joinpath(network_path, filename * ".gml")
#)
#NG.plot_spatial_network(network)
