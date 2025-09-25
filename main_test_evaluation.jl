
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Measurements
import Statistics

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\targeted\shell_nr_4\\"

all_order_metrics_dict = GU.load_h5_dict(analysis_data_path*"all_order_metrics.h5")

# extract the network type, that is srs from a string like "C:\\Users\\HemmannF\\OneDrive - Université de Fribourg\\structure_analysis\\analysis_data\\neural_network_targeted\\test_networks\\run_1\\srs_beta_0.2500_t_max_0.1709_t_gradient_0.1367_order_metrics.h5"

filenames_vec = all_order_metrics_dict["filenames_vec"]
network_type_vec = [match(r"/([^/]+)_nr_vertices", path).captures[1] for path in filenames_vec]
nr_vertices_vec = [parse(Int, match(r"nr_vertices_(\d+)", path).captures[1]) for path in filenames_vec]
bond_bending_vec = all_order_metrics_dict["bond_bending_const_vec"]
t_max_vec = all_order_metrics_dict["t_max_vec"]
t_gradient_vec = all_order_metrics_dict["t_gradient_vec"]

core_filenames_vec = [replace(split(path, '/') |> last, "_order_metrics.h5" => "") for path in filenames_vec] 
unique_core_filenames_vec = unique(core_filenames_vec)
#println(unique_core_filenames_vec)

# remove all unique core filenames that contain the string "lcs"
#unique_core_filenames_vec = filter(name -> !occursin("lcs", name), unique_core_filenames_vec)

order_metrics_vec = [
        "network_type",
        "nr_vertices_vec",
        "bond_bending_const_vec",
        "t_max_vec",
        "t_gradient_vec",
        "bond_length_std_vec",
        "bond_angle_std_vec",
        "dihedral_angle_entropy_vec",
        "bond_orientation_entropy_vec",
        "coordination_nr_mean_vec",
        "coordination_nr_std_vec",
        "vertex_homogeneity_metric_vec",
        "ring_size_mean_vec",
        "ring_size_std_vec",
        "ring_radius_mean_vec",
        "ring_radius_std_vec",
        "critical_pore_radius_vec",
        "anisotropy_metric_from_structure_factor_vec",
        "anisotropy_metric_from_structure_factor_bonds_vec",
        "hyperuniformity_alpha_vec_values",
        "hyperuniformity_alpha_vec_uncertainties",]

all_order_metrics_dict["network_type"] = network_type_vec   
all_order_metrics_dict["nr_vertices_vec"] = nr_vertices_vec
all_order_metrics_dict["hyperuniformity_alpha_vec_values"] = Measurements.value.(all_order_metrics_dict["hyperuniformity_alpha_vec"])
all_order_metrics_dict["hyperuniformity_alpha_vec_uncertainties"] = Measurements.uncertainty.(all_order_metrics_dict["hyperuniformity_alpha_vec"])

means_filename = "measured_order_metric_means.txt"
stds_filename = "measured_order_metric_stds.txt"

open(analysis_data_path*means_filename, "w") do io
    println(io, join(order_metrics_vec, '\t'))
end

open(analysis_data_path*stds_filename, "w") do io
    println(io, join(order_metrics_vec, '\t'))
end

# each core filename appears 20 times and has different values for the order
# metrics. Create a table that, for each filename contains the mean values of 
# all order metrics
for unique_core_filename in unique_core_filenames_vec
    println("Processing core filename: ", unique_core_filename)
    # get the mask of the current core filename
    mask = core_filenames_vec .== unique_core_filename
    means = [all_order_metrics_dict["network_type"][mask][1]]
    append!(means, [string(all_order_metrics_dict["nr_vertices_vec"][mask][1])])
    append!(means, [string(Statistics.mean(all_order_metrics_dict[metric][mask])) for metric in order_metrics_vec[3:end]])

    open(analysis_data_path*means_filename, "a") do io
        println(io, join(means, '\t'))
    end

    stds = [all_order_metrics_dict["network_type"][mask][1],
    string(all_order_metrics_dict["nr_vertices_vec"][mask][1]),
    all_order_metrics_dict["bond_bending_const_vec"][mask][1],
    all_order_metrics_dict["t_max_vec"][mask][1],
    all_order_metrics_dict["t_gradient_vec"][mask][1]]
    append!(stds, [string(Statistics.std(all_order_metrics_dict[metric][mask])) for metric in order_metrics_vec[6:end]])
    
    open(analysis_data_path*stds_filename, "a") do io
        println(io, join(stds, '\t'))
    end
end
