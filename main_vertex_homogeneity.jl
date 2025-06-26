
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

function my_func_1()
    

    analysis_data_path = raw"../analysis_data/neural_network_networks/lcs/"
    network_path = replace(analysis_data_path, "analysis_data" => "structures")
    
    all_filenames = []
    all_roots = []
    for (root, dirs, fs) in walkdir(analysis_data_path)
        for file in fs
            push!(all_filenames, file)
            push!(all_roots, root*"//")
        end
    end
    
    order_metrics_filenames = all_filenames[(occursin.("order_metrics", all_filenames)
            .&& .! occursin.("all_order_metrics", all_filenames)
            .&& endswith.(all_filenames, ".h5"))]
    
    roots = all_roots[(occursin.("order_metrics", all_filenames)
            .&& .! occursin.("all_order_metrics", all_filenames)
            .&& endswith.(all_filenames, ".h5"))]
    
    # print progress in percentage
    count = 0
    total = length(roots)
    println("Processing ", total, " files...")
    
    for i in eachindex(roots)
        order_metric_dict = GU.load_h5_dict(roots[i]*order_metrics_filenames[i])
    
        network_filename = replace(order_metrics_filenames[i], "_order_metrics.h5" => ".gml")
        network_root = replace(roots[i], "analysis_data" => "structures")
    
        spatial_network = NG.load_spatial_network_from_gml(network_root*network_filename)
    
        filename = replace(network_filename, ".gml" => "")
    
        correlation_functions_dict = NA.get_correlation_functions(
        spatial_network;
        save_result = true,
        save_path = roots[i]*filename,
        label = nothing)
    
        vertex_homogeneity_metric = NA.get_vertex_homogeneity_metric(correlation_functions_dict)
    
        order_metric_dict["vertex_homogeneity_metric_vec"] = vertex_homogeneity_metric
    
        GU.save_dict_to_h5(order_metric_dict, roots[i]*order_metrics_filenames[i])
    
        #if i == 1
        #    println(order_metric_dict)
        #end
    
        count += 1
    
        if count % 10 == 0
            println("Processed ", count, " files (", round(count / total * 100, digits=2), "%)")
        end
    end
end 

# Call the function to execute the analysis
my_func_1()


function my_func_2()
    

    analysis_data_path = raw"../analysis_data/neural_network_networks/srs/"
    network_path = replace(analysis_data_path, "analysis_data" => "structures")
    
    all_filenames = []
    all_roots = []
    for (root, dirs, fs) in walkdir(analysis_data_path)
        for file in fs
            push!(all_filenames, file)
            push!(all_roots, root*"//")
        end
    end
    
    order_metrics_filenames = all_filenames[(occursin.("order_metrics", all_filenames)
            .&& .! occursin.("all_order_metrics", all_filenames)
            .&& endswith.(all_filenames, ".h5"))]
    
    roots = all_roots[(occursin.("order_metrics", all_filenames)
            .&& .! occursin.("all_order_metrics", all_filenames)
            .&& endswith.(all_filenames, ".h5"))]
    
    # print progress in percentage
    count = 0
    total = length(roots)
    println("Processing ", total, " files...")
    
    for i in eachindex(roots)
        order_metric_dict = GU.load_h5_dict(roots[i]*order_metrics_filenames[i])
    
        network_filename = replace(order_metrics_filenames[i], "_order_metrics.h5" => ".gml")
        network_root = replace(roots[i], "analysis_data" => "structures")
    
        spatial_network = NG.load_spatial_network_from_gml(network_root*network_filename)
    
        filename = replace(network_filename, ".gml" => "")
    
        correlation_functions_dict = NA.get_correlation_functions(
        spatial_network;
        save_result = true,
        save_path = roots[i]*filename,
        label = nothing)
    
        vertex_homogeneity_metric = NA.get_vertex_homogeneity_metric(correlation_functions_dict)
    
        order_metric_dict["vertex_homogeneity_metric_vec"] = vertex_homogeneity_metric
    
        GU.save_dict_to_h5(order_metric_dict, roots[i]*order_metrics_filenames[i])
    
        #if i == 1
        #    println(order_metric_dict)
        #end
    
        count += 1
    
        if count % 10 == 0
            println("Processed ", count, " files (", round(count / total * 100, digits=2), "%)")
        end
    end
end 

# Call the function to execute the analysis
my_func_2()