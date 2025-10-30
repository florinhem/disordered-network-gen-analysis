"""
these functions provide utilities for analyzing networks
"""

"""
Convert cartesian to spherical coordinates
"""
function convert_cartesian_to_spherical(cartesian_vec::Vector)

    # check if vector is 3d
    if length(cartesian_vec) !== 3
        @error "conversion to spherical coordinates only implemented in 3d"
        return []
    end

    # calculate r, theta and phi 
    r_length = LinearAlgebra.norm(cartesian_vec)
    theta = atan( sqrt(cartesian_vec[1]^2 + cartesian_vec[2]^2), 
        cartesian_vec[3] )
    phi = atan(cartesian_vec[2], cartesian_vec[1] )

    return [r_length, theta, phi]
end


"""
Average y(x) over different network realizations. The y_arr contains values for
the same network and different x along the 1 dimension and for the same x and 
different networks along the 2 dimension
"""
function get_multiple_network_average_vec(y_arr::Array{Real})

    # get the number of x values, for which y(x) was evaluated
    nr_x_values = length(y_arr[:,1])

    # initialize vector for average values of y
    y_average_vec = Vector{Measurements.Measurement{Float64}}(undef, 
        nr_x_values)

    # loop through index corresponding to different x values and average over
    # different realizations
    for x_index in 1:nr_x_values

        y_average_vec[x_index] = Measurements.measurement(
            Statistics.mean(y_arr[x_index,:]),
            Statistics.std(y_arr[x_index,:]))
    end

    return y_average_vec
end


"""
Average the quantity y over different network realizations. The y_vec contains
this quantity for different networks
"""
function get_multiple_network_average(y_vec::Vector{Real})

    # calculate mean and standard deviation of quantity
    y_average = Measurements.measurement(
            Statistics.mean(y_vec),
            Statistics.std(y_vec)
        )

    return y_average
end


"""
Get the minimal and maximal vertex coordinates along all three directions
"""
function get_min_max_vertex_coords(spatial_network)
    # get the minimal and maximal vertex positions along the three axes
    min_vertex_coords = [Inf, Inf, Inf]
    max_vertex_coords = [-Inf, -Inf, -Inf]
    for vertex in MetaGraphsNext.labels(spatial_network)
        for i in 1:3
            if spatial_network[vertex]["position"][i] < min_vertex_coords[i]
                min_vertex_coords[i] = spatial_network[vertex]["position"][i]
            end
            if spatial_network[vertex]["position"][i] > max_vertex_coords[i]
                max_vertex_coords[i] = spatial_network[vertex]["position"][i]
            end
        end
    end
    return min_vertex_coords, max_vertex_coords
end


"""
Convert a dictionary of Steinhardt local bond order parameters q_l with
0 <= l <= l_max into a vector of length l_max+1
"""
function convert_q_l_dict_to_vec(q_l_dict::Dict, l_max::Int)

    # initialize vector for q_l values with same data type as dictionary
    # which could be Measurements.Measurement{Float64} or Float64
    q_l_vec = Vector{typeof(q_l_dict[1])}(undef, l_max+1)

    for l in 0:l_max
        q_l_vec[l+1] = q_l_dict[l]
    end

    return q_l_vec
end


"""
Calculate the average nr of monte carlo steps for quenching for all evolution
dicts. I got the result 13.7 +- 6.74 for 216 vertices
"""
function get_monte_carlo_steps_for_quenching_vec(
    evolution_dicts_directory_path::String; 
    nr_runs = 5)

    nr_monte_carlo_steps_for_quenching_vec = Vector{Float64}()

    total_nr_dicts = 5*136
    counter = 0

    for i in 1:nr_runs
        current_path = evolution_dicts_directory_path*"run_"*string(i)*"\\"

        # get all files in directory
        filenames = readdir(current_path)

        # get filenames of all evolution dicts
        filenames_evolution_dicts = filter(filename 
            -> endswith(filename, "_evolution.h5"), filenames)

        for filename in filenames_evolution_dicts

            println("Progress: "*string(counter/total_nr_dicts*100)*" %")
            counter += 1

            # load evolution dict
            evolution_dict = GU.load_h5_dict(current_path*filename)

            evolution_dict["nr_monte_carlo_steps_per_temperature_vec"]

            nr_monte_carlo_moves_per_step = 18*evolution_dict["nr_vertices"]

            nr_monte_carlo_moves_before_quenching = (
                nr_monte_carlo_moves_per_step*sum(evolution_dict[
                    "nr_monte_carlo_steps_per_temperature_vec"][1:end-1]))

            nr_monte_carlo_moves_for_quenching = (length(
                    evolution_dict["move_accepted_vec"])
                -nr_monte_carlo_moves_before_quenching)

            nr_monte_carlo_steps_for_quenching = (
                nr_monte_carlo_moves_for_quenching
                /nr_monte_carlo_moves_per_step)

            if nr_monte_carlo_steps_for_quenching > 49
                println(i, filename)
            end

            push!(nr_monte_carlo_steps_for_quenching_vec,
                nr_monte_carlo_steps_for_quenching)
        end

    end

    sort!(nr_monte_carlo_steps_for_quenching_vec)

    return nr_monte_carlo_steps_for_quenching_vec
end


"""
For a single file, calculate several functions that characterize the structure
"""
function get_all_dicts_from_network_single_file(
    filename::String,
    spatial_network_path::String,
    analysis_data_path::String;
    digital_sphere_mask_path 
        = raw"..\analysis_data\random_networks\digital_sphere_masks\\",
    pore_size_sampling_grid_size = 0.2,
    max_pore_radius = 3.0,
    hyperuniformity_min_wavenumber_to_consider::Float64 = 0.0,
    periodic_boundary_conditions::Bool = true,
    print_progress::Bool = false,
    print_lock = Threads.ReentrantLock())

    # get network
    spatial_network = NG.load_spatial_network_from_gml(
        spatial_network_path*filename*".gml")
    
    # get ring size distribution
    ring_size_distribution_dict = get_ring_size_distribution(
        spatial_network;
        periodic_boundary_conditions = periodic_boundary_conditions,
        save_result = true,
        save_path = analysis_data_path*filename)

    # get ring radius distribution

    ring_radius_distribution_dict = get_ring_radius_distribution(
        spatial_network,
        ring_size_distribution_dict;
        save_result = true,
        save_path = analysis_data_path*filename)

    # get structure factor by wavevector array for vertices
    structure_factor_dict = get_structure_factor_by_wavevector_array(
        spatial_network;
        consider_bonds = false,
        periodic_boundary_conditions = periodic_boundary_conditions,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing,
        print_progress = print_progress,
        thread_nr = Threads.threadid(),
        print_lock = print_lock)
    
    # get structure factor by wavevector array for bonds
    structure_factor_bonds_dict = get_structure_factor_by_wavevector_array(
        spatial_network;
        consider_bonds = true,
        periodic_boundary_conditions = periodic_boundary_conditions,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing,
        print_progress = print_progress,
        thread_nr = Threads.threadid(),
        print_lock = print_lock)

    # get angle averaged structure factor for vertices
    structure_factor_angle_averaged_dict = get_structure_factor_angle_averaged(
        structure_factor_dict;
        consider_bonds = false,
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing)

    # get angle averaged structure factor for bonds
    structure_factor_bonds_angle_averaged_dict = (
        get_structure_factor_angle_averaged(
            structure_factor_bonds_dict;
            consider_bonds = true,
            gaussian_filter = true,
            gaussian_filter_sigma_x = 2*pi/25, 
            gaussian_filter_filtered_data_x_step_length = 2*pi/25,
            save_result = true,
            save_path = analysis_data_path*filename,
            label = nothing))

    # get correlation functions
    correlation_functions_dict = get_correlation_functions(
        spatial_network;
        distance_histogram_bin_width = 0.02,
        periodic_boundary_conditions = periodic_boundary_conditions,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing)
    
    # get pore size distribution
    pore_size_distribution_dict = get_pore_size_distribution(
        spatial_network;
        sampling_grid_size = pore_size_sampling_grid_size,
        max_pore_radius = max_pore_radius,
        periodic_boundary_conditions = periodic_boundary_conditions,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing,
        digital_sphere_mask_path = digital_sphere_mask_path,
        print_progress = print_progress,
        thread_nr = Threads.threadid(),
        print_lock = print_lock)
    
    # get all order metrics for the network
    order_metrics_dict = get_order_metrics(
        filename,
        spatial_network_path,
        analysis_data_path;
        l_max_steinhardt_q_l = 12,
        hyperuniformity_min_wavenumber_to_consider
            = hyperuniformity_min_wavenumber_to_consider,
        save_result = true,
        )

    return
end


"""
For each structure dict in a list of run and filenames, calculate
several functions that characterize the structure
"""
function get_all_dicts_from_networks_single_thread(run_and_filename_chunk,
    spatial_networks_path::String,
    analysis_data_path::String;
    digital_sphere_mask_path 
        = raw"..\analysis_data\random_networks\digital_sphere_masks\\",
    pore_size_sampling_grid_size = 0.2,
    max_pore_radius = 3.0,
    periodic_boundary_conditions::Bool = true,
    print_progress::Bool = false,
    print_lock = Threads.ReentrantLock())

    # loop through files
    for run_and_filename in run_and_filename_chunk

        # print current thread id and run/filename if desired
        if print_progress
            lock(print_lock) do
                Format.printfmtln("Thread {1} is handling file {2}",
                    Threads.threadid(), run_and_filename)
            end
        end

        get_all_dicts_from_network_single_file(
            run_and_filename[7:end],
            spatial_networks_path*run_and_filename[1:6],
            analysis_data_path*run_and_filename[1:6];
            digital_sphere_mask_path = digital_sphere_mask_path,
            pore_size_sampling_grid_size = pore_size_sampling_grid_size,
            max_pore_radius = max_pore_radius,
            periodic_boundary_conditions = periodic_boundary_conditions,
            print_progress = print_progress,
            print_lock = print_lock)
        
    end

    return
end


"""
For all structure dicts in a folder, calculate several functions that
characterize the structures using multithreading
"""
function get_all_dicts_from_networks_multithreading(
    spatial_networks_path::String,
    analysis_data_path::String;
    digital_sphere_mask_path 
        = raw"..\analysis_data\random_networks\digital_sphere_masks\\",
    pore_size_sampling_grid_size = 0.2,
    max_pore_radius = 3.0,
    periodic_boundary_conditions::Bool = true,
    print_progress::Bool = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())

    dicts_to_save = ["_correlation_functions.h5",
        "_pore_size_distribution.h5",
        "_ring_radius_distribution.h5",
        "_ring_size_distribution.h5",
        "_structure_factor_angle_averaged.h5",
        "_structure_factor_bonds_angle_averaged.h5",
        "_structure_factor_array.h5",
        "_structure_factor_bonds_array.h5",
        "_order_metrics.h5"]

    # vector for run and filename
    run_and_filename_vec = []

    # loop through runs
    for i in runs_vec

        spatial_networks_path_current_run = (
            spatial_networks_path*"run_"*string(i)*"/")

        filenames = readdir(spatial_networks_path_current_run)
        filenames_evolution_dicts = filter(filename 
            -> endswith(filename, "_evolution.h5"), filenames)
        filenames_spatial_networks = [filename[1:end-13] 
            for filename in filenames_evolution_dicts]

        # filter out filenames that are already fully analyzed
        filenames_spatial_networks_not_analyzed = filter(filename 
            -> !all(dict -> isfile(analysis_data_path*"run_"*string(i)*"/"*
                filename*dict), dicts_to_save), filenames_spatial_networks)

        # print number of files not analyzed
        if print_progress
            lock(print_lock) do
                println(length(filenames_spatial_networks_not_analyzed), 
                " files not analyzed in run ", i)
            end
        end

        # append to list of all structure dict paths
        append!(run_and_filename_vec, "run_"*string(i)*"/" 
            .* filenames_spatial_networks_not_analyzed)

    end

    # split filenames into chunks for multi-threading
    run_and_filename_chunks = Iterators.partition(run_and_filename_vec, 
        length(run_and_filename_vec) ÷ Threads.nthreads())

    # run all filename chunk in parallel in different threads
    map(run_and_filename_chunks) do run_and_filename_chunk

        Threads.@spawn get_all_dicts_from_networks_single_thread(
            run_and_filename_chunk,
            spatial_networks_path,
            analysis_data_path;
            digital_sphere_mask_path = digital_sphere_mask_path,
            pore_size_sampling_grid_size = pore_size_sampling_grid_size,
            max_pore_radius = max_pore_radius,
            periodic_boundary_conditions = periodic_boundary_conditions,
            print_progress = print_progress,
            print_lock = print_lock)
    end

    return
end


"""
Extract the parameters bond_bending_const, t_max and t_gradient of the Monte 
Carlo algorithm from a filename
"""
function extract_parameters(filename::String)
    pattern = r"beta_([0-9.]+)_t_max_([0-9.]+)_t_gradient_([0-9.]+)"
    m = match(pattern, filename)
    if m !== nothing
        bond_bending_const = parse(Float64, m.captures[1])
        t_max = parse(Float64, m.captures[2])
        t_gradient = parse(Float64, m.captures[3])
        return bond_bending_const, t_max, t_gradient
    else
        error("Pattern not found in filename: $filename")
    end
end


"""
For given data paths and filename, calculate all order metrics
"""
function get_order_metrics(filename::String,
    network_path::String,
    analysis_data_path::String;
    l_max_steinhardt_q_l::Int64 = 12,
    hyperuniformity_min_wavenumber_to_consider::Float64 = 0.0,
    save_result = false)

    # load network
    spatial_network = NG.load_spatial_network_from_gml(
        network_path*filename*".gml")

    # if the key "total_energy" is present in the spatial network,
    # save its
    if haskey(spatial_network[], "total_energy")
        total_keating_energy = spatial_network[]["total_energy"]
    end

    # Measure the standard deviation of bond lengths
    bond_length_std, bond_length_vec = get_bond_length_std(spatial_network)

    # Measure the standard deviation of bond angles
    bond_angle_std, bond_angle_vec = get_bond_angle_std(spatial_network)

    # Measure the entropy of dihedral angles
    dihedral_angle_entropy = get_dihedral_angle_entropy(
        spatial_network)

    # Mesure the entropy of bond orientations
    bond_orientation_entropy = get_bond_orientation_entropy(spatial_network)

    # Get the coordination number statistics
    coordination_nr_mean, coordination_nr_std = get_coordination_nr_statistics(
        spatial_network)

    # get Steinhardt local bond order parameters and store them in a vector
    q_l_total_network_mean_dict = get_q_l_total_network_mean_dict(
        spatial_network, l_max_steinhardt_q_l)

    q_l_vec = convert_q_l_dict_to_vec(q_l_total_network_mean_dict, 
        l_max_steinhardt_q_l)

    # TODO: check if q_l_uncertainties are properly saved, because neural 
    # networks predicts them to be the same as q_l_values

    # load ring size distribution
    ring_size_distribution_dict = GU.load_h5_dict(
        analysis_data_path*filename*"_ring_size_distribution.h5")

    # get ring size statistics
    ring_size_mean, ring_size_std = get_ring_statistics(
        ring_size_distribution_dict)

    # load ring radius distribution
    ring_radius_distribution_dict = GU.load_h5_dict(
        analysis_data_path*filename*"_ring_radius_distribution.h5")

    # get ring radius statistics
    ring_radius_mean, ring_radius_std = get_ring_radius_statistics(
        ring_radius_distribution_dict)

    # load correlation functions
    correlation_functions_dict = GU.load_h5_dict(
        analysis_data_path*filename*"_correlation_functions.h5")

    # get vertex homogeneity metric
    vertex_homogeneity_metric = get_vertex_homogeneity_metric(
        correlation_functions_dict, consider_uncoordinated_vertices=false)

    # get uncoordinated neighbor distance
    uncoordinated_neighbor_distance = get_vertex_homogeneity_metric(
        correlation_functions_dict, consider_uncoordinated_vertices=true)

    # load pore size distribution
    pore_size_distribution_dict = GU.load_h5_dict(
        analysis_data_path*filename*"_pore_size_distribution.h5")

    # get critical pore radius
    critical_pore_radius = (
        get_critical_pore_radius(pore_size_distribution_dict))

    # load angle averaged structure factor for vertices
    structure_factor_angle_averaged_dict = GU.load_h5_dict(
        analysis_data_path*filename*"_structure_factor_angle_averaged.h5")

    # get anisotropy metric from structure factor for vertices
    anisotropy_metric_from_structure_factor = (
        get_anisotropy_metric_from_structure_factor(
            structure_factor_angle_averaged_dict))

    # load angle averaged structure factor for bonds
    structure_factor_bonds_angle_averaged_dict = GU.load_h5_dict(
        analysis_data_path*filename
        *"_structure_factor_bonds_angle_averaged.h5")

    # get anisotropy metric from structure factor for bonds
    anisotropy_metric_from_structure_factor_bonds = (
        get_anisotropy_metric_from_structure_factor(
            structure_factor_bonds_angle_averaged_dict))

    # get the alpha value that captures whether the network is hyperuniform 
    hyperuniformity_alpha = (
        get_hyperuniformity_alpha(structure_factor_bonds_angle_averaged_dict;
            min_wavenumber_to_consider =
                hyperuniformity_min_wavenumber_to_consider))

    # create dict to save
    order_metrics_dict = Dict(
        "bond_length_std" => bond_length_std,
        "bond_angle_std" => bond_angle_std,
        "dihedral_angle_entropy" => dihedral_angle_entropy,
        "bond_orientation_entropy" => bond_orientation_entropy,
        "coordination_nr_mean" => coordination_nr_mean,
        "coordination_nr_std" => coordination_nr_std,
        "q_l_vec" => q_l_vec,
        "ring_size_mean" => ring_size_mean,
        "ring_size_std" => ring_size_std,
        "ring_radius_mean" => ring_radius_mean,
        "ring_radius_std" => ring_radius_std,
        "vertex_homogeneity_metric" => vertex_homogeneity_metric,
        "uncoordinated_neighbor_distance" => uncoordinated_neighbor_distance,
        "critical_pore_radius" => critical_pore_radius,
        "anisotropy_metric_from_structure_factor" 
            => anisotropy_metric_from_structure_factor,
        "anisotropy_metric_from_structure_factor_bonds" 
            => anisotropy_metric_from_structure_factor_bonds,
        "hyperuniformity_alpha" => hyperuniformity_alpha,
    )

    if haskey(spatial_network[], "total_energy")
        order_metrics_dict["total_keating_energy"] = (total_keating_energy)
    end

    if save_result
        GU.save_dict_to_h5(order_metrics_dict, 
            analysis_data_path*filename*"_order_metrics.h5")
    end

    return order_metrics_dict
end


function extract_folder_and_filename(path::AbstractString)
    # Normalize and split the path into components (handles / and \ correctly)
    parts = splitpath(path)

    # Get last two components: the folder (e.g., "run_1") and filename
    folder = parts[end - 1]
    filename = parts[end]

    # Strip known suffix from filename
    name = replace(filename, "_order_metrics.h5" => "")

    return "$folder/$name"
end


function get_order_metrics_all_files(
    analysis_data_path::String;
    l_max_steinhardt_q_l::Int64 = 12,
    save_result::Bool = false,
    save_algorithm_parameters_from_filename::Bool = false)

    # get all dictionaries containing order metrics
    all_filenames = []
    for (root, dirs, fs) in walkdir(analysis_data_path)
        for file in fs
            push!(all_filenames, joinpath(root, file))
        end
    end
    # filter for order metrics files with type h5
    order_metrics_filenames = [filename for filename in all_filenames
        if (occursin("order_metrics", filename) 
            && !occursin("all_order_metrics", filename)
            && endswith(filename, ".h5"))]

    # initialize vectors for all order metrics
    total_keating_energy_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    bond_length_std_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    bond_angle_std_vec = Vector{Float64}(undef, 
        length(order_metrics_filenames))
    dihedral_angle_entropy_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    bond_orientation_entropy_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    coordination_nr_mean_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    coordination_nr_std_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    q_l_mat = Matrix{Measurements.Measurement{Float64}}(undef, 
        l_max_steinhardt_q_l+1, length(order_metrics_filenames))
    vertex_homogeneity_metric_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    uncoordinated_neighbor_distance_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    ring_size_mean_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    ring_size_std_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    ring_radius_mean_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    ring_radius_std_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    critical_pore_radius_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    anisotropy_metric_from_structure_factor_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    anisotropy_metric_from_structure_factor_bonds_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    hyperuniformity_alpha_vec = Vector{Measurements.Measurement{Float64}}(
        undef, length(order_metrics_filenames))

    # loop through order metric filenames
    for i in eachindex(order_metrics_filenames)

        # load order metrics
        order_metrics_dict = GU.load_h5_dict(
            order_metrics_filenames[i])

        # get all order metrics and save them to the corresponding vectors
        if haskey(order_metrics_dict, "total_keating_energy")
            total_keating_energy_vec[i] = (
                order_metrics_dict["total_keating_energy"])
        else
            total_keating_energy_vec[i] = 0.0
        end
        bond_length_std_vec[i] = order_metrics_dict["bond_length_std"]
        bond_angle_std_vec[i] = order_metrics_dict["bond_angle_std"]
        dihedral_angle_entropy_vec[i] = (
            order_metrics_dict["dihedral_angle_entropy"])
        bond_orientation_entropy_vec[i] = (
            order_metrics_dict["bond_orientation_entropy"])
        coordination_nr_mean_vec[i] = (
            order_metrics_dict["coordination_nr_mean"])
        coordination_nr_std_vec[i] = (
            order_metrics_dict["coordination_nr_std"])
        q_l_mat[:,i] = order_metrics_dict["q_l_vec"]
        vertex_homogeneity_metric_vec[i] = (
            order_metrics_dict["vertex_homogeneity_metric"])
        uncoordinated_neighbor_distance_vec[i] = (
            order_metrics_dict["uncoordinated_neighbor_distance"])
        ring_size_mean_vec[i] = order_metrics_dict["ring_size_mean"]
        ring_size_std_vec[i] = order_metrics_dict["ring_size_std"]
        ring_radius_mean_vec[i] = order_metrics_dict["ring_radius_mean"]
        ring_radius_std_vec[i] = order_metrics_dict["ring_radius_std"]
        critical_pore_radius_vec[i] = (
            order_metrics_dict["critical_pore_radius"])
        anisotropy_metric_from_structure_factor_vec[i] = (
            order_metrics_dict["anisotropy_metric_from_structure_factor"])
        anisotropy_metric_from_structure_factor_bonds_vec[i] = (
            order_metrics_dict[
                "anisotropy_metric_from_structure_factor_bonds"])
        hyperuniformity_alpha_vec[i] = (
            order_metrics_dict["hyperuniformity_alpha"])

    end

    # sort all vectors with respect to the keating energy
    order_metrics_filenames = order_metrics_filenames[
        sortperm(total_keating_energy_vec)]
    bond_length_std_vec = bond_length_std_vec[
        sortperm(total_keating_energy_vec)]
    bond_angle_std_vec = bond_angle_std_vec[
        sortperm(total_keating_energy_vec)]
    dihedral_angle_entropy_vec = dihedral_angle_entropy_vec[
        sortperm(total_keating_energy_vec)]
    bond_orientation_entropy_vec = bond_orientation_entropy_vec[
        sortperm(total_keating_energy_vec)]
    coordination_nr_mean_vec = coordination_nr_mean_vec[
        sortperm(total_keating_energy_vec)]
    coordination_nr_std_vec = coordination_nr_std_vec[
        sortperm(total_keating_energy_vec)]
    q_l_mat = q_l_mat[:, sortperm(total_keating_energy_vec)]
    vertex_homogeneity_metric_vec = (
        vertex_homogeneity_metric_vec[sortperm(total_keating_energy_vec)])
    uncoordinated_neighbor_distance_vec = (
        uncoordinated_neighbor_distance_vec[sortperm(total_keating_energy_vec)]
    )
    ring_size_mean_vec = ring_size_mean_vec[
        sortperm(total_keating_energy_vec)]
    ring_size_std_vec = ring_size_std_vec[
        sortperm(total_keating_energy_vec)]
    ring_radius_mean_vec = ring_radius_mean_vec[
        sortperm(total_keating_energy_vec)]
    ring_radius_std_vec = ring_radius_std_vec[
        sortperm(total_keating_energy_vec)]
    critical_pore_radius_vec = (critical_pore_radius_vec[
            sortperm(total_keating_energy_vec)])
    anisotropy_metric_from_structure_factor_vec = (
        anisotropy_metric_from_structure_factor_vec[
            sortperm(total_keating_energy_vec)])
    anisotropy_metric_from_structure_factor_bonds_vec = (
        anisotropy_metric_from_structure_factor_bonds_vec[
            sortperm(total_keating_energy_vec)])
    hyperuniformity_alpha_vec = hyperuniformity_alpha_vec[
        sortperm(total_keating_energy_vec)]

    sort!(total_keating_energy_vec)

    # cut the filenames vector to only contain the filenames and the parent 
    # directory
    order_metrics_filenames = [
        extract_folder_and_filename(filename) for filename in 
        order_metrics_filenames]

    # create dict to save
    order_metrics_dict = Dict(
        "total_keating_energy_vec" => total_keating_energy_vec,
        "bond_length_std_vec" => bond_length_std_vec,
        "bond_angle_std_vec" => bond_angle_std_vec,
        "dihedral_angle_entropy_vec" => dihedral_angle_entropy_vec,
        "bond_orientation_entropy_vec" => bond_orientation_entropy_vec,
        "coordination_nr_mean_vec" => coordination_nr_mean_vec,
        "coordination_nr_std_vec" => coordination_nr_std_vec,
        "q_l_mat" => q_l_mat,
        "vertex_homogeneity_metric_vec" => vertex_homogeneity_metric_vec,
        "uncoordinated_neighbor_distance_vec" 
            => uncoordinated_neighbor_distance_vec,
        "ring_size_mean_vec" => ring_size_mean_vec,
        "ring_size_std_vec" => ring_size_std_vec,
        "ring_radius_mean_vec" => ring_radius_mean_vec,
        "ring_radius_std_vec" => ring_radius_std_vec,
        "critical_pore_radius_vec" => critical_pore_radius_vec,
        "anisotropy_metric_from_structure_factor_vec" 
            => anisotropy_metric_from_structure_factor_vec,
        "anisotropy_metric_from_structure_factor_bonds_vec" 
            => anisotropy_metric_from_structure_factor_bonds_vec,
        "hyperuniformity_alpha_vec" => hyperuniformity_alpha_vec,
        "filenames_vec" => order_metrics_filenames
    )

    if save_algorithm_parameters_from_filename
        # get algorithm parameters from filename
        bond_bending_const_vec = Vector{Float64}(undef,
            length(order_metrics_filenames))
        t_max_vec = Vector{Float64}(undef,
            length(order_metrics_filenames))
        t_gradient_vec = Vector{Float64}(undef,
            length(order_metrics_filenames))
        for i in eachindex(order_metrics_filenames)
            bond_bending_const_vec[i], t_max_vec[i], t_gradient_vec[i] = (
                extract_parameters(order_metrics_filenames[i]))
        end

        # add algorithm parameters to dict
        order_metrics_dict["bond_bending_const_vec"] = (
            bond_bending_const_vec)
        order_metrics_dict["t_max_vec"] = t_max_vec
        order_metrics_dict["t_gradient_vec"] = t_gradient_vec
    end

    # save dict to h5 file if desired
    if save_result
        GU.save_dict_to_h5(order_metrics_dict,
        analysis_data_path*"all_order_metrics.h5")
    end

    return order_metrics_dict
end


"""
Save the dictionary containing all order metrics to a csv file
"""
function save_order_metrics_dict_to_csv(
    order_metrics_dict::Dict,
    save_path::String)
    
    # delete the key "q_mat"
    delete!(order_metrics_dict, "q_l_mat")

    # convert dict to DataFrame
    df = DataFrames.DataFrame(order_metrics_dict)

    # save to csv file
    CSV.write(
        save_path*"all_order_metrics.csv",
        df;
        writeheader=true,
        delim=',',
        quotechar='"',
        stringtype=String
    )

    return
end


"""
Print the melting temperatures for various network configurations.
"""
function print_melting_temperatures(
    ;
    network_type_vec,
    nr_vertices_vec,
    bond_bending_const_vec,
    theta_ground_state_vec,
    acceptance_probability_vec,
    relax_globally_after_threshold_cycle_vec,
    shell_nr_vec
    )

    println(Threads.nthreads())
    # Pair nr_vertices and network_type using zip
    vertex_network_pairs = zip(network_type_vec, nr_vertices_vec)

    # Create the iterator for all combinations
    Iter = collect(Iterators.product(
        vertex_network_pairs,
        bond_bending_const_vec,
        theta_ground_state_vec,
        acceptance_probability_vec,
        relax_globally_after_threshold_cycle_vec,
        shell_nr_vec
    ))
    data=[]
    data_lock = Threads.ReentrantLock()
    Threads.@threads for (
        (network_type, nr_vertices),
        bond_bending_const,
        theta_ground_state,
        acceptance_probability,
        relax_globally_after_threshold_cycle,
        shell_nr) in Iter
        println(
            "$network_type"*", "*
            "$nr_vertices"*", "*
            "$bond_bending_const"*", "*
            "$theta_ground_state"*", "*
            "$acceptance_probability"*", "*
            "$relax_globally_after_threshold_cycle"*", "*
            "$shell_nr" )
        evolution_dict = get_evolution_dict(;
            nr_vertices = nr_vertices,
            network_type=network_type,
            bond_bending_const=bond_bending_const,
            min_ring_size=3,
            theta_ground_state=theta_ground_state,
            relax_globally_after_threshold_cycle
                =relax_globally_after_threshold_cycle,
            shell_nr=shell_nr
        )
        spatial_network = NG.get_periodic_network(evolution_dict)
        spatial_network_copy=deepcopy(spatial_network)

        # E start
        E_start=NG.get_total_energy_keating(spatial_network_copy)
        T_melt_vec=[]

        #We choose nr_vertices so that on average we get once every vertex
        for i in collect(range(0, 0, length=nr_vertices)) 
            spatial_network=deepcopy(spatial_network_copy)
            # E end
            switched_chain = NG.get_random_chain(
                spatial_network;
                declined_chains = [],
                remaining_chains = [],
                min_ring_size = evolution_dict["min_ring_size"])
            # switch bonds
            spatial_network = NG.switch_chain!(spatial_network, switched_chain)
            spatial_network = NG.relax_network_keating!(
                spatial_network,
                switched_chain,
                evolution_dict)
            E_end=NG.get_total_energy_keating(spatial_network)
            dE=E_end-E_start
            T_melt=-dE/log(acceptance_probability)
            append!(T_melt_vec,T_melt)
        end
        
        result=(network_type, nr_vertices, bond_bending_const, 
            theta_ground_state, acceptance_probability, 
            relax_globally_after_threshold_cycle, shell_nr, 
            minimum(T_melt_vec))
        lock(data_lock)
        try
            push!(data, result)
        finally
            # Always unlock the lock, even if an error occurs
            unlock(data_lock)
        end
    end

    sort!(data, by = x -> (x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]))
    for (nt, nv, b, t, acceptance_probability, 
        relax_globally_after_threshold_cycle, shell_nr, T) in data
        println("[\"$nt\", $nv, $b, $t, $acceptance_probability, $relax_globally_after_threshold_cycle, $shell_nr, $T ],")
    end

    return
end


"""
Convert a network with periodic boundary conditions into a network without
periodic boundary conditions.
"""
function convert_periodic_to_non_periodic(spatial_network)

    # convert to a spatial_network_without periodic boundaries
    spatial_network_no_pbc = NG.cut_bonds_out_of_supercell!(
        deepcopy(spatial_network);
        vector_out_of_supercell_length = 1)

    # get the minimal
    min_vertex_coords, max_vertex_coords = get_min_max_vertex_coords(
                spatial_network_no_pbc)

    new_supercell_edge_length = (
        maximum(max_vertex_coords .- min_vertex_coords) * 2)

    spatial_network_no_pbc[]["supercell_edge_length"] = (
        new_supercell_edge_length)

    # shift all vertices by half the new supercell edge length
    for vertex in MetaGraphsNext.labels(spatial_network_no_pbc)
        spatial_network_no_pbc[vertex]["position"] .+= (
            new_supercell_edge_length/4 .- min_vertex_coords)
    end

    return spatial_network_no_pbc
end
