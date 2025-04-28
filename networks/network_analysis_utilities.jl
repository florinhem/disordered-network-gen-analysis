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
    evolution_dicts_directory_path::String; nr_runs = 5)

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
            # TODO convert 18 (hardcoded) into 4*3*3/2 with coordination_nr

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
    structure_dict_path::String,
    analysis_data_path::String;
    bond_radius = 0.35,
    voxel_edge_length = 0.1,
    pore_size_sampling_grid_size = 0.2,
    print_progress::Bool = false,
    print_lock = Threads.ReentrantLock())

    # get network
    spatial_network = NG.load_spatial_network_from_gml(
        spatial_network_path*filename*".gml")

    # create structure dict
    structure_dict = get_binary_data_from_spatial_network(
        spatial_network;
        bond_radius = bond_radius,
        voxel_edge_length = voxel_edge_length,
        save_path = structure_dict_path*filename,
        filename = filename,
        save_result=true)

    # get ring size distribution
    ring_size_distribution_dict = get_ring_size_distribution(
        spatial_network;
        save_result = true,
        save_path = analysis_data_path*filename)

    # get ring radius distribution
    ring_radius_distribution_dict = get_ring_radius_distribution(
        spatial_network,
        ring_size_distribution_dict;
        save_result = true,
        save_path = analysis_data_path*filename)

    # get autocovariance function as a function of direction
    autocovariance_fct_direction_dict = (
        get_autocovariance_fct_by_sampling_indices_array(
            structure_dict;
            save_result = true,
            save_path = analysis_data_path*filename,
            print_progress = print_progress,
            thread_nr = Threads.threadid() ,
            print_lock = print_lock))

    # get spectral density by wavevector array
    spectral_density_dict = get_spectral_density_by_wavevector_array_fft(
        structure_dict;
        save_autocovariance_fct_direction_dict = false,
        save_result = true,
        save_path = analysis_data_path*filename,
        autocovariance_fct_direction_dict = autocovariance_fct_direction_dict)

    # get angle averaged spectral density
    spectral_density_angle_averaged_dict = get_spectral_density_angle_averaged(
        spectral_density_dict;
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = true,
        save_path = analysis_data_path*filename)

    # get volume fraction variance
    volume_fract_variance_dict = get_volume_fract_variance(
        autocovariance_fct_direction_dict;
        save_result = true,
        save_path = analysis_data_path*filename)

    # get structure factor by wavevector array
    structure_factor_dict = get_structure_factor_by_wavevector_array(
        spatial_network;
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing)

    # get angle averaged structure factor
    structure_factor_angle_averaged_dict = get_structure_factor_angle_averaged(
        structure_factor_dict::Dict;
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing)

    # get correlation functions
    correlation_functions_dict = get_correlation_functions(
        spatial_network;
        distance_histogram_bin_width = 0.02,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing)

    # get pore size distribution
    pore_size_distribution_dict = get_pore_size_distribution(
        spatial_network;
        sampling_grid_size = pore_size_sampling_grid_size,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = nothing,
        print_progress = print_progress,
        thread_nr = Threads.threadid() ,
        print_lock = print_lock)

    # get all order metrics that contain information about small length scales
    small_scale_order_metrics_dict = get_small_length_scale_order_metrics(
        filename,
    spatial_network_path,
    analysis_data_path;
    l_max_steinhardt_q_l = 12,
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
    structure_dicts_path::String,
    analysis_data_path::String;
    bond_radius = 0.35,
    voxel_edge_length = 0.1,
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels::Int64 = 20000,
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

        get_all_dicts_from_network_single_file(run_and_filename[7:end],
    spatial_networks_path*run_and_filename[1:6],
    structure_dicts_path*run_and_filename[1:6],
    analysis_data_path*run_and_filename[1:6];
    bond_radius = bond_radius,
    voxel_edge_length = voxel_edge_length,
    structure_factor_diamond_std_value_ratio 
        = structure_factor_diamond_std_value_ratio,
    spectral_density_diamond_std_value_ratio 
        = spectral_density_diamond_std_value_ratio,
    pore_size_distribution_nr_sampled_voxels 
        = pore_size_distribution_nr_sampled_voxels,
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
    spatial_networks_path,
    structure_dicts_path,
    analysis_data_path::String;
    bond_radius = 0.35,
    voxel_edge_length = 0.1,
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels::Int64 = 20000,
    print_progress::Bool = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())

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

        # filter out filenames that are already in the analysis data path
        filenames_spatial_networks_not_analyzed = filter(filename 
            -> !isfile(analysis_data_path*"run_"*string(i)*"/"*filename
                *"_small_scale_order_metrics.h5"), filenames_spatial_networks)

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
            structure_dicts_path,
            analysis_data_path;
            bond_radius = bond_radius,
            voxel_edge_length = voxel_edge_length,
            structure_factor_diamond_std_value_ratio 
                = structure_factor_diamond_std_value_ratio,
            spectral_density_diamond_std_value_ratio 
                = spectral_density_diamond_std_value_ratio,
            pore_size_distribution_nr_sampled_voxels 
                = pore_size_distribution_nr_sampled_voxels,
            print_progress = print_progress,
            print_lock = print_lock)
    end

    return
end


"""
For given data paths and filename, calculate all order metrics, which can be
calculated from small length scales
"""
function get_small_length_scale_order_metrics(filename::String,
    network_path::String,
    analysis_data_path::String;
    l_max_steinhardt_q_l::Int64 = 12,
    save_result = false,
    )

    # load network
    spatial_network = NG.load_spatial_network_from_gml(
        network_path*filename*".gml")

    # get total keating energy of final network
    total_keating_energy = spatial_network[]["total_energy"]

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

    # get cluster metric
    vertex_homogeneity_metric = get_vertex_homogeneity_metric(
        correlation_functions_dict)

    # load pore size distribution
    pore_size_distribution_dict = GU.load_h5_dict(
        analysis_data_path*filename*"_pore_size_distribution.h5")

    # get second moment of the pore size distribution
    pore_size_distribution_second_moment = (
        get_pore_size_distribution_second_moment(pore_size_distribution_dict))

    # load angle averaged structure factor
    structure_factor_dict = GU.load_h5_dict(
        analysis_data_path*filename*"_structure_factor.h5")

    # get anisotropy metric from structure factor
    anisotropy_metric_from_structure_factor = (
        get_anisotropy_metric_from_structure_factor(structure_factor_dict))

    # load angle averaged spectral density
    spectral_density_angle_averaged_dict = GU.load_h5_dict(
        analysis_data_path*filename*"_spectral_density_angle_averaged.h5")

    # get anisotropy metric from spectral density
    anisotropy_metric_from_spectral_density = (
        get_anisotropy_metric_from_spectral_density(
            spectral_density_angle_averaged_dict;
            diamond_std_value_ratio=spectral_density_diamond_std_value_ratio))

    # create dict to save
    small_scale_order_metrics_dict = Dict(
        "total_keating_energy" => total_keating_energy,
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
        "pore_size_distribution_second_moment" 
            => pore_size_distribution_second_moment,
        "anisotropy_metric_from_structure_factor" 
            => anisotropy_metric_from_structure_factor,
        "anisotropy_metric_from_spectral_density" 
            => anisotropy_metric_from_spectral_density
    )

    if save_result
        GU.save_dict_to_h5(small_scale_order_metrics_dict, 
            analysis_data_path*filename*"_small_scale_order_metrics.h5")
    end

    return small_scale_order_metrics_dict
end


function get_small_length_scale_order_metrics_all_files(
    analysis_data_path::String;
    l_max_steinhardt_q_l::Int64 = 12,
    save_result::Bool = false,)

    # get filenames of small scale order
    # anisotropy_metric_from_spectral_density
    all_filenames = readdir(analysis_data_path)
    order_metrics_filenames = [filename for filename in all_filenames
        if occursin("small_scale_order_metrics", filename)]

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
    ring_size_mean_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    ring_size_std_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    ring_radius_mean_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    ring_radius_std_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    pore_size_distribution_second_moment_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    anisotropy_metric_from_structure_factor_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))
    anisotropy_metric_from_spectral_density_vec = Vector{Float64}(undef,
        length(order_metrics_filenames))

    # loop through order metric filenames
    for i in eachindex(order_metrics_filenames)

        # load order metrics
        order_metrics_dict = GU.load_h5_dict(
            analysis_data_path*order_metrics_filenames[i])

        # get all order metrics and save them to the corresponding vectors
        total_keating_energy_vec[i] = (
            order_metrics_dict["total_keating_energy"])
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
        ring_size_mean_vec[i] = order_metrics_dict["ring_size_mean"]
        ring_size_std_vec[i] = order_metrics_dict["ring_size_std"]
        ring_radius_mean_vec[i] = order_metrics_dict["ring_radius_mean"]
        ring_radius_std_vec[i] = order_metrics_dict["ring_radius_std"]
        pore_size_distribution_second_moment_vec[i] = (
            order_metrics_dict["pore_size_distribution_second_moment"])
        anisotropy_metric_from_structure_factor_vec[i] = (
            order_metrics_dict["anisotropy_metric_from_structure_factor"])
        anisotropy_metric_from_spectral_density_vec[i] = (
            order_metrics_dict["anisotropy_metric_from_spectral_density"])

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
    ring_size_mean_vec = ring_size_mean_vec[
        sortperm(total_keating_energy_vec)]
    ring_size_std_vec = ring_size_std_vec[
        sortperm(total_keating_energy_vec)]
    ring_radius_mean_vec = ring_radius_mean_vec[
        sortperm(total_keating_energy_vec)]
    ring_radius_std_vec = ring_radius_std_vec[
        sortperm(total_keating_energy_vec)]
    pore_size_distribution_second_moment_vec = (
        pore_size_distribution_second_moment_vec[
            sortperm(total_keating_energy_vec)])
    anisotropy_metric_from_structure_factor_vec = (
        anisotropy_metric_from_structure_factor_vec[
            sortperm(total_keating_energy_vec)])
    anisotropy_metric_from_spectral_density_vec = (
        anisotropy_metric_from_spectral_density_vec[
            sortperm(total_keating_energy_vec)])

    sort!(total_keating_energy_vec)

    # create dict to save
    order_metrics_dict = Dict(
        "total_keating_energy_vec" => total_keating_energy_vec,
        "bond_length_std_vec" => bond_length_std_vec,
        "bond_angle_std_vec" => bond_angle_std_vec,
        "dihedral_angle_entropy_vec" => dihedral_angle_entropy_vec,
        "bond_orientation_entropy_vec" => bond_orientation_entropy_vec,
        "coordination_nr_mean_vec" => coordination_nr_mean_vec,
        "coordination_nr_std_vec" => coordination_nr_std_vec,
        # TODO: "q_l_mat" => q_l_mat, creates error so far
        "vertex_homogeneity_metric_vec" => vertex_homogeneity_metric_vec,
        "ring_size_mean_vec" => ring_size_mean_vec,
        "ring_size_std_vec" => ring_size_std_vec,
        "ring_radius_mean_vec" => ring_radius_mean_vec,
        "ring_radius_std_vec" => ring_radius_std_vec,
        "pore_size_distribution_second_moment_vec" 
            => pore_size_distribution_second_moment_vec,
        "anisotropy_metric_from_structure_factor_vec" 
            => anisotropy_metric_from_structure_factor_vec,
        "anisotropy_metric_from_spectral_density_vec" 
            => anisotropy_metric_from_spectral_density_vec,
        "filenames_vec" => order_metrics_filenames
    )

    # save dict to h5 file if desired
    if save_result
        GU.save_dict_to_h5(order_metrics_dict,
        analysis_data_path*"all_order_metrics.h5")
    end

    return order_metrics_dict
end