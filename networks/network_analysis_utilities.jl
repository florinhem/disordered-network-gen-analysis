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
    theta = atan( sqrt(cartesian_vec[1]^2 + cartesian_vec[2]^2), cartesian_vec[3] )
    phi = atan(cartesian_vec[2], cartesian_vec[1] )

    return [r_length, theta, phi]
end


"""
Average y(x) over different graph realizations. The y_arr contains
values for the same graph and different x along the 1 dimension
and for the same x and different graphs along the 2 dimension
"""
function get_multiple_graph_average_vec(
    y_arr::Array{Real})

    # get the number of x values, for which y(x) was evaluated
    nr_x_values = length(y_arr[:,1])

    # initialize vector for average values of y
    y_average_vec = Vector{Measurements.Measurement{Float64}}(undef, nr_x_values)

    # loop through index corresponding to different x values and average
    # over different realizations
    for x_index in 1:nr_x_values

        y_average_vec[x_index] = Measurements.measurement(
            Statistics.mean(y_arr[x_index,:]),
            Statistics.std(y_arr[x_index,:])
        )
    end

    return y_average_vec
end


"""
Average the quantity y over different graph realizations.
The y_vec contains this quantity for different graphs
"""
function get_multiple_graph_average(y_vec::Vector{Real})

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
function convert_q_l_dict_to_vec(q_l_dict::Dict, 
                                l_max::Int)

    # initialize vector for q_l values with same data type as dictionary
    # which could be Measurements.Measurement{Float64} or Float64
    q_l_vec = Vector{typeof(q_l_dict[1])}(undef, l_max+1)

    # loop through all l values and store q_l values in vector
    for l in 0:l_max

        q_l_vec[l+1] = q_l_dict[l]
    end

    return q_l_vec
end


"""
Calculate the average nr of monte carlo steps for quenching for all evolution dicts.
I got the result 13.7 +- 6.74 for 216 vertices
"""
function get_monte_carlo_steps_for_quenching_vec(evolution_dicts_directory_path::String; nr_runs = 5)

    nr_monte_carlo_steps_for_quenching_vec = Vector{Float64}()

    total_nr_dicts = 5*136
    counter = 0

    for i in 1:nr_runs

        current_path = evolution_dicts_directory_path*"run_"*string(i)*"\\"

        # get all files in directory
        filenames = readdir(current_path)

        # get filenames of all evolution dicts
        filenames_evolution_dicts = filter(filename -> endswith(filename, "_evolution.h5"), filenames)

        for filename in filenames_evolution_dicts

            println("Progress: "*string(counter/total_nr_dicts*100)*" %")
            counter += 1

            # load evolution dict
            evolution_dict = GU.load_h5_dict(current_path*filename)

            evolution_dict["nr_monte_carlo_steps_per_temperature_vec"]

            nr_monte_carlo_moves_per_step = 18*evolution_dict["nr_vertices"]

            nr_monte_carlo_moves_before_quenching = nr_monte_carlo_moves_per_step*sum(evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][1:end-1])

            nr_monte_carlo_moves_for_quenching = length(evolution_dict["move_accepted_vec"])-nr_monte_carlo_moves_before_quenching

            nr_monte_carlo_steps_for_quenching = nr_monte_carlo_moves_for_quenching/nr_monte_carlo_moves_per_step

            if nr_monte_carlo_steps_for_quenching > 49
                println(i, filename)
            end

            push!(nr_monte_carlo_steps_for_quenching_vec, nr_monte_carlo_steps_for_quenching)
        end

    end

    sort!(nr_monte_carlo_steps_for_quenching_vec)

    return nr_monte_carlo_steps_for_quenching_vec
end


"""
For given data paths and filename, calculate all order metrics, which can be
calculated from small length scales
"""
function get_small_length_scale_order_metrics(filename::String,
    graph_path::String,
    analysis_data_path::String;
    l_max_steinhardt_q_l::Int64 = 12,
    diamond_structure_factor_peak_wavenumber::Float64 = 4.676810,
    diamond_structure_factor_peak_std::Float64 = 8.5573,
    diamond_spectral_density_peak_wavenumber::Float64 = 4.680517,
    diamond_spectral_density_peak_std::Float64 = 521.88398,
    save_result = false,
    )

    # load graph
    graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename)

    # get total keating energy of final network
    total_keating_energy = graph_dict["total_energy"]

    # Measure the standard deviation of bond lengths
    bond_length_std, bond_length_vec = get_bond_length_std(graph_dict)

    # Measure the standard deviation of bond angles
    bond_angle_std, bond_angle_vec = get_bond_angle_std(graph_dict)

    # Measure the standard deviation of dihedral angles
    dihedral_angle_std, dihedral_angle_vec = get_dihedral_angle_std(graph_dict)

    # get Steinhardt local bond order parameters and store them in a vector
    q_l_total_network_mean_dict = get_q_l_total_network_mean_dict(graph_dict,
    l_max_steinhardt_q_l)

    q_l_vec = convert_q_l_dict_to_vec(q_l_total_network_mean_dict, 
        l_max_steinhardt_q_l)

    # load correlation functions
    correlation_functions_dict = GU.load_h5_dict(analysis_data_path*filename*"_correlation_functions.h5")

    # get cluster metric
    cluster_metric = get_cluster_metric(correlation_functions_dict)

    # load angle averaged structure factor
    structure_factor_angle_averaged_dict = GU.load_h5_dict(analysis_data_path*filename*"_structure_factor_angle_averaged.h5")

    # get anisotropy metric from structure factor
    anisotropy_metric_from_structure_factor = get_anisotropy_metric_from_structure_factor(
        structure_factor_angle_averaged_dict,
        diamond_structure_factor_peak_wavenumber,
        diamond_structure_factor_peak_std
    )

    # load angle averaged spectral density
    spectral_density_angle_averaged_dict = GU.load_h5_dict(analysis_data_path*filename*"_spectral_density_angle_averaged.h5")

    # get anisotropy metric from spectral density
    anisotropy_metric_from_spectral_density = get_anisotropy_metric_from_spectral_density(
        spectral_density_angle_averaged_dict,
        diamond_spectral_density_peak_wavenumber,
        diamond_spectral_density_peak_std
    )

    # create dict to save
    small_scale_order_metrics_dict = Dict(
        "total_keating_energy" => total_keating_energy,
        "bond_length_std" => bond_length_std,
        "bond_angle_std" => bond_angle_std,
        "dihedral_angle_std" => dihedral_angle_std,
        "q_l_vec" => q_l_vec,
        "cluster_metric" => cluster_metric,
        "anisotropy_metric_from_structure_factor" => anisotropy_metric_from_structure_factor,
        "anisotropy_metric_from_spectral_density" => anisotropy_metric_from_spectral_density
    )

    if save_result
        GU.save_dict_to_h5(small_scale_order_metrics_dict, analysis_data_path*filename*"_small_scale_order_metrics.h5")
    end

    return small_scale_order_metrics_dict
end