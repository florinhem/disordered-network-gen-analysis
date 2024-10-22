"""
Functions to calculate the autocovariance function and the spectral
density for 3d data assuming statistical isotropy and only working with 
absolute values of distance and momentum

This code requires the following Packages:
import Statistics   # for statistical operations like mean()
import Measurements # for handling data with uncertainty and error propagation
"""

"""
Check if coordinate lies within data array
"""
function coordinate_is_inside_array(coordinate::Tuple, size_data::Tuple )

    # this variable stores the boolean whether coodinate lies within data array
    inside_array = true

    # loop through dimensions to check coordinates
    for i in 1:3

        if !(coordinate[i] in 1:size_data[i])
            inside_array = false
        end

    end

    return inside_array
end



"""
Get 2 random coordinates that lie inside the structure at a given distance 
to each other
"""
function get_random_coordinates(sampling_distance, size_data::Tuple )

    # get first random coordinate x1
    x1_vec = Vector{Int64}(undef, 3)

    for i in 1:3
        x1_vec[i] = rand(1:size_data[i])

    end

    # convert vector to tuple
    x1 = Tuple(x1_vec)

    # get second random coordinate at distance r from x1 by getting a 
    # random direction according to the equations 1 and 2 in this website
    # https://mathworld.wolfram.com/SpherePointPicking.html

    theta = 2*pi*rand()
    phi = acos( 2 * rand() - 1 )

    x2 = Int.( round.( x1 .+ sampling_distance .* ( 
        cos(phi) *sin(theta),
        sin(phi)*sin(theta),
        cos(theta) )  ) ) 
    

    # if x2 lies within data array, return the two coordinates
    if coordinate_is_inside_array(x2, size_data)
        return [x1, x2]

    # otherwise run the function again
    else
        return get_random_coordinates(sampling_distance, size_data)

    end

end


"""
Get the autocovariance function for statistically homogeneous and isotropic
3d media
"""
function get_autocovariance_fct_isotrope(
    sampling_distance,
    size_data::Tuple,
    volume_fract_tot::Float64,
    data_binary::Array{Bool}
    ;
    nr_measurements_per_distance::Int64 = 10000)

    # initialize vector from which the two point prob. fct. 
    # will be calculated later
    two_point_prob_fct_summand_vec = Vector{Float64}(undef, 
        nr_measurements_per_distance)

    # generate given number of random vectors to sample 2 point prob. function
    for i in 1:nr_measurements_per_distance

        # get two random coordinates at given sampling distance to each other
        x1, x2 = get_random_coordinates(sampling_distance, size_data)

        # calculate the contribution to the two point prob. fct. 
        # from these coodinates
        two_point_prob_fct_summand_vec[i] = data_binary[x1...] * data_binary[x2...]

    end

    # calculate ensemble averaged 2 point prob. function
    # the error is calculated according to the equation given at
    # https://en.wikipedia.org/wiki/Checking_whether_a_coin_is_fair
    # below the sentence
    # "In statistics, the estimate of a proportion of a sample "
    # (denoted by p) has a standard error given by
    two_point_prob_fct = Measurements.measurement( 
        Statistics.mean( two_point_prob_fct_summand_vec ),
        sqrt( volume_fract_tot^2 * (1 - volume_fract_tot^2) 
                / nr_measurements_per_distance  ) )

    # determine autocovariance function
    autocovariance_fct = two_point_prob_fct - volume_fract_tot^2

    return autocovariance_fct
    
end


"""
get the number of sampling distances
"""
function get_nr_sampling_distances(mean_edge_length_data::Int64)

    # get the length of the vector with entries 
    # 0:1:Int(round(mean_edge_length_data/2))
    nr_sampling_distances = Int(round(mean_edge_length_data/2)) + 1

    return nr_sampling_distances
    
end


"""
get a vector of sampling distances at which the 2 point prob. function and 
then the autocovariance function are measured. Based on this vector the 
discrete Fourier transformation will be determined to calculate the 
spectral density
"""
function get_sampling_distance_vec(
    mean_edge_length_data::Int64; 
    nr_sampling_distances::Int64 = 
        get_nr_sampling_distances(mean_edge_length_data) )

    # get vector of sampling distances
    sampling_distance_vec = collect(LinRange(
        0, 
        Int(round(mean_edge_length_data/2)), 
        nr_sampling_distances)) 

    return sampling_distance_vec

end


"""
determine autocovariance function as a function of sampling distance for 
isotrope 3d binary media
"""
function get_autocovariance_fct_isotrope_by_sampling_distance_vec(
    structure_dict::Dict
    ;
    nr_sampling_distances::Int64 = 
        get_nr_sampling_distances(structure_dict["mean_edge_length_data"] ),
    nr_measurements_per_distance::Int64 = 10000,
    sampling_distance_vec = 
        get_sampling_distance_vec(
            structure_dict["mean_edge_length_data"]; 
            nr_sampling_distances = nr_sampling_distances),
    save_result = false,
    save_path = raw"..\analysis_data\sample_name",
    voxel_edge_length = nothing,
    label = nothing
    )

    # create vector where for each sampling distance the 
    # autocovariance function and its uncertainty will be stored
    autocovariance_fct_vec = Vector{Measurements.Measurement}(undef, 
        nr_sampling_distances)

    # for each sampling distance get autocovariance function and its uncertainty
    for i in eachindex(sampling_distance_vec)

        # get vector of local volume fractions and the number of 
        # independent samples
        autocovariance_fct_vec[i] = get_autocovariance_fct_isotrope(
            sampling_distance_vec[i],
            structure_dict["size_data"],
            structure_dict["volume_fract_tot"],
            structure_dict["data_binary"];
            nr_measurements_per_distance = nr_measurements_per_distance)

        # print current calculation status
        println("sampling distance "*string(sampling_distance_vec[i])*" done")

    end

    # create dictionary for current plot
    autocovariance_fct_dict = Dict(
        "sampling_distance_vec" => sampling_distance_vec,
        "autocovariance_fct_vec" => autocovariance_fct_vec,
        "nr_measurements_per_distance" => nr_measurements_per_distance,
        "voxel_edge_length" => structure_dict["voxel_edge_length"] ,
        "label" => structure_dict["label"])

    # if desired, adjust voxel edge length and label
    autocovariance_fct_dict = modify_keys_in_dict(
        autocovariance_fct_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        # save the plot_dict to a H5 file
        GU.save_dict_to_h5(copy(autocovariance_fct_dict), save_path*
            "_autocovariance_fct.h5")

    end

    return autocovariance_fct_dict
    
end


"""
Calculate the spectral density for statistically homogeneous and isotropic
3d media
"""
function get_spectral_density_isotrope(
    autocovariance_fct_vec::Vector,
    wavenumber::Float64;
    sampling_distance_vec = get_sampling_distance_vec(
        mean_edge_length_data; 
        nr_sampling_distances = nr_sampling_distances),
    sampling_distance_cutoff = 1000,
    voxel_edge_length = 10)

    # get sampling distance step length 
    # (inverval between adjacent sampling distances)
    sampling_distance_step_length = sampling_distance_vec[2] - 
        sampling_distance_vec[1]

    # find index of sampling distance that is the closest to the cutoff
    sampling_distance_cutoff_index = argmin( abs.(sampling_distance_vec 
        .- (sampling_distance_cutoff/voxel_edge_length ) ) )

    # from sampling distances and corresponding autocovariance_fcts, 
    # calculate spectral density, which is the Fourier transform of the 
    # autocovariance fct
    spectral_density = (4*pi*sampling_distance_step_length/wavenumber 
        * sum( sampling_distance_vec[1:sampling_distance_cutoff_index]
            .* autocovariance_fct_vec[1:sampling_distance_cutoff_index]
            .* sin.( wavenumber 
                .* sampling_distance_vec[1:sampling_distance_cutoff_index] ) ) )

    return spectral_density

end


"""
get wavenumber vector such that the discrete Fourier transform can be 
determined properly based on the sampled distances.
The info about minimal and maximal wavenumbers is taken from the DFT intro PDF
"""
function get_wavenumber_vec(sampling_distance_vec::Vector)
    
    # determine the length between two adjacent sampling distances
    sampling_distance_step_length = sampling_distance_vec[2] - 
        sampling_distance_vec[1]

    # determine the nr of sampling distances
    nr_sampling_distances = length( sampling_distance_vec )

    # determine fundamental wavenumber
    fundamental_wavenumber = 2*pi/(sampling_distance_step_length*
        nr_sampling_distances)

    # determine sampled wavenumbers vec
    wavenumber_vec = collect(1:floor(nr_sampling_distances/2)) .* 
        fundamental_wavenumber

    return wavenumber_vec
    
end



"""
determine spectral density as a function of the wavenumber k for 
isotrope 3d binary media.
The argument sampling_distance_cutoff is the sampling distance in units of nm, 
up to which the autocovariance function will be Fourier transformed to yield 
the spectral density.
This cutoff is necessary, because above a certain threshold, the 
autocovariance function is dominated by noise, which falsifies the 
Fourier transform
"""
function get_spectral_density_isotrope_by_wavenumber_vec(
    structure_dict::Dict
    ;
    nr_sampling_distances::Int64 = 
        get_nr_sampling_distances(structure_dict["mean_edge_length_data"]),
    nr_measurements_per_distance::Int64 = 10000,
    sampling_distance_vec = 
        get_sampling_distance_vec(structure_dict["mean_edge_length_data"]; 
            nr_sampling_distances = nr_sampling_distances),
    autocovariance_fct_dict = 
        get_autocovariance_fct_isotrope_by_sampling_distance_vec(
            structure_dict;
            nr_sampling_distances = nr_sampling_distances,
            nr_measurements_per_distance = nr_measurements_per_distance,
            sampling_distance_vec = sampling_distance_vec),
    sampling_distance_cutoff = 1000,
    save_result = false,
    save_path = raw"..\analysis_data\sample_name",
    voxel_edge_length = nothing,
    label = nothing)

    # set voxel edge length if not specified in the arguments
    if voxel_edge_length === nothing
        voxel_edge_length = autocovariance_fct_dict["voxel_edge_length"]
    end

    # get vector of window edge lengths that will be measured
    wavenumber_vec = get_wavenumber_vec(
        autocovariance_fct_dict["sampling_distance_vec"] )

    # determine number of wavenumbers for which Fourier transform is calculated
    nr_wavenumbers = length(wavenumber_vec)

    # create vector where for each wavenumber the spectral density and its
    # uncertainty will be stored
    spectral_density_vec = Vector{Measurements.Measurement}(undef, 
        nr_wavenumbers)

    # for each wavenumber get spectral density and its uncertainty
    for i in eachindex(wavenumber_vec)

        # get vector of local volume fractions and the number of 
        # independent samples
        spectral_density_vec[i] = get_spectral_density_isotrope(
            autocovariance_fct_dict["autocovariance_fct_vec"],
            wavenumber_vec[i];
            sampling_distance_vec=autocovariance_fct_dict[
                "sampling_distance_vec"],
            sampling_distance_cutoff = sampling_distance_cutoff,
            voxel_edge_length = voxel_edge_length)

        # print current calculation status
        println("wavenumber "*string(wavenumber_vec[i])*" done")

    end

    # create dictionary for current plot
    spectral_density_dict = Dict(
        "wavenumber_vec" => wavenumber_vec,
        "spectral_density_vec" => spectral_density_vec,
        "nr_measurements_per_distance" => 
                autocovariance_fct_dict["nr_measurements_per_distance"],
        "sampling_distance_cutoff" => sampling_distance_cutoff,
        "voxel_edge_length" => voxel_edge_length,
        "label" => autocovariance_fct_dict["label"] )

    # if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(
        spectral_density_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict), save_path*
            "_spectral_density.h5")

    end

    return spectral_density_dict
    
end