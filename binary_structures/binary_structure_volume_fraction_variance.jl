"""
These functions are used to analyse binary structure data by scanning it with 
spherical or cubic measuring windows and measuring the local volume fraction
variance. The way this variance scales as a function of the measuring window
size is a measure of hyperuniformity.
"""



"""
for a given window center and dimensions and the total data dimensions
determine index ranges of the data array that lie within the window
"""
function get_window_index_domains(nr_dimensions_data::Int64,
                                    edge_length_window::Int64, random_window_center::Tuple)

        # create empty arrays where upper and lower index bounds of the measuring window
        # on the data array will be stored
        lower_indices_vector = zeros(Int64, nr_dimensions_data)
        upper_indices_vector = zeros(Int64, nr_dimensions_data)

        # for each coordinate calculate lower and upper index bound
        for i in 1:nr_dimensions_data
            lower_indices_vector[i] = random_window_center[i] - Int( (edge_length_window - 1)/2  )
            upper_indices_vector[i] = random_window_center[i] + Int( (edge_length_window - 1)/2  )
    
        end

    return [lower_indices_vector, upper_indices_vector]

end



"""
determine index coordinates that lie within a sphere of given center position and diameter
"""
function get_coord_inside_window_vec(nr_dimensions_data::Int64,
                                    edge_length_window::Int64, 
                                    window_center::Tuple;
                                    window_shape::String="spherical")

    # create vector where coordinates will be stored
    coord_inside_window_vec = []

    # get upper and lower indices of cubic domain where sphere lies within
    lower_indices_vector, upper_indices_vector = get_window_index_domains(nr_dimensions_data,
                                                                        edge_length_window,
                                                                        window_center)

    # check if data is 3d, otherwise this function does not work
    if nr_dimensions_data !== 3
        @error "Determination of coodinates within sphere only works for 3d data."
    end

    # loop through dimensions and determine coordinates within window
    for i in lower_indices_vector[1]:upper_indices_vector[1]
        for j in lower_indices_vector[2]:upper_indices_vector[2]
            for k in lower_indices_vector[3]:upper_indices_vector[3]

                # differetiate between spherical and cubic windows
                if cmp(window_shape, "spherical") == 0

                    # calculate squared distance of coordinate from window center
                    sq_window_center_distance = ((i-window_center[1])^2 
                                                + (j-window_center[2])^2
                                                + (k-window_center[3])^2  )

                    # check if coordinate lies whithin sphere
                    if sq_window_center_distance < (edge_length_window/2)^2
                        push!(coord_inside_window_vec, (i,j,k))

                    end

                elseif cmp(window_shape, "cubic") == 0
                    push!(coord_inside_window_vec, (i,j,k))

                else
                    @error "Window shape is neither spherical nor cubic."
                end

            end
        end
    end

    return coord_inside_window_vec
end


"""
get window volume depending on its shape
"""
function get_window_volume(nr_dimensions_data::Int64, 
                            edge_length_window::Int64;
                            window_shape::String="spherical")
    
    # check if structure is 3d, otherwise this function does not work
    if nr_dimensions_data !== 3
        @error "Number of independent samples can currently
                 only be calculated for 3d structures"
    end

    # calculate window volume depending on its shape
    if cmp(window_shape, "spherical") == 0
        window_volume = 4/3 * pi * (edge_length_window/2)^3 

    elseif cmp(window_shape, "cubic") == 0
        window_volume = edge_length_window^nr_dimensions_data 

    else
        @error "Window shape is neither spherical nor cubic."
    end

    return window_volume

end


"""
calculate the number of independent samples from the array of sampled voxels
"""
function get_nr_independent_samples(nr_dimensions_data::Int64, 
                                        edge_length_window::Int64, 
                                        sampled_voxels_array::Array{Int64};
                                        window_shape::String="spherical")
    
    # get window volume depending on its shape
    window_volume = get_window_volume(nr_dimensions_data, 
                                        edge_length_window;
                                        window_shape=window_shape)

    # calculate number of independent samples
    nr_independent_samples = 1/window_volume * sum(sampled_voxels_array)

    return nr_independent_samples

end


"""
determine the number of measurements per window size depending on positioning and shape
"""
function get_nr_measurements(nr_dimensions_data::Int64, size_data::Tuple,
                            edge_length_window::Int64; 
                            constraints_nr_measurements::Tuple{Int64, Int64} = (100,10000),
                            window_positioning::String="random",
                            window_shape::String="spherical")

    # differetiate between random and scanned window positioning
    if cmp(window_positioning, "random") == 0

        # determine volume of data and window
        volume_data = prod(size_data)
        volume_window =  get_window_volume(nr_dimensions_data, 
                                            edge_length_window;
                                            window_shape=window_shape)

        # calculate unconstrained number of measurements
        nr_measurements_unconstrained = Int(floor( volume_data/(2*volume_window) ) )

        # constrain number measurements by given lower and upper bounds
        nr_measurements = minimum( [ maximum( [ nr_measurements_unconstrained, 
                                                constraints_nr_measurements[1] ] ), 
                                    constraints_nr_measurements[2] ] )  

    elseif cmp(window_positioning, "scanned") == 0

        # differetiate between cubic and spherical window window_shape
        if cmp(window_shape, "spherical") == 0

            # spherical windows are placed such that they overlap by one radius
            nr_measurements = Int( prod( 2 .* floor.( size_data ./ (edge_length_window + 1) )  .-1 ) )

        elseif cmp(window_shape, "cubic") == 0

            # cubic windows are placed one next to the other
            nr_measurements = Int( prod( floor.( size_data ./ edge_length_window ) ) )

        else
            @error "Window shape is neither spherical nor cubic."
        end

    else
        @error "Window positioning is neither random nor scanned."
    end

    return nr_measurements

end


"""
perform a single measurement on local volume fraction and update the sampled voxels array
"""
function get_local_volume_fract_and_update_sampled_voxels!(nr_dimensions_data::Int64,
                                                            edge_length_window::Int64,
                                                            data_binary::Array{Bool},
                                                            sampled_voxels_array::Array{Int64},
                                                            window_center::Tuple;
                                                            window_shape::String="spherical" )

    # determine index coordinates that lie within the window depending on its shape
    coord_inside_window_vec = get_coord_inside_window_vec(nr_dimensions_data,
                                                        edge_length_window, 
                                                        window_center;
                                                        window_shape = window_shape)

    # in this variable the nr of voxels with value 1 will be stored
    nr_1_voxels = 0

    # sum over data entries within window and continuously update sampled voxels array
    for coord_tuple in coord_inside_window_vec

        # count data entry at current coordinate within window
        nr_1_voxels += data_binary[coord_tuple...]

        # update sampled voxels array
        sampled_voxels_array[coord_tuple...] = 1

    end

    # determine local volume fraction in window array
    local_volume_fract = nr_1_voxels / length(coord_inside_window_vec)
    
    return [local_volume_fract, sampled_voxels_array]
    
end


"""
get random window center position
"""
function get_random_window_center(nr_dimensions_data::Int64, size_data::Tuple,
                                    edge_length_window::Int64)

    # determine lower bounds on coordinates
    lower_bound = Int(ceil(edge_length_window/2))

    # create vector where random coordinates will be stored
    random_window_center = zeros(Int64, nr_dimensions_data)

    # assign a random number to each window center coordinate
    for i in 1:nr_dimensions_data

        # get upper bounds on coordinates
        upper_bound = size_data[i] - lower_bound + 1

        random_window_center[i] = rand((lower_bound:upper_bound))

    end

    return Tuple(random_window_center)

end


"""
get vector of window centers depending on their positioning and shape
"""
function get_window_centers_vec(nr_dimensions_data::Int64, 
                                edge_length_window::Int64,
                                size_data::Tuple,
                                nr_measurements::Int64;
                                window_positioning::String="random",
                                window_shape::String="spherical")

    # initialize vector where window centers will be stored in
    window_centers_vec = Vector{Tuple}(undef, nr_measurements)

    # differetiate between random and scanned window positioning
    if cmp(window_positioning, "random") == 0

        # get random window centers 
        for i in eachindex(window_centers_vec)
            window_centers_vec[i] = get_random_window_center(nr_dimensions_data, 
                                                            size_data,
                                                            edge_length_window)
        end

    elseif cmp(window_positioning, "scanned") == 0

        if nr_dimensions_data !==3
            @error "Scanned window centers are only implemented for 3d data."
        end

        # differetiate between cubic and spherical window window_shape
        if cmp(window_shape, "spherical") == 0

            # determine nr of measurements along all dimensions
            # the spherical windows are shifted such that the centers of neighboring windows do not overlap 
            # therefore +1 is added to the window edge length
            nr_measurements_along_dimensions = Int.( 2 .* floor.( size_data ./ (edge_length_window+1) ) .-1 )

            # count the current nr of measurement window
            current_nr_measurement = 1

            # spherical windows are placed such that they overlap by one radius
            for i in 1:nr_measurements_along_dimensions[1]
                for j in 1:nr_measurements_along_dimensions[2]
                    for k in 1:nr_measurements_along_dimensions[3]

                        window_centers_vec[current_nr_measurement] = ( i * Int( (edge_length_window+1)  /2),
                                                                    j * Int( (edge_length_window+1) /2),
                                                                    k * Int( (edge_length_window+1) /2) )

                        # update the current nr of measurement window
                        current_nr_measurement += 1
                    end
                end
            end

        elseif cmp(window_shape, "cubic") == 0

            # cubic windows are placed one next to the other

            # determine nr of measurements along all dimensions
            nr_measurements_along_dimensions = Int.( floor.( size_data ./ edge_length_window ) )

            # count the current nr of measurement window
            current_nr_measurement = 1

            # spherical windows are placed such that they overlap by one radius
            for i in 0:(nr_measurements_along_dimensions[1] - 1 )
                for j in 0:(nr_measurements_along_dimensions[2] - 1 )
                    for k in 0:(nr_measurements_along_dimensions[3] - 1 )

                        window_centers_vec[current_nr_measurement] = ( Int( (1/2 + i) * edge_length_window + 1/2 ),
                                                                        Int( (1/2 + j) * edge_length_window + 1/2 ),
                                                                        Int( (1/2 + k) * edge_length_window + 1/2 ) )

                        # update the current nr of measurement window
                        current_nr_measurement += 1
                    end
                end
            end

        else
            @error "Window shape is neither spherical nor cubic."
        end

    else
        @error "Window positioning is neither random nor scanned."
    end

    return window_centers_vec

end


"""
for a given window size, shape and positioning get the local volume fraction 
    variance by using windows of random position
"""
function get_local_volume_fract_vec(nr_dimensions_data::Int64, 
                                            edge_length_window::Int64,
                                            size_data::Tuple, 
                                            data_binary::Array{Bool};
                                            window_positioning::String="random",
                                            window_shape::String="spherical",
                                            constraints_nr_measurements::Tuple{Int64, Int64} = (100,10000) ) 

    # create auxiliary array which records all voxels that have been sampled
    # by at least one measurement
    sampled_voxels_array = zeros(Int64, size_data)

    # determine number of measurements based on window positioning and shape
    nr_measurements = get_nr_measurements(nr_dimensions_data, size_data,
                                        edge_length_window; 
                                        constraints_nr_measurements = constraints_nr_measurements,
                                        window_positioning=window_positioning,
                                        window_shape=window_shape)

    # get vector of window centers depending on their positioning and shape
    window_centers_vec = get_window_centers_vec(nr_dimensions_data, 
                                                edge_length_window,
                                                size_data,
                                                nr_measurements;
                                                window_positioning=window_positioning,
                                                window_shape=window_shape)


    # create vector where local volume fractions will be stored
    local_volume_fract_vec = Vector{Float64}(undef, nr_measurements)  #  zeros(Float64, nr_measurements)

    # for each measurement, calculate local volume fraction and update sampled voxels array
    for i in 1:nr_measurements

        # perform a single measurement on local volume fraction and update the sampled voxels array
        local_volume_fract, sampled_voxels_array = get_local_volume_fract_and_update_sampled_voxels!(nr_dimensions_data,
                                                                                                edge_length_window,
                                                                                                data_binary,
                                                                                                sampled_voxels_array,
                                                                                                window_centers_vec[i];
                                                                                                window_shape=window_shape )

        # save local volume fraction
        local_volume_fract_vec[i] = local_volume_fract

    end
    
    # get number of independent samples from the array of sampled voxels
    nr_independent_samples = get_nr_independent_samples(nr_dimensions_data, 
                                                    edge_length_window, 
                                                    sampled_voxels_array)


    return [local_volume_fract_vec, nr_independent_samples]

end



"""
get a vector with odd measuring window edge lengths
"""
function get_window_edge_length_vec(mean_edge_length_data::Int64; nr_window_sizes::Int64 = 100 )

    # set maximal window edge length
    max_window_edge_length = Int( floor( mean_edge_length_data/4 ) )

    # determine by which even number the window edge length should be incremented at each step
    # depending on the size of the total edge length and the desired number of window sizes
    window_edge_length_increment = maximum([2, 
                                    Int64( 2*ceil( mean_edge_length_data/(4*nr_window_sizes) ) ) ] )  

    # create vector with measuring window edge lengths according to determined increment
    window_edge_length_vec = collect(3:window_edge_length_increment:max_window_edge_length)

    return window_edge_length_vec

end



"""
from vector of local volume fractions, get their variance
"""
function get_local_volume_fract_variance(volume_fract_tot::Float64,
                                        local_volume_fract_vec::Vector{Float64} )

    # calculate local volume fraction variance
    local_volume_fract_variance = Statistics.mean( (local_volume_fract_vec .- volume_fract_tot ) .^2 )

    return local_volume_fract_variance

end



"""
determine uncertainty on the local volume fraction using the equation given in this discussion:
https://math.stackexchange.com/questions/72975/variance-of-sample-variance
"""
function get_uncertainty_local_volume_fract_variance(nr_measurements,
                                                    volume_fract_tot::Float64,
                                                    local_volume_fract_vec::Vector{Float64},
                                                    local_volume_fract_variance::Float64)

    # first, calculate the fourth central moment of the local volume fraction
    fourth_central_moment = Statistics.mean( (local_volume_fract_vec .- volume_fract_tot ) .^4 )

    # calculate uncertainty on the local volume fraction
    uncertainty_local_volume_fract_variance = sqrt( fourth_central_moment / nr_measurements
                                                    - local_volume_fract_variance^2
                                                        *(nr_measurements-3)
                                                        /(nr_measurements*(nr_measurements-1))
                                                    )

    return uncertainty_local_volume_fract_variance

end



"""
determine local volume fraction variance as a function of measuring window size.
This can be done for ordered or spherical windows at random or scanned positions.
"""
function get_local_volume_fract_variance_by_window_vec(structure_dict::Dict;
        nr_window_sizes::Int64 = 100,
        window_positioning::String="random",
        window_shape::String="spherical",
        constraints_nr_measurements::Tuple{Int64, Int64} = (100,10000),
        save_result = false,
        save_path = raw"..\analysis_data\sample_name",
        voxel_edge_length = nothing,
        label = nothing)

    # get vector of window edge lengths that will be measured
    window_edge_length_vec = get_window_edge_length_vec(structure_dict["mean_edge_length_data"]; nr_window_sizes )

    # create vector where for each window size the local volume fraction variance and its
    # uncertainty will be stored
    local_volume_fract_variance_vec = zeros(Measurements.Measurement, length(window_edge_length_vec))

    # for each measuring window size get local volume fraction variance and its uncertainty
    for i in eachindex(window_edge_length_vec)

        # get vector of local volume fractions and the number of independent samples
        local_volume_fract_vec, nr_independent_samples = get_local_volume_fract_vec(
                                        structure_dict["nr_dimensions_data"],
                                        window_edge_length_vec[i],
                                        structure_dict["size_data"], 
                                        structure_dict["data_binary"];
                                        window_positioning=window_positioning,
                                        window_shape=window_shape,
                                        constraints_nr_measurements=constraints_nr_measurements )

        # determine variance of local volume fraction
        local_volume_fract_variance = get_local_volume_fract_variance(structure_dict["volume_fract_tot"],
                                                                    local_volume_fract_vec )

        # determine uncertainty on variance of local volume fraction using the number of independent samples
        uncertainty_local_volume_fract_variance = get_uncertainty_local_volume_fract_variance(
                                                                    nr_independent_samples,
                                                                    structure_dict["volume_fract_tot"],
                                                                    local_volume_fract_vec,
                                                                    local_volume_fract_variance)

        # save variance of local volume fraction and its uncertainty
        local_volume_fract_variance_vec[i] = Measurements.measurement(local_volume_fract_variance, 
                                                                    uncertainty_local_volume_fract_variance)

        # print current calculation status
        println("window length "*string(window_edge_length_vec[i])*" done")

    end

    # create dictionary for current plot
    local_volume_fract_variance_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "window_positioning" => window_positioning,
                    "window_shape" => window_shape,
                    "voxel_edge_length" => voxel_edge_length,
                    "label" => label )

    # if desired, adjust voxel edge length and label
    local_volume_fract_variance_dict = modify_keys_in_dict(local_volume_fract_variance_dict, 
                                                            voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(local_volume_fract_variance_dict), 
            save_path*"_volume_fraction_variance.h5")

    end

    return local_volume_fract_variance_dict
    
end



