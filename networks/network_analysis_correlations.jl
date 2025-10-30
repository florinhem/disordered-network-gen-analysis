"""
these functions can be used to characterize networks by means of order metrics
measuring correlations
"""


"""
Get array of wavevectors along positive z direction for which the structure
factor is calculated when isotropy is not assumed. Each wavevector coordinate
is spanned by k_j = 2*pi*n_j/L_j where n_j is an integer within
-maximal_wavevector_int*supercell_edge_length <= n_j 
<= maximal_wavevector_int*supercell_edge_length for j = 1,2 and
0 <= n_j <= maximal_wavevector_int*supercell_edge_length for j = 3
"""
function get_wavevector_array_positive_z(
    spatial_network::MetaGraphsNext.MetaGraph;
    maximal_wavevector_int::Int64 = 5,
    periodic_boundary_conditions::Bool = true)

    if periodic_boundary_conditions
        network_edge_length = spatial_network[]["supercell_edge_length"]
    else
        min_vertex_coords, max_vertex_coords = get_min_max_vertex_coords(
            spatial_network)
        network_edge_length = minimum(max_vertex_coords .- min_vertex_coords)
    end

    # get nr of wavevectors per positive direction
    nr_wavevectors_per_pos_direction = Int(ceil(maximal_wavevector_int
        *network_edge_length))

    # generate empty array for wavevectors positions
    wavevector_array_positive_z = Array{Float64}(undef, 
        2*nr_wavevectors_per_pos_direction + 1,
        2*nr_wavevectors_per_pos_direction + 1,
        nr_wavevectors_per_pos_direction + 1,
        spatial_network[]["nr_dimensions"])

    # loop through Cartesian Indices
    for i in CartesianIndices(wavevector_array_positive_z)

        # fill wavevector array depending on dimension of current index
        wavevector_array_positive_z[i] =  (
            2*pi/network_edge_length*(
                (i[1] - nr_wavevectors_per_pos_direction - 1)* ==(i[4], 1) 
                + (i[2] - nr_wavevectors_per_pos_direction - 1)* ==(i[4], 2)
                + (i[3] - 1)* ==(i[4], 3)))

    end

    return wavevector_array_positive_z
end


"""
Evaluate a separable N-D Hann window at a point. For each dimension `d`, uses a
1D Hann window supported on `[min_coords[d], max_coords[d]]`:
w_d(x) = 0.5 * (1 - cospi(2 * t)), where t = (x - xmin) / (xmax - xmin), 
for x ∈ [xmin, xmax], and w_d(x) = 0 otherwise.
Returns the product ∏_d w_d(position_vec[d]).
"""
function hann_apodization_fct(position_vec::Vector{Float64},
    min_coords::Vector{Float64},
    max_coords::Vector{Float64},
    unused_parameter::Float64
    )::Float64

    length(position_vec) == length(min_coords) == length(max_coords) ||
        throw(ArgumentError("All vectors must have the same length."))

    window_fct_value = 1.0
    n = length(position_vec)

    @inbounds @simd for d in 1:n
        xmin = min_coords[d]
        xmax = max_coords[d]
        xmin < xmax || throw(ArgumentError(
            "min_coords[$d] must be < max_coords[$d]."))

        x = position_vec[d]
        if x < xmin || x > xmax
            return 0.0
        end

        # Normalized position in [0,1]
        t = (x - xmin) / (xmax - xmin)

        # 1D Hann in this dimension
        wd = 0.5 * (1.0 - cospi(2.0 * t))

        # Early exit if we hit an exact zero (e.g., exactly at a boundary)
        if wd == 0.0
            return 0.0
        end

        window_fct_value *= wd
    end

    # for each dimension, a normalization factor of ?? needs to be multiplied
    # TODO: find correct normalization factor
    normalization_factor = (8/3)^n
    window_fct_value *= normalization_factor

    return window_fct_value
end


"""
Hamming window: w_d(x) = 0.54 - 0.46 * cospi(2t), where 
t = (x - xmin) / (xmax - xmin), for x ∈ [xmin, xmax], and w_d(x) = 0 otherwise.
Less broadening than Hann, but more sidelobes.
"""
function hamming_apodization_fct(position_vec::Vector{Float64},
    min_coords::Vector{Float64},
    max_coords::Vector{Float64},
    unused_parameter::Float64
    )::Float64

    length(position_vec) == length(min_coords) == length(max_coords) ||
        throw(ArgumentError("All vectors must have the same length."))

    window_fct_value = 1.0
    n = length(position_vec)

    @inbounds @simd for d in 1:n
        xmin, xmax = min_coords[d], max_coords[d]
        x = position_vec[d]
        if x < xmin || x > xmax
            return 0.0
        end
        t = (x - xmin) / (xmax - xmin)
        wd = 0.54 - 0.46 * cospi(2.0 * t)
        window_fct_value *= wd
    end

    # for each dimension, a normalization factor of 1/0.54 needs to be 
    # multiplied
    normalization_factor = (1/0.54)^n
    window_fct_value *= normalization_factor

    return window_fct_value
end


"""
Gaussian window: w_d(x) = exp(-0.5 * ((t - 0.5) / sigma)^2),
with sigma controlling trade-off (smaller sigma → stronger tapering).
"""
function gaussian_apodization_fct(position_vec::Vector{Float64},
    min_coords::Vector{Float64},
    max_coords::Vector{Float64},
    sigma::Float64=0.5
    )  # default: fairly broad

    length(position_vec) == length(min_coords) == length(max_coords) ||
        throw(ArgumentError("All vectors must have the same length."))

    window_fct_value = 1.0
    n = length(position_vec)

    @inbounds @simd for d in 1:n
        xmin, xmax = min_coords[d], max_coords[d]
        x = position_vec[d]
        if x < xmin || x > xmax
            return 0.0
        end
        t = (x - xmin) / (xmax - xmin)
        wd = exp(-0.5 * ((t - 0.5)/sigma)^2)
        window_fct_value *= wd
    end

    # for each dimension, a normalization factor needs to be multiplied which
    # is the 1/(integral of the 1D Gaussian window over [-1/2, 1/2] for the
    # given sigma). 1/0.855624 is the normalization factor for sigma=0.5
    normalization_factor = (1/0.855624)^n
    window_fct_value *= normalization_factor

    return window_fct_value
end


"""
Tukey apodization window in nD (product of 1D Tukey windows per dimension).

Args:
  position_vec, min_coords, max_coords : vectors of same length (n dims)
Kwargs:
  alpha::Float64 = 0.5   # fraction in cosine tapers 
  (0 -> rectangular, 1 -> Hann)

Returns:
  window value at the given position (Float64)
"""
function tukey_apodization_fct(position_vec::Vector{Float64},
    min_coords::Vector{Float64},
    max_coords::Vector{Float64},
    alpha::Float64 = 0.5
    )

    length(position_vec) == length(min_coords) == length(max_coords) ||
        throw(ArgumentError("All vectors must have the same length."))

    # Clamp alpha gently to a sensible range to avoid surprises
    alpha = max(0.0, min(alpha, 1.0))

    window_fct_value = 1.0
    n = length(position_vec)

    @inbounds @simd for d in 1:n
        xmin, xmax = min_coords[d], max_coords[d]
        x = position_vec[d]
        if x < xmin || x > xmax
            return 0.0
        end
        t = (x - xmin) / (xmax - xmin)  # t ∈ [0,1]

        wd = if alpha <= 0.0
            1.0                            # rectangular
        elseif alpha >= 1.0
            0.5 * (1.0 - cos(2π * t))      # Hann
        else
            if t < alpha/2
                0.5 * (1.0 + cospi(2t/alpha - 1))
            elseif t <= 1 - alpha/2
                1.0
            else
                0.5 * (1.0 + cospi(2t/alpha - 2/alpha + 1))
            end
        end

        window_fct_value *= wd
    end

    # normalize using the 1D area = 1 - alpha/2  (analytical)
    normalization_factor = (1.0 / (1.0 - alpha/2))^n
    window_fct_value *= normalization_factor

    return window_fct_value
end


"""
Measure structure factor as a function of wavevector using the scattering
intensity estimator as described in equation 24 of 10.1007/s11222-023-10219-1.
Here, only the vertices of the network are considered.
"""
function get_structure_factor(
    wavevector::Vector{Float64},
    spatial_network::MetaGraphsNext.MetaGraph;
    periodic_boundary_conditions::Bool=true,
    apodization_fct = hann_apodization_fct,
    apodization_fct_parameter::Float64 = 0.5,
    min_coords=[spatial_network[]["supercell_edge_length"]/4,
        spatial_network[]["supercell_edge_length"]/4,
        spatial_network[]["supercell_edge_length"]/4],
    max_coords=[3/4*spatial_network[]["supercell_edge_length"],
        3/4*spatial_network[]["supercell_edge_length"],
        3/4*spatial_network[]["supercell_edge_length"]])

    # initialize the sum of the scattering field
    scattering_field_sum = 0.0 + 0.0*im
    
    # perform sum over all vertices
    for vertex in MetaGraphsNext.labels(spatial_network)

        # get vertex position
        vertex_pos = spatial_network[vertex]["position"]

        # calculate structure factor contribution of current vertex and
        # wavevector
        scattering_field = exp(-im*LinearAlgebra.dot(wavevector, 
                vertex_pos))
        
        if !periodic_boundary_conditions
            # apply apodization window to reduce effects of sharp edges
            apodization_window_value = apodization_fct(
                vertex_pos, min_coords, max_coords, apodization_fct_parameter)

            scattering_field *= apodization_window_value
        end

        scattering_field_sum += scattering_field
    end

    # calculate structure factor
    structure_factor = 1/spatial_network[]["nr_vertices"] * abs2(
        scattering_field_sum)

    return structure_factor
end


"""
Measure structure factor as a function of wavevector using the scattering
intensity estimator as described in equation 24 of 10.1007/s11222-023-10219-1
and considering the bonds of the network.
"""
function get_structure_factor_bonds(
    wavevector::Vector{Float64},
    spatial_network::MetaGraphsNext.MetaGraph;
    periodic_boundary_conditions::Bool=true,
    apodization_fct = hann_apodization_fct,
    apodization_fct_parameter::Float64 = 0.5,
    min_coords=[spatial_network[]["supercell_edge_length"]/4,
        spatial_network[]["supercell_edge_length"]/4,
        spatial_network[]["supercell_edge_length"]/4],
    max_coords=[3/4*spatial_network[]["supercell_edge_length"],
        3/4*spatial_network[]["supercell_edge_length"],
        3/4*spatial_network[]["supercell_edge_length"]],
    no_div_by_zero_value = 1e-10)

    # initialize the sum of the scattering field
    scattering_field_sum = 0.0 + 0.0*im
    
    # perform sum over all bonds
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # get the mid-point of the bond considering periodic boundary
        # conditions
        bond_vector = spatial_network[bond...]["vector"] 
        bond_mid_point = (((spatial_network[bond[1]]["position"] 
                .+ bond_vector ./ 2) 
            .% spatial_network[]["supercell_edge_length"]))

        # get scalar product of the wavevector with the bond mid-point and the 
        # bond vector
        scalar_prod_mid_point = LinearAlgebra.dot(wavevector, bond_mid_point)
        scalar_prod_vector = LinearAlgebra.dot(wavevector, bond_vector)
    
        # calculate structure factor contribution of current bond and 
        # wavevector
        scattering_field = (
            2/(scalar_prod_vector + no_div_by_zero_value)
            * exp(-im*scalar_prod_mid_point)
            * sin(scalar_prod_vector/2))

        if !periodic_boundary_conditions
            # apply apodization window to reduce effects of sharp edges
            apodization_window_value = apodization_fct(
                bond_mid_point, min_coords, max_coords, 
                apodization_fct_parameter)

            scattering_field *= apodization_window_value
        end

        scattering_field_sum += scattering_field
    end

    nr_bonds = length(MetaGraphsNext.edge_labels(spatial_network))

    # calculate structure factor
    structure_factor = 1/nr_bonds * abs2(scattering_field_sum)

    return structure_factor
end



"""
Measure structure factor for an array of wavevectors using the scattering
intensity estimator as described in equation 24 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_by_wavevector_array(
    spatial_network::MetaGraphsNext.MetaGraph;
    consider_bonds::Bool = false,
    maximal_wavevector_int::Int64 = 5,
    periodic_boundary_conditions::Bool = true,
    apodization_fct = hann_apodization_fct,
    apodization_fct_parameter::Float64 = 0.5,
    wavevector_array_positive_z::Array{Float64} 
        = get_wavevector_array_positive_z(spatial_network; 
            maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=periodic_boundary_conditions),
    save_result = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing,
    print_progress::Bool = false,
    thread_nr::Int64 = 0,
    print_lock = Threads.ReentrantLock())

    # initialize structure factor array
    structure_factor_array = Array{Float64}(undef, 
        size(wavevector_array_positive_z)[1:3]...)

    # get min and max coordinates of the network to apply apodization in the 
    # case of non-periodic boundary conditions
    min_vertex_coords, max_vertex_coords = get_min_max_vertex_coords(
        spatial_network)
    # modify the min and max coords to create the largest possible cubic window
    # that fully contains network
    min_vertex_coords = fill(maximum(min_vertex_coords), 3)
    max_vertex_coords = fill(minimum(max_vertex_coords), 3)

    if print_progress
        # get number of wavevectors to sample
        nr_sampled_wavevectors = size(wavevector_array_positive_z)[1] *
            size(wavevector_array_positive_z)[2] *
            size(wavevector_array_positive_z)[3]

        # initialize progress counter
        wavevector_count = 0

        if consider_bonds
            bonds_string = "bonds"
        else
            bonds_string = "vertices"
        end
    end

    # get structure factor for all wavevectors
    for i in 1:size(wavevector_array_positive_z)[1], 
        j in 1:size(wavevector_array_positive_z)[2], 
        k in 1:size(wavevector_array_positive_z)[3]

        # get structure factor either for bonds or for vertices of the network
        if consider_bonds
            # get structure factor for current wavevector and bond network
            structure_factor_array[i,j,k] = get_structure_factor_bonds(
                wavevector_array_positive_z[i,j,k,:], spatial_network;
                periodic_boundary_conditions=periodic_boundary_conditions,
                apodization_fct=apodization_fct,
                apodization_fct_parameter=apodization_fct_parameter,
                min_coords=min_vertex_coords,
                max_coords=max_vertex_coords)

        else
            # get structure factor for current wavevector and spatial network
            structure_factor_array[i,j,k] = get_structure_factor(
                wavevector_array_positive_z[i,j,k,:], spatial_network;
                periodic_boundary_conditions=periodic_boundary_conditions,
                apodization_fct=apodization_fct,
                apodization_fct_parameter=apodization_fct_parameter,
                min_coords=min_vertex_coords,
                max_coords=max_vertex_coords)
        end

        # print progress
        if print_progress
            wavevector_count += 1

            # print every 100th voxel
            if wavevector_count % 10000 == 0
                progress_percentage = (wavevector_count
                    /nr_sampled_wavevectors*100)

                lock(print_lock) do
                    Format.printfmtln("Structure factor calculation for {1:s} 
                    thread nr {2:d}: {3:.1f} %", bonds_string, thread_nr, 
                    progress_percentage)
                end
            end
        end
    end

    # extend structure factor to negative z direction
    structure_factor_array = cat(dims=3, 
        structure_factor_array[end:-1:1,end:-1:1,end:-1:2],
        structure_factor_array)

    # extend wavevector array to negative z direction
    wavevector_array_negative_z = reverse(
        wavevector_array_positive_z[:,:,2:end,:], dims=3)
    wavevector_array_negative_z[:,:,:,3] .*= (-1)
    wavevector_array = cat(dims=3, 
        wavevector_array_negative_z, wavevector_array_positive_z)

    # get vector of vector of wavenumbers along all directions
    wavenumber_vec_vec = [wavevector_array[:,1,1,1], 
                          wavevector_array[1,:,1,2], 
                          wavevector_array[1,1,:,3]]

    # create dict to save
    structure_factor_dict = Dict{String, Any}(
        "wavevector_array" => wavevector_array, 
        "wavenumber_vec_vec" => wavenumber_vec_vec,
        "structure_factor_array" => structure_factor_array)

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_dict["label"] = label
    end

    # save results if desired
    if save_result
        if consider_bonds
            bonds_string = "_bonds"
        else
            bonds_string = ""
        end
        GU.save_dict_to_h5(copy(structure_factor_dict), 
                save_path*"_structure_factor"*bonds_string*"_array.h5")
    end
                                            
    return structure_factor_dict
end


"""
From the dictionary containing the 3d structure factor, get a dictionary that,
for each index vector squared contains all values of the structure factor. The 
index vector squared is proportional to the wavenumber squared. The returned 
dictionary can be used to calculate the angle averaged structure factor 
"""
function get_structure_factor_by_index_vec_squared_dict(
    structure_factor_dict::Dict;
    consider_spectral_density::Bool = false)

    if consider_spectral_density
        data_type_string = "spectral_density"
    else
        data_type_string = "structure_factor"
    end

    structure_factor_by_index_vec_squared_dict = Dict{Int64, Vector{Float64}}()

    # determine origin vector of cartesian coordinates
    index_vector_origin = Int.((size(
        structure_factor_dict[data_type_string*"_array"]) .+ 1 ) ./ 2)

    # get lattice constant of reciprocal lattice
    reciprocal_lattice_constant = LinearAlgebra.norm(
        structure_factor_dict["wavevector_array"][(index_vector_origin 
            .+ [1,0,0])...,:])

    # loop through Cartesian Indices
    for i in CartesianIndices(structure_factor_dict[data_type_string*"_array"])

        # determine actual index vector
        index_vector = i.I .- index_vector_origin

        # calculate length squared of index vector
        index_vector_length_squared = sum(index_vector.^2)

        # check if key exists in dictionary
        if index_vector_length_squared in keys(
            structure_factor_by_index_vec_squared_dict)

            # add structure factor to dictionary
            push!(structure_factor_by_index_vec_squared_dict[
                    index_vector_length_squared], 
                abs(structure_factor_dict[data_type_string*"_array"][i]))

        else

            # add key and structure factor to dictionary
            structure_factor_by_index_vec_squared_dict[index_vector_length_squared]=(
                [abs(structure_factor_dict[data_type_string*"_array"][i])])

        end

    end

    return [structure_factor_by_index_vec_squared_dict, 
        reciprocal_lattice_constant]
end


"""
Calculate angle averaged structure factor from 3d array of structure factor or
equivalenty for the spectral density
"""
function get_structure_factor_angle_averaged(
    structure_factor_dict::Dict;
    consider_spectral_density::Bool = false,
    consider_bonds::Bool = false,
    gaussian_filter::Bool = true,
    gaussian_filter_sigma_x::Float64 = 2*pi/25, 
    gaussian_filter_filtered_data_x_step_length::Float64 = 2*pi/25,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing)

    if consider_spectral_density
        data_type_string = "spectral_density"
    else
        data_type_string = "structure_factor"
    end

    # get the dictionary that, for each index vector squared contains all 
    # values of the structure factor. The index vector squared is proportional 
    # to the square of the wavenumber
    structure_factor_by_index_vec_squared_dict, reciprocal_lattice_constant = (
        get_structure_factor_by_index_vec_squared_dict(structure_factor_dict; 
            consider_spectral_density = consider_spectral_density))
    
    # initialize wavenumber vector and structure factor vector
    unfiltered_wavenumber_vec = Vector{Float64}()
    unfiltered_structure_factor_vec = Vector{
        Measurements.Measurement{Float64}}()

    # get vector of wavenumbers and angle averaged structure factor including
    # their uncertainty
    for key in keys(structure_factor_by_index_vec_squared_dict)

        # get wavenumber from wavevector
        push!(unfiltered_wavenumber_vec, reciprocal_lattice_constant*sqrt(key))

        # get angle averaged structure factor
        push!(unfiltered_structure_factor_vec, 
            Measurements.measurement(Statistics.mean(
                structure_factor_by_index_vec_squared_dict[key]),
            Statistics.std(structure_factor_by_index_vec_squared_dict[key])))

    end

    # sort wavenumber vector and structure factor vector
    unfiltered_structure_factor_vec = unfiltered_structure_factor_vec[
            sortperm(unfiltered_wavenumber_vec)]
        sort!(unfiltered_wavenumber_vec)

    # create dict to save
    structure_factor_angle_averaged_dict = Dict{String, Any}(
        "unfiltered_wavenumber_vec" => unfiltered_wavenumber_vec, 
        "unfiltered_"*data_type_string*"_vec" 
            => unfiltered_structure_factor_vec)

    # apply gaussian filter if desired
    if gaussian_filter
        filtered_data_x, filtered_data_y = GU.gaussian_filter_1d(
            unfiltered_wavenumber_vec[2:end], 
            unfiltered_structure_factor_vec[2:end]; 
            sigma_x=gaussian_filter_sigma_x, 
            filtered_data_x_step_length
                =gaussian_filter_filtered_data_x_step_length)

        structure_factor_angle_averaged_dict["wavenumber_vec"] = (
            filtered_data_x)
        structure_factor_angle_averaged_dict[data_type_string*"_vec"] = (
            filtered_data_y)
        structure_factor_angle_averaged_dict["gaussian_filter_sigma_x"] = (
            gaussian_filter_sigma_x)
    end

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_angle_averaged_dict["label"] = label
    end

    if save_result
        if consider_bonds
            bonds_string = "_bonds"
        else
            bonds_string = ""
        end
        GU.save_dict_to_h5(copy(structure_factor_angle_averaged_dict),
            save_path*"_"*data_type_string*bonds_string*"_angle_averaged.h5")
    end
                                            
    return structure_factor_angle_averaged_dict
end


function get_correlation_functions_from_vertex_distances(
    distance_vec::Vector{Float64},
    nr_vertices,
    maximal_vertex_distance::Float64,
    vertex_density::Float64;
    distance_histogram_bin_width::Float64 = 0.02)

    # get histogram of distance vector
    distance_histogram = StatsBase.fit(
        StatsBase.Histogram, distance_vec, 
        0.0:distance_histogram_bin_width:maximal_vertex_distance, 
        closed=:left)

    vertex_distance_vec = collect(distance_histogram.edges[1][2:end] 
        .- distance_histogram_bin_width/2)
    vertex_nr_vec = distance_histogram.weights

    # calculate pair correlaion function
    # (eq. 3 in 10.1016/j.physrep.2018.03.001)
    pair_correlation_fct_vec = ((1/(4*pi*vertex_density * nr_vertices) ) 
        .* vertex_nr_vec ./ vertex_distance_vec.^2) .* (
            2/distance_histogram_bin_width)

    # set the first entry of the pair correlation function to zero, to correct
    # for weird behavior of the histogram for crystalline networks
    pair_correlation_fct_vec[1] = 0.0

    # calculate total correlation function
    # (eq. 4 in 10.1016/j.physrep.2018.03.001)
    total_correlation_fct_vec = pair_correlation_fct_vec .- 1

    # calculate cumulative coordination number 
    # (eq. 5 in 10.1016/j.physrep.2018.03.001)
    cumulative_coord_nr_vec = (4*pi*vertex_density*distance_histogram_bin_width
        * cumsum(vertex_distance_vec.^2 .* pair_correlation_fct_vec ))

    return [vertex_distance_vec, cumulative_coord_nr_vec, 
        pair_correlation_fct_vec, total_correlation_fct_vec]
end


"""
Calculate pair correlation function, total correlation function and
cumulative coordination number as defined in equations 3-5 of
10.1016/j.physrep.2018.03.001
"""
function get_correlation_functions(
    spatial_network::MetaGraphsNext.MetaGraph;
    distance_histogram_bin_width::Float64 = 0.02,
    periodic_boundary_conditions::Bool = true,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing)

    # create vector of all vertex positions
    vertex_position_mat = Matrix{Float64}(undef, 3, 
        spatial_network[]["nr_vertices"])

    for i in 1:spatial_network[]["nr_vertices"]
        vertex_position_mat[:, i] =spatial_network[i]["position"]
    end

    # get vector of distance between all pairs of vertices by considering
    # periodic boundary conditions
    distance_vec = Vector{Float64}(undef, Int(spatial_network[]["nr_vertices"]
        *(spatial_network[]["nr_vertices"]-1)/2))

    # also get vector of distances between all uncoordinated pairs of vertices
    distance_uncoordinated_vec = Vector{Float64}(undef, 0)

    current_index = 1

    for i in 1:spatial_network[]["nr_vertices"]-1
        neighbors = collect(MetaGraphsNext.neighbor_labels(spatial_network, i))

        for j in i+1:spatial_network[]["nr_vertices"]
            if periodic_boundary_conditions
                distance_vec[current_index] = LinearAlgebra.norm(
                    NG.get_distance_vector_pbc(
                        vertex_position_mat[:, i], vertex_position_mat[:, j], 
                        spatial_network[]["supercell_edge_length"]))
                
                if j ∉ neighbors
                    push!(distance_uncoordinated_vec, 
                        distance_vec[current_index])
                end
            else
                distance_vec[current_index] = LinearAlgebra.norm(
                    vertex_position_mat[:, i] 
                    - vertex_position_mat[:, j])

                if j ∉ neighbors
                    push!(distance_uncoordinated_vec, 
                        distance_vec[current_index])
                end
            end

            current_index += 1
        end
    end

    if periodic_boundary_conditions
        maximal_vertex_distance = spatial_network[]["supercell_edge_length"]/2
    else
        min_vertex_coords, max_vertex_coords = get_min_max_vertex_coords(
            spatial_network)
        maximal_vertex_distance = minimum(max_vertex_coords 
            .- min_vertex_coords)/2
    end

    # get_vertex density
    if periodic_boundary_conditions
        vertex_density = (spatial_network[]["nr_vertices"] 
            / spatial_network[]["supercell_edge_length"]^3)
    else
        vertex_density = spatial_network[]["nr_vertices"] / (
            (max_vertex_coords[1] - min_vertex_coords[1])
            *(max_vertex_coords[2] - min_vertex_coords[2])
            *(max_vertex_coords[3] - min_vertex_coords[3]))
    end

    # get correlation functions for all vertex distances
    (vertex_distance_vec, cumulative_coord_nr_vec, pair_correlation_fct_vec, 
        total_correlation_fct_vec) = (
        get_correlation_functions_from_vertex_distances(
            distance_vec, spatial_network[]["nr_vertices"], 
            maximal_vertex_distance, vertex_density; 
            distance_histogram_bin_width=distance_histogram_bin_width))

    # get correlation functions for uncoordinated vertex distances
    coordination_nr_mean, coordination_nr_std = get_coordination_nr_statistics(
        spatial_network)
    nr_vertices_uncoordinated = (
        spatial_network[]["nr_vertices"] - coordination_nr_mean)
    uncoordinated_vertex_density = (vertex_density 
        * nr_vertices_uncoordinated / spatial_network[]["nr_vertices"])
    (uncoordinated_vertex_distance_vec, cumulative_coord_nr_uncoordinated_vec, 
        pair_correlation_fct_uncoordinated_vec, 
        total_correlation_fct_uncoordinated_vec) = (
        get_correlation_functions_from_vertex_distances(
            distance_uncoordinated_vec, nr_vertices_uncoordinated, 
            maximal_vertex_distance, uncoordinated_vertex_density; 
            distance_histogram_bin_width=distance_histogram_bin_width))

    # create dict to save
    correlation_functions_dict = Dict{String, Any}(
        "vertex_distance_vec" => vertex_distance_vec, 
        "cumulative_coord_nr_vec" => cumulative_coord_nr_vec,
        "pair_correlation_fct_vec" => pair_correlation_fct_vec,
        "total_correlation_fct_vec" => total_correlation_fct_vec,
        "vertex_density" => vertex_density,
        "cumulative_coord_nr_uncoordinated_vec" 
            => cumulative_coord_nr_uncoordinated_vec,
        "pair_correlation_fct_uncoordinated_vec" 
            => pair_correlation_fct_uncoordinated_vec,
        "total_correlation_fct_uncoordinated_vec" 
            => total_correlation_fct_uncoordinated_vec,
        "uncoordinated_vertex_density" => uncoordinated_vertex_density
        )

    # add label to dictionary if label is not nothing
    if label !== nothing
        correlation_functions_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(correlation_functions_dict),
            save_path*"_correlation_functions.h5")
    end

    return correlation_functions_dict
end


"""
Define a vertex_homogeneity metric as the average distance to the nearest
distance to the next vertex independently whether it is a neighbor or not. This
distance is smaller than 1 if the network is clustered
"""
function get_vertex_homogeneity_metric(correlation_functions_dict::Dict;
    consider_uncoordinated_vertices::Bool = false)

    # get the smallest vertex distance where the cumulative coordination number
    # is above 1
    if consider_uncoordinated_vertices
        cumulative_coord_nr = (
            correlation_functions_dict[
                "cumulative_coord_nr_uncoordinated_vec"])
    else
        cumulative_coord_nr = (
            correlation_functions_dict["cumulative_coord_nr_vec"])
    end
    argmin_index = findfirst(cumulative_coord_nr .> 1)
    vertex_homogeneity_metric = (
        correlation_functions_dict["vertex_distance_vec"][argmin_index])

    return vertex_homogeneity_metric
end


"""
Define an anisotropy metric based on the structure factor or the 
spectral density that ranges between 0 and 1 (high anisotropy). The metric is
based on the coefficient of variation (std over mean) of the structure factor
as a function of direction. In reality, values around 0.45 represent low
anisotropy and values around 0.55 represent high anisotropy.
"""
function get_anisotropy_metric_from_structure_factor(
    structure_factor_angle_averaged_dict::Dict;
    maximal_length_to_check = 3.0,
    nr_closest_wavenumbers = 3,
    normalization_parameter = 1.0)

    # set the wavenumbers where structure factor will be checked
    wavenumbers_to_check_vec = (2*pi) ./ collect(
        0.5:0.5:maximal_length_to_check+0.01)

    # for each wavenumber to check, calculate a normalized coefficient of
    # variation (std over mean)
    anisotropy_metric_vec = Vector{Float64}(undef, 
        length(wavenumbers_to_check_vec))

    for i in eachindex(wavenumbers_to_check_vec)
        # get the three indices of the closest wavenumbers to the current
        # wavenumber to check
        structure_factor_indices_vec = sortperm(abs.(
            structure_factor_angle_averaged_dict["unfiltered_wavenumber_vec"]
            .- wavenumbers_to_check_vec[i]))[1:nr_closest_wavenumbers]

        # get the structure factor for the three wavenumbers
        structure_factor_vec = structure_factor_angle_averaged_dict[
                "unfiltered_structure_factor_vec"][structure_factor_indices_vec[1]]
        for j in 2:nr_closest_wavenumbers
            structure_factor_vec = vcat(structure_factor_vec, 
                structure_factor_angle_averaged_dict[
                    "unfiltered_structure_factor_vec"][
                        structure_factor_indices_vec[j]])
        end

        # calculate the coefficient of variation
        mean_structure_factor = Statistics.mean(structure_factor_vec)
        coefficient_of_variation = (
            Measurements.uncertainty(mean_structure_factor)
            /Measurements.value(mean_structure_factor))

        # calculate the anisotropy metric that ranges between 0 and 1
        anisotropy_metric_vec[i] = coefficient_of_variation/(
            normalization_parameter + coefficient_of_variation)
    end

    # calculate the average of the anisotropy entropy metric
    anisotropy_metric = Statistics.mean(anisotropy_metric_vec)

    return anisotropy_metric
end


"""
Get local number variance for a spherical window of given radius from the
structure factor according to eq 58 and 60 in 10.1016/j.physrep.2018.03.001
"""
function get_local_nr_variance(
    spatial_network::MetaGraphsNext.MetaGraph,
    structure_factor_dict::Dict,
    window_radius;
    periodic_boundary_conditions::Bool = true)

    # check if system is 3d
    if spatial_network[]["nr_dimensions"] != 3
        @error "Calculation of number variance is only implemented for 3d
            systems"
    end

    # get_vertex density
    if periodic_boundary_conditions
        vertex_density = (spatial_network[]["nr_vertices"] 
            / spatial_network[]["supercell_edge_length"]^spatial_network[][
                "nr_dimensions"])
    else
        min_vertex_coords, max_vertex_coords = get_min_max_vertex_coords(
            spatial_network)
        vertex_density = (spatial_network[]["nr_vertices"] / (
            (max_vertex_coords[1] - min_vertex_coords[1])
            *(max_vertex_coords[2] - min_vertex_coords[2])
            *(max_vertex_coords[3] - min_vertex_coords[3])))
    end

    # calculate wavenumber sampling step length
    wavenumber_step_length = (structure_factor_dict["wavenumber_vec"][2] 
        - structure_factor_dict["wavenumber_vec"][1])

    # calculate local nr variance
    local_nr_variance = (vertex_density * 32 * pi^2 * window_radius^2
        * wavenumber_step_length * sum(
            structure_factor_dict["structure_factor_vec"] 
            .* sin.( structure_factor_dict["wavenumber_vec"] 
                .* window_radius) .^2
            ./ structure_factor_dict["wavenumber_vec"].^2
        )    
    )
    
    return local_nr_variance
end


"""
Measure local nr variance as a function of window radius from structure factor 
according to equations 47 and 58 in 10.1016/j.physrep.2018.03.001
"""
function get_local_nr_variance_by_window_radius_vec(
    spatial_network,
    structure_factor_angle_averaged_dict::Dict;
    window_radius_step_length::Float64 = 0.1,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\random_networks\sample_name",
    label = nothing)

    # get vector of sphere radii
    if periodic_boundary_conditions
        max_window_radius = spatial_network[]["supercell_edge_length"]/4
    else
        min_vertex_coords, max_vertex_coords = get_min_max_vertex_coords(
            spatial_network)
        max_window_radius = minimum(max_vertex_coords 
            .- min_vertex_coords)/4
    end

    sphere_radius_vec = collect(window_radius_step_length:
        window_radius_step_length:max_window_radius)

    # initialize local nr variance vector
    local_nr_variance_vec = Vector{Measurements.Measurement{Float64}}(undef, 
        length(sphere_radius_vec))

    # get local nr variance vector as a window radius
    for i in eachindex(sphere_radius_vec)
        local_nr_variance_vec[i] = get_local_nr_variance(spatial_network,
        structure_factor_angle_averaged_dict,
        sphere_radius_vec[i])

    end

    # get local nr variance over sphere volume which indicates hyperuniformity
    # at large window radii
    local_nr_variance_over_sphere_volume_vec = (
        local_nr_variance_vec ./ (4/3*pi*sphere_radius_vec.^3))

    local_nr_variance_dict = Dict("sphere_radius_vec" => sphere_radius_vec,
        "local_nr_variance_vec" => local_nr_variance_vec,
        "local_nr_variance_over_sphere_volume_vec" => local_nr_variance_over_sphere_volume_vec)

    if label !== nothing
        local_nr_variance_dict["label"] = label
    end

    if save_result
        GU.save_dict_to_h5(copy(local_nr_variance_dict), 
            save_path*"_local_nr_variance.h5")
    end

    return local_nr_variance_dict
end


"""
Get the excess spreadability for a given t value according to eq 23 of
10.1103/PhysRevE.104.054102.
"""
function get_excess_spreadability(
    structure_factor_array_bonds_dict::Dict,
    time::Float64;
    min_wavenumber_to_consider::Float64 = 0.0,
    consider_spectral_density::Bool = false)

    if consider_spectral_density
        data_type_string = "spectral_density"
    else
        data_type_string = "structure_factor"
    end

    wavevector_array = structure_factor_array_bonds_dict["wavevector_array"]

    # only consider wavevectors above a certain threshold length to avoid 
    # artifacts, e. g. from the apodization
    wavenumber_array = norms = mapslices(x -> LinearAlgebra.norm(x),
        structure_factor_array_bonds_dict["wavevector_array"]; dims=4)
    mask = wavenumber_array .> min_wavenumber_to_consider

    wavevector_volume_element = prod(
        wavevector_array[2,2,2,:] .- wavevector_array[1,1,1,:])

    # calculate excess spreadability according to eq 17 of 
    # 10.1103/PhysRevE.104.054102
    excess_spreadability = 1/wavevector_volume_element * sum(
        (structure_factor_array_bonds_dict[
            data_type_string*"_array"] .*  exp.(-wavenumber_array.^2 .* time)
            )[mask])

    return excess_spreadability
end



"""
Get the exponent alpha, that determines the scaling of the structure factor
or the spectral density with the wavenumber k for k -> 0 as in eq 25 of
10.1103/PhysRevE.109.064108. If this exponent is 0, the system is
nonhyperuniform, if it is >0, the system is hyperuniform. Three classes of
hyperuniformity are defined: 0 < alpha < 1, alpha = 1 and alpha > 1 as written
below eq 21 of 10.1103/PhysRevE.109.064108.
"""
function get_hyperuniformity_alpha(structure_factor_array_bonds_dict::Dict;
    consider_spectral_density::Bool = false,
    t_range = (1e-1, 1),
    min_wavenumber_to_consider::Float64 = 0.0)

    # calculate the excess spreadability for each t value
    excess_spreadability_values = [get_excess_spreadability(
        structure_factor_array_bonds_dict, time; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) 
        for time in t_range]

    # get logarithm of time and excess spreadability
    log_time_values = log10.(t_range)
    log_excess_spreadability_values = log10.(excess_spreadability_values)

    # get the slope of the log-log curve from the two data points
    # (log(time), log(excess_spreadability))
    slope = (
        (log_excess_spreadability_values[2] 
            - log_excess_spreadability_values[1]) 
        / (log_time_values[2] - log_time_values[1]))

    # get alpha from the minimal slope according to section IIIB of 
    # 10.1103/PhysRevE.109.064108
    hyperuniformity_alpha = -2*slope - 3
    
    return hyperuniformity_alpha
end
