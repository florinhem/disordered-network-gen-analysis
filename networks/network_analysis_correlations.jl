"""
these functions can be used to characterize networks by means of order metrics
measuring correlations
"""

"""
Measure structure factor for a given wavenumber, averaged over angles according
to Barlett's isotropic estimator as described in equation 40 of 
10.1007/s11222-023-10219-1
"""
function get_structure_factor_bartlett_isotrope(
    spatial_network::MetaGraphsNext.MetaGraph, 
    wavenumber)

    # check if structure is 3d
    if spatial_network[]["nr_dimensions"] !== 3
        @error "Structure factor calculation is, so far, only implemented for 
            3d."
    end

    # initialize double sum
    double_sum = 0.0

    # perform sum over all combinations of vertices
    for vertex_1 in MetaGraphsNext.labels(spatial_network)

        # get first vertex position
        vertex_1_pos = spatial_network[vertex_1]["position"]

        for vertex_2 in vertex_1+1:spatial_network[]["nr_vertices"]

            # get second vertex position
            vertex_2_pos = spatial_network[vertex_2]["position"]

            # distance between vertices
            vertex_distance = LinearAlgebra.norm(vertex_2_pos .- vertex_1_pos)

            # calculate structure factor contribution of current vertex
            # combination
            double_sum += sin(wavenumber*vertex_distance) / (
                wavenumber*vertex_distance)

        end
    end

    # calculate structure factor
    structure_factor = 1 + 2/spatial_network[]["nr_vertices"] * double_sum

    return structure_factor
end


"""
Get vector of wavenumbers, for which the structure factor is calculated
"""
function get_wavenumber_vec(
    spatial_network::MetaGraphsNext.MetaGraph;
    sampling_distance_step_length = 0.05,
    maximal_sampling_distance 
        = spatial_network[]["supercell_edge_length"]/sqrt(3))

    # determine virtual nr of sampling distances (in reality I don't sample in
    # direct space anywhere)
    nr_sampling_distances = floor(maximal_sampling_distance
        /(sampling_distance_step_length))

    # get nr of wavenumbers which is half the nr of sampling distances
    nr_wavenumbers = Int(floor(nr_sampling_distances/2))

    # get fundamental wavenumber
    fundamental_wavenumber = 2*pi/maximal_sampling_distance

    # get vector of wavenumbers
    wavenumber_vec = collect(1:nr_wavenumbers) * fundamental_wavenumber

    return wavenumber_vec
end


"""
Measure structure factor as a function of wavenumber averaged over angles
according to Barlett's isotropic estimator as described in equation 26 of 
10.1007/s11222-023-10219-1
"""
function get_structure_factor_bartlett_isotrope_by_wavenumber_vec(
    spatial_network::MetaGraphsNext.MetaGraph;
    sampling_distance_step_length = 0.1,
    maximal_sampling_distance = spatial_network[]["supercell_edge_length"]/2,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\random_networks\sample_name",
    print_progress::Bool = false,
    label = nothing)

    # get vector of wavenumbers
    wavenumber_vec = get_wavenumber_vec(spatial_network; 
        sampling_distance_step_length = sampling_distance_step_length,
        maximal_sampling_distance = maximal_sampling_distance)

    # initialize structure factor vector
    structure_factor_bartlett_vec = Vector{Float64}(undef, 
        length(wavenumber_vec))

    # get vector of structure factor as a function of wavenumber
    for i in eachindex(wavenumber_vec)
        structure_factor_bartlett_vec[i] = (
            get_structure_factor_bartlett_isotrope(
                spatial_network, wavenumber_vec[i]))

        # print progress
        if print_progress
            println("Progress: ", i/length(wavenumber_vec)*100, "%")
        end

    end

    structure_factor_bartlett_dict = Dict("wavenumber_vec" => wavenumber_vec,
        "unfiltered_structure_factor_vec" => structure_factor_bartlett_vec,
        "sampling_distance_step_length" => sampling_distance_step_length,
        "maximal_sampling_distance" => maximal_sampling_distance )

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_bartlett_dict["label"] = label
    end

    if save_result
        GU.save_dict_to_h5(copy(structure_factor_bartlett_dict),
            save_path*"_structure_factor_bartlett_isotrope.h5")

    end

    return structure_factor_bartlett_dict
end


"""
Measure structure factor as a function of wavenumber averaged over angles
according to scattering intensity estimator as described in equation 40 of
10.1007/s11222-023-10219-1
"""
function get_structure_factor_isotrope(
    spatial_network::MetaGraphsNext.MetaGraph, 
    wavenumber;
    nr_wavevector_samples::Int = 10000)

    # check if structure is 3d
    if spatial_network[]["nr_dimensions"] !== 3
        @error "Structure factor calculation is, so far, only implemented for 
            3d."
    end

    # get desired number of wavevector samples
    theta_vec = 2*pi*rand(nr_wavevector_samples)
    phi_vec = acos.(2*rand(nr_wavevector_samples) .- 1)
    wavevector_mat =  wavenumber .* stack( [sin.(phi_vec).*cos.(theta_vec), 
                                            sin.(phi_vec).*sin.(theta_vec), 
                                            cos.(phi_vec)] , 
                                            dims=1)

    # initialize structure factor sum
    structure_factor = 0.0
    
    # perform sum over all wavevector samples
    for wavevector in eachcol(wavevector_mat)

        # initialize the sum of the scattering field
        scattering_field_sum = 0.0 + 0.0*im

        # perform sum over all vertices
        for vertex in MetaGraphsNext.labels(spatial_network)

            # get vertex position
            vertex_pos = spatial_network[vertex]["position"]
        
            # calculate structure factor contribution of current vertex and wavevector
            scattering_field_sum += exp(-im*LinearAlgebra.dot(wavevector, 
                vertex_pos))

        end

        # calculate structure factor
        structure_factor = 1/spatial_network[]["nr_vertices"] * abs2(
            scattering_field_sum)

        # add structure factor to sum
        structure_factor += structure_factor/nr_wavevector_samples

    end

    return structure_factor
end


"""
Measure structure factor for a given wavenumber, averaged over angles according
to the scattering intensity estimator as described in equation 26 of
10.1007/s11222-023-10219-1
"""
function get_structure_factor_isotrope_by_wavenumber_vec(
    spatial_network::MetaGraphsNext.MetaGraph;
    sampling_distance_step_length = 0.1,
    maximal_sampling_distance = spatial_network[]["supercell_edge_length"]/2,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\random_networks\sample_name",
    print_progress::Bool = false,
    label = nothing)

    # get vector of wavenumbers
    wavenumber_vec = get_wavenumber_vec(spatial_network; 
        sampling_distance_step_length = sampling_distance_step_length,
        maximal_sampling_distance = maximal_sampling_distance)

    # initialize structure factor vector
    unfiltered_structure_factor_vec = Vector{Float64}(undef, 
        length(wavenumber_vec))

    # get vector of structure factor as a function of wavenumber
    for i in eachindex(wavenumber_vec)
        unfiltered_structure_factor_vec[i] = get_structure_factor_isotrope(
            spatial_network, wavenumber_vec[i])

        # print progress
        if print_progress
            println("Progress: ", i/length(wavenumber_vec)*100, "%")
        end

    end

    # create dictionary for current plot
    structure_factor_dict = Dict("wavenumber_vec" => wavenumber_vec,
        "unfiltered_structure_factor_vec" => unfiltered_structure_factor_vec,
        "sampling_distance_step_length" => sampling_distance_step_length,
        "maximal_sampling_distance" => maximal_sampling_distance,
        "label" => label )

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(structure_factor_dict),
            save_path*"_structure_factor_isotrope.h5")

    end

    return structure_factor_dict
end


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
    maximal_wavevector_int::Int64 = 5)

    # get nr of wavevectors per positive direction
    nr_wavevectors_per_pos_direction = Int(ceil(maximal_wavevector_int
        *spatial_network[]["supercell_edge_length"]))

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
            2*pi/spatial_network[]["supercell_edge_length"]*(
                (i[1] - nr_wavevectors_per_pos_direction - 1)* ==(i[4], 1) 
                + (i[2] - nr_wavevectors_per_pos_direction - 1)* ==(i[4], 2)))

    end

    return wavevector_array_positive_z
end


"""
Measure structure factor as a function of wavevector using the scattering
intensity estimator as described in equation 24 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor(
    wavevector::Vector{Float64},
    spatial_network::MetaGraphsNext.MetaGraph)

    # initialize the sum of the scattering field
    scattering_field_sum = 0.0 + 0.0*im
    
    # perform sum over all vertices
    for vertex in MetaGraphsNext.labels(spatial_network)

        # get vertex position
        vertex_pos = spatial_network[vertex]["position"]
    
        # calculate structure factor contribution of current vertex and
        # wavevector
        scattering_field_sum += exp(-im*LinearAlgebra.dot(wavevector, 
            vertex_pos))

    end

    # calculate structure factor
    structure_factor = 1/spatial_network[]["nr_vertices"] * abs2(
        scattering_field_sum)

    return structure_factor
end


"""
Measure structure factor for an array of wavevectors using the scattering
intensity estimator as described in equation 24 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_by_wavevector_array(
    spatial_network::MetaGraphsNext.MetaGraph;
    maximal_wavevector_int::Int64 = 5,
    wavevector_array_positive_z::Array{Float64} 
        = get_wavevector_array_positive_z(spatial_network; 
            maximal_wavevector_int=maximal_wavevector_int),
    save_result = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing)

    # initialize structure factor array
    structure_factor_array = Array{Float64}(undef, 
        size(wavevector_array_positive_z)[1:3]...)

    # get structure factor for all wavevectors
    for i in 1:size(wavevector_array_positive_z)[1], 
        j in 1:size(wavevector_array_positive_z)[2], 
        k in 1:size(wavevector_array_positive_z)[3]

        structure_factor_array[i,j,k] = get_structure_factor(
            wavevector_array_positive_z[i,j,k,:], spatial_network)

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
        GU.save_dict_to_h5(copy(structure_factor_dict), save_path*"_structure_factor_array.h5")
    end
                                            
    return structure_factor_dict
end


"""
From the dictionary containing the 3d structure factor, get a dictionary that,
for each index vector squared contains all values of the structure factor. The 
index vector squared is proportional to the wavenumber squared. The returned 
dictionary can be used to calculate the angle averaged structure factor or an
anisotropy metric
"""
function get_structure_factor_by_index_vec_squared_dict(
    structure_factor_dict::Dict;)

    structure_factor_by_index_vec_squared_dict = Dict{Int64, Vector{Float64}}()

    # determine origin vector of cartesian coordinates
    index_vector_origin = Int.((size(
        structure_factor_dict["structure_factor_array"]) .+ 1 ) ./ 2)

    # get lattice constant of reciprocal lattice
    reciprocal_lattice_constant = LinearAlgebra.norm(
        structure_factor_dict["wavevector_array"][(index_vector_origin 
            .+ [1,0,0])...,:])

    # loop through Cartesian Indices
    for i in CartesianIndices(structure_factor_dict["structure_factor_array"])

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
                structure_factor_dict["structure_factor_array"][i])

        else

            # add key and structure factor to dictionary
            structure_factor_by_index_vec_squared_dict[index_vector_length_squared]=(
                [structure_factor_dict["structure_factor_array"][i]])

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
    gaussian_filter::Bool = true,
    gaussian_filter_sigma_x::Float64 = 2*pi/25, 
    gaussian_filter_filtered_data_x_step_length::Float64 = 2*pi/25,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing)

    # get the dictionary that, for each index vector squared contains all 
    # values of the structure factor. The index vector squared is proportional 
    # to the square of the wavenumber
    structure_factor_by_index_vec_squared_dict, reciprocal_lattice_constant = (
        get_structure_factor_by_index_vec_squared_dict(
        structure_factor_dict))
    
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
        "unfiltered_structure_factor_vec" => unfiltered_structure_factor_vec)

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
        structure_factor_angle_averaged_dict["structure_factor_vec"] = (
            filtered_data_y)
        structure_factor_angle_averaged_dict["gaussian_filter_sigma_x"] = (
            gaussian_filter_sigma_x)
    end

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_angle_averaged_dict["label"] = label
    end

    if save_result
        GU.save_dict_to_h5(copy(structure_factor_angle_averaged_dict),
            save_path*"_structure_factor_angle_averaged.h5")
    end
                                            
    return structure_factor_angle_averaged_dict
end


"""
Calculate pair correlation function, total correlation function and
cumulative coordination number as defined in equations 3-5 of
10.1016/j.physrep.2018.03.001
"""
function get_correlation_functions(
    spatial_network::MetaGraphsNext.MetaGraph;
    distance_histogram_bin_width::Float64 = 0.02,
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

    current_index = 1

    for i in 1:spatial_network[]["nr_vertices"]-1
        for j in i+1:spatial_network[]["nr_vertices"]
            distance_vec[current_index] = LinearAlgebra.norm(
                NG.get_distance_vector_pbc(
                    vertex_position_mat[:, i], vertex_position_mat[:, j], 
                    spatial_network[]["supercell_edge_length"]))

            current_index += 1
        end
    end

    # get histogram of distance vector
    distance_histogram = StatsBase.fit(
        StatsBase.Histogram, distance_vec, 
        0.0:distance_histogram_bin_width:spatial_network[][
            "supercell_edge_length"]/2, 
        closed=:left)

    vertex_distance_vec = collect(distance_histogram.edges[1][2:end] 
        .- distance_histogram_bin_width/2)
    vertex_nr_vec = distance_histogram.weights

    # get_vertex density
    vertex_density = spatial_network[]["nr_vertices"] / spatial_network[][
        "supercell_edge_length"]^3

    # calculate pair correlaion function
    # (eq. 3 in 10.1016/j.physrep.2018.03.001)
    pair_correlation_fct_vec = ((1/(4*pi*vertex_density 
            * spatial_network[]["nr_vertices"]) ) 
        .* vertex_nr_vec ./ vertex_distance_vec.^2) .* (
            2/distance_histogram_bin_width)

    # calculate total correlation function
    # (eq. 4 in 10.1016/j.physrep.2018.03.001)
    total_correlation_fct_vec = pair_correlation_fct_vec .- 1

    # calculate cumulative coordination number 
    # (eq. 5 in 10.1016/j.physrep.2018.03.001)
    cumulative_coord_nr_vec = (4*pi*vertex_density*distance_histogram_bin_width
        * cumsum(vertex_distance_vec.^2 .* pair_correlation_fct_vec ))

    # create dict to save
    correlation_functions_dict = Dict{String, Any}(
        "vertex_distance_vec" => vertex_distance_vec, 
        "cumulative_coord_nr_vec" => cumulative_coord_nr_vec,
        "pair_correlation_fct_vec" => pair_correlation_fct_vec,
        "total_correlation_fct_vec" => total_correlation_fct_vec)

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
Define a homogeneity metric as the average sphere radius within which there is 
the given number of vertices. If the homogeneity metric is close to 0, the
network is clustered
"""
function get_vertex_homogeneity_metric(correlation_functions_dict::Dict;
    nr_vertices_within_sphere::Int = 10)

    # get distance where cumulative coordination nr equals the given number of 
    # vertices. To account for the fact that the central vertex is not
    # contained in the cumulative coordination number, 1 has to be subtracted
    vertex_homogeneity_metric = (
        correlation_functions_dict["vertex_distance_vec"][
            findfirst(x -> x > nr_vertices_within_sphere-1, 
            correlation_functions_dict["cumulative_coord_nr_vec"])])

    return vertex_homogeneity_metric
end


"""
Get the second moment of the pore size distribution which, according to 
10.1103/PhysRevE.104.014127 is an estimate of the critical pore radius squared
(although they use another definition of the pore size distribution)
"""
function get_pore_size_distribution_second_moment(
    pore_size_distribution_dict::Dict)

    pore_size_step_length = (pore_size_distribution_dict["pore_size_vec"][2] 
        - pore_size_distribution_dict["pore_size_vec"][1])

    # get second moment of pore size distribution
    pore_size_distribution_second_moment = pore_size_step_length * sum(
        pore_size_distribution_dict["pore_size_vec"].^2 
        .* pore_size_distribution_dict["pore_size_distribution"])

    return pore_size_distribution_second_moment
end


"""
Define an anisotropy metric based on the structure factor or the 
spectral density that ranges between 0 and 1 (high anisotropy). The metric is
based on the coefficient of variation (std over mean) of the structure factor
as a function of direction. In reality, values around 0.45 represent low
anisotropy and values around 0.55 represent high anisotropy.
"""
function get_anisotropy_metric_from_structure_factor(
    structure_factor_dict::Dict;
    maximal_length_to_check = 3.0,
    nr_closest_wavenumbers = 3,
    normalization_parameter = 1.0)

    # set the wavenumbers where structure factor will be checked
    wavenumbers_to_check_vec = (2*pi) ./ collect(
        0.5:0.5:maximal_length_to_check+0.01)

    # get the dictionary that, for each index vector squared contains all 
    # values of the structure factor. The index vector squared is proportional 
    # to the square of the wavenumber
    structure_factor_by_index_vec_squared_dict, reciprocal_lattice_constant = (
        get_structure_factor_by_index_vec_squared_dict(
        structure_factor_dict))

    # first, create a vector of all index vectors squared
    index_vec_squared_vec = collect(keys(
        structure_factor_by_index_vec_squared_dict))

    # get the three index vectors squared that are closest to the wavenumbers
    # to check
    index_vec_squared_to_check_arr = Array{Int64}(undef, 
    nr_closest_wavenumbers, 
        length(wavenumbers_to_check_vec))
    for i in eachindex(wavenumbers_to_check_vec)
        # get the three index vectors squared that are closest to the wavenumber
        # to check
        index_vec_squared_to_check_arr[:, i] = index_vec_squared_vec[sortperm(
            abs.(index_vec_squared_vec .- (wavenumbers_to_check_vec[i]/
                reciprocal_lattice_constant)^2))][1:nr_closest_wavenumbers]
    end

    # for each index vector squared, calculate a normalized coefficient of
    # variation (std over mean)
    anisotropy_metric_vec = Vector{Float64}(undef, 
        length(wavenumbers_to_check_vec))

    for i in eachindex(wavenumbers_to_check_vec)
        # get the structure factor for the three index vectors squared
        structure_factor_vec =
            structure_factor_by_index_vec_squared_dict[
                index_vec_squared_to_check_arr[1, i]] 
        for j in 2:nr_closest_wavenumbers
            structure_factor_vec = vcat(
                structure_factor_vec, 
                structure_factor_by_index_vec_squared_dict[
                    index_vec_squared_to_check_arr[j, i]])
        end

        # calculate the coefficient of variation
        coefficient_of_variation = (Statistics.std(structure_factor_vec)
            /(Statistics.mean(structure_factor_vec)))

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
    window_radius)

    # check if system is 3d
    if spatial_network[]["nr_dimensions"] != 3
        @error "Calculation of number variance is only implemented for 3d
            systems"
    end

    # calculate vertex density
    vertex_density = (spatial_network[]["nr_vertices"] 
        / spatial_network[]["supercell_edge_length"]^spatial_network[][
            "nr_dimensions"])

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
Measure local nr variance as a function of window radius according from
structure factor according to equations 47 and d58 in 
10.1016/j.physrep.2018.03.001
"""
function get_local_nr_variance_by_window_radius_vec(
    spatial_network,
    structure_factor_angle_averaged_dict::Dict;
    window_radius_step_length::Float64 = 0.1,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\random_networks\sample_name",
    label = nothing)

    # get vector of sphere radii
    sphere_radius_vec = collect(window_radius_step_length:
    window_radius_step_length
    :spatial_network[]["supercell_edge_length"]/4
        )

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
Get hyperuniformity metric which is the structure factor at zero momentum
normalized by the height of the first peak in the structure factor as defined
in equation 251 in 10.1016/j.physrep.2018.03.001
"""
function get_hyperuniformity_metric(structure_factor_dict::Dict)

    # locate peaks of structure factor
    pks, vals = Peaks.findmaxima(structure_factor_dict["structure_factor_vec"])

    # cut structure factor data at momentum just above first peak
    structure_factor_cut_vec = structure_factor_dict[
        "structure_factor_vec"][1:pks[2]-1]
    wavenumber_cut_vec = structure_factor_dict["wavenumber_vec"][1:pks[2]-1]

    # set the order of the fitted polynomial
    polynomial_order = 5

    # fit polynomial of given order to cut data
    polynomial_fit = Polynomials.fit(wavenumber_cut_vec, 
                                    structure_factor_cut_vec,
                                    polynomial_order)

    # get extrapolated structure factor at zero momentum
    structure_factor_zero_momentum = polynomial_fit(0)

    # get first derivative of polynomial
    polynomial_derivative = Polynomials.derivative(polynomial_fit)

    # in case the structure factor is provided with uncertainty, obtain the
    # values
    polynomial_derivative_values = Polynomials.Polynomial(Measurements.value.( 
        collect(polynomial_derivative) ))
    
    # get critical momenta which is roots of first derivative of polynomial
    critical_momenta = Polynomials.roots(polynomial_derivative_values)

    # get real critical momenta
    critical_momenta_real = real.(
        critical_momenta[imag.(critical_momenta) .== 0])

    # get fitted structure factor at highest peak
    structure_factor_first_peak = maximum( 
        polynomial_fit.(critical_momenta_real) )

    # get hyperuniformity metric
    hyperuniformity_metric = (structure_factor_zero_momentum
        /structure_factor_first_peak)

    return [hyperuniformity_metric, polynomial_fit]
end


