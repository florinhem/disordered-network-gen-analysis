"""
these functions can be used to characterize networks
by means of order metrics measuring correlations
"""

"""
Measure structure factor for a given wavenumber,
averaged over angles according to Barlett's isotropic estimator
as described in equation 40 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_bartlett_isotrope(graph_dict::Dict, wavenumber)

    # check if structure is 3d
    if graph_dict["nr_dimensions"] !== 3
        @error "Structure factor calculation is, so far, only implemented for 3d."
    end

    # initialize double sum
    double_sum = 0.0

    # perform sum over all combinations of vertices
    for vertex_1 in MetaGraphsNext.labels(graph_dict["spatial_network"])

        # get first vertex position
        vertex_1_pos = graph_dict["spatial_network"][vertex_1]["position"]

        for vertex_2 in vertex_1+1:graph_dict["nr_vertices"]

            # get second vertex position
            vertex_2_pos = graph_dict["spatial_network"][vertex_2]["position"]

            # distance between vertices
            vertex_distance = LinearAlgebra.norm(vertex_2_pos .- vertex_1_pos)

            # calculate structure factor contribution of current vertex combination
            double_sum += sin(wavenumber*vertex_distance) / (
                                    wavenumber*vertex_distance)

        end
    end

    # calculate structure factor
    structure_factor = 1 + 2/graph_dict["nr_vertices"] * double_sum

    return structure_factor
end


"""
Get vector of wavenumbers, for which the structure factor is calculated
"""
function get_wavenumber_vec(graph_dict::Dict;
    sampling_distance_step_length = 0.05,
    maximal_sampling_distance = graph_dict["supercell_edge_length"]/sqrt(3))

    # determine virtual nr of sampling distances
    # (in reality I don't sample in direct space anywhere)
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
Measure structure factor as a function of wavenumber
averaged over angles according to Barlett's isotropic estimator
as described in equation 26 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_bartlett_isotrope_by_wavenumber_vec(
    graph_dict::Dict;
    sampling_distance_step_length = 0.1,
    maximal_sampling_distance = graph_dict["supercell_edge_length"]/2,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\random_networks\sample_name",
    print_progress::Bool = false,
    label = nothing)

    # get vector of wavenumbers
    wavenumber_vec = get_wavenumber_vec(graph_dict; 
        sampling_distance_step_length = sampling_distance_step_length,
        maximal_sampling_distance = maximal_sampling_distance)

    # initialize structure factor vector
    structure_factor_bartlett_vec = Vector{Float64}(undef, length(wavenumber_vec))

    # get vector of structure factor as a function of wavenumber
    for i in eachindex(wavenumber_vec)
        structure_factor_bartlett_vec[i] = get_structure_factor_bartlett_isotrope(
                                            graph_dict, wavenumber_vec[i])

        # print progress
        if print_progress
            println("Progress: ", i/length(wavenumber_vec)*100, "%")
        end

    end

    # create dictionary for current plot
    structure_factor_bartlett_dict = Dict("wavenumber_vec" => wavenumber_vec,
                            "structure_factor_vec" => structure_factor_bartlett_vec,
                            "sampling_distance_step_length" => 
                            sampling_distance_step_length,
                            "maximal_sampling_distance" => maximal_sampling_distance )

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_bartlett_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(structure_factor_bartlett_dict);
                        save_path=save_path*"_structure_factor_bartlett_isotrope.h5")

    end

    return structure_factor_bartlett_dict
end


"""
Measure structure factor as a function of wavenumber
averaged over angles according to scattering intensity estimator
as described in equation 40 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_isotrope(graph_dict::Dict, wavenumber;
    nr_wavevector_samples::Int = 10000)

    # check if structure is 3d
    if graph_dict["nr_dimensions"] !== 3
        @error "Structure factor calculation is, so far, only implemented for 3d."
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
        for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

            # get vertex position
            vertex_pos = graph_dict["spatial_network"][vertex]["position"]
        
            # calculate structure factor contribution of current vertex and wavevector
            scattering_field_sum += exp(-im*LinearAlgebra.dot(wavevector, vertex_pos))

        end

        # calculate structure factor
        structure_factor = 1/graph_dict["nr_vertices"] * abs2(scattering_field_sum)

        # add structure factor to sum
        structure_factor += structure_factor/nr_wavevector_samples

    end

    return structure_factor
end


"""
Measure structure factor for a given wavenumber,
averaged over angles according to the scattering intensity estimator
as described in equation 26 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict::Dict;
    sampling_distance_step_length = 0.1,
    maximal_sampling_distance = graph_dict["supercell_edge_length"]/2,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\random_networks\sample_name",
    print_progress::Bool = false,
    label = nothing)

    # get vector of wavenumbers
    wavenumber_vec = get_wavenumber_vec(graph_dict; 
        sampling_distance_step_length = sampling_distance_step_length,
        maximal_sampling_distance = maximal_sampling_distance)

    # initialize structure factor vector
    structure_factor_vec = Vector{Float64}(undef, length(wavenumber_vec))

    # get vector of structure factor as a function of wavenumber
    for i in eachindex(wavenumber_vec)
        structure_factor_vec[i] = get_structure_factor_isotrope(
                                            graph_dict, wavenumber_vec[i])

        # print progress
        if print_progress
            println("Progress: ", i/length(wavenumber_vec)*100, "%")
        end

    end

    # create dictionary for current plot
    structure_factor_dict = Dict("wavenumber_vec" => wavenumber_vec,
                            "structure_factor_vec" => structure_factor_vec,
                            "sampling_distance_step_length" => 
                            sampling_distance_step_length,
                            "maximal_sampling_distance" => maximal_sampling_distance,
                            "label" => label )

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(structure_factor_dict);
                        save_path=save_path*"_structure_factor_isotrope.h5")

    end

    return structure_factor_dict
end


"""
Get array of wavevectors along positive z direction for which the structure factor is calculated
when isotropy is not assumed.
Each wavevector coordinate is spanned by k_j = 2*pi*n_j/L_j where n_j is an integer within
-maximal_wavevector_int*supercell_edge_length <= n_j <= maximal_wavevector_int*supercell_edge_length
for j = 1,2 and
0 <= n_j <= maximal_wavevector_int*supercell_edge_length for j = 3
"""
function get_wavevector_array_positive_z(graph_dict::Dict;
    maximal_wavevector_int::Int64 = 5)

    # get nr of wavevectors per positive direction
    nr_wavevectors_per_pos_direction = Int(ceil(maximal_wavevector_int*graph_dict["supercell_edge_length"]))

    # generate empty array for wavevectors positions
    wavevector_array_positive_z = Array{Float64}(undef, 
                                        2*nr_wavevectors_per_pos_direction + 1,
                                        2*nr_wavevectors_per_pos_direction + 1,
                                        nr_wavevectors_per_pos_direction + 1,
                                        graph_dict["nr_dimensions"])

    # loop through Cartesian Indices
    for i in CartesianIndices(wavevector_array_positive_z)

        # fill wavevector array depending on dimension of current index
        wavevector_array_positive_z[i] =  2*pi/graph_dict["supercell_edge_length"]*(
            (i[1] - nr_wavevectors_per_pos_direction - 1)* ==(i[4], 1) 
            + (i[2] - nr_wavevectors_per_pos_direction - 1)* ==(i[4], 2) 
            + (i[3] - 1)* ==(i[4], 3)
           ) 

    end

    return wavevector_array_positive_z
end


"""
Measure structure factor as a function of wavevector
using the scattering intensity estimator
as described in equation 24 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor(wavevector::Vector{Float64},
    graph_dict::Dict)

    # initialize the sum of the scattering field
    scattering_field_sum = 0.0 + 0.0*im
    
    # perform sum over all vertices
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        # get vertex position
        vertex_pos = graph_dict["spatial_network"][vertex]["position"]
    
        # calculate structure factor contribution of current vertex and wavevector
        scattering_field_sum += exp(-im*LinearAlgebra.dot(wavevector, vertex_pos))

    end

    # calculate structure factor
    structure_factor = 1/graph_dict["nr_vertices"] * abs2(scattering_field_sum)

    return structure_factor
end


"""
Measure structure factor for an array of wavevectors
using the scattering intensity estimator
as described in equation 24 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_by_wavevector_array(graph_dict::Dict;
    maximal_wavevector_int::Int64 = 5,
    wavevector_array_positive_z::Array{Float64} = 
        get_wavevector_array_positive_z(graph_dict; maximal_wavevector_int=maximal_wavevector_int),
    save_result = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing)

    # initialize structure factor array
    structure_factor_array = Array{Float64}(undef, size(wavevector_array_positive_z)[1:3]...)

    # get structure factor for all wavevectors
    for i in 1:size(wavevector_array_positive_z)[1], j in 1:size(wavevector_array_positive_z)[2], 
        k in 1:size(wavevector_array_positive_z)[3]

        structure_factor_array[i,j,k] = get_structure_factor(wavevector_array_positive_z[i,j,k,:], graph_dict)

    end

    # extend structure factor to negative z direction
    structure_factor_array = cat(dims=3, 
            structure_factor_array[end:-1:1,end:-1:1,end:-1:2], structure_factor_array)

    # extend wavevector array to negative z direction
    wavevector_array_negative_z = reverse(wavevector_array_positive_z[:,:,2:end,:], dims=3)
    wavevector_array_negative_z[:,:,:,3] .*= (-1)
    wavevector_array = cat(dims=3, 
        wavevector_array_negative_z, wavevector_array_positive_z)

    # get vector of vector of wavenumbers along all directions
    wavenumber_vec_vec = [wavevector_array[:,1,1,1], 
                        wavevector_array[1,:,1,2], 
                        wavevector_array[1,1,:,3]]

    # create dict to save
    structure_factor_dict = Dict("wavevector_array" => wavevector_array, 
                                "wavenumber_vec_vec" => wavenumber_vec_vec,
                                "structure_factor_array" => structure_factor_array)

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(structure_factor_dict);
                        save_path=save_path*"_structure_factor_array.h5")

    end
                                            
    return structure_factor_dict
end


"""
Calculate angle averaged structure factor from 3d array of structure factor
"""
function get_structure_factor_angle_averaged(structure_factor_dict::Dict;
    gaussian_filter::Bool = true,
    gaussian_filter_sigma_x::Float64 = 2*pi/10, 
    gaussian_filter_filtered_data_x_step_length::Float64 = 2*pi/10,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing)

    # create a dictionary for angle averaged structure factor
    # with the length squared of the index vector as the key,
    # because it is proportional to the square of the wavenumber
    structure_factor_angle_averaged_dict = Dict{Int64, Vector{Float64}}()

    # determine origin vector of cartesian coordinates
    index_vector_origin = Int.((size(structure_factor_dict["structure_factor_array"]) .+ 1 ) ./ 2)

    # loop through Cartesian Indices
    for i in CartesianIndices(structure_factor_dict["structure_factor_array"])

        # determine actual index vector
        index_vector = i.I .- index_vector_origin

        # calculate length squared of index vector
        index_vector_length_squared = sum(index_vector.^2)

        # check if key exists in dictionary
        if index_vector_length_squared in keys(structure_factor_angle_averaged_dict)

            # add structure factor to dictionary
            push!(structure_factor_angle_averaged_dict[index_vector_length_squared], 
                structure_factor_dict["structure_factor_array"][i])

        else

            # add key and structure factor to dictionary
            structure_factor_angle_averaged_dict[index_vector_length_squared] = 
                [structure_factor_dict["structure_factor_array"][i]]

        end

    end
    
    # initialize wavenumber vector and structure factor vector
    wavenumber_vec = Vector{Float64}()
    structure_factor_angle_averaged_vec = Vector{Measurements.Measurement{Float64}}()

    # get lattice constant of reciprocal lattice
    reciprocal_lattice_constant = LinearAlgebra.norm(
        structure_factor_dict["wavevector_array"][(index_vector_origin .+ [1,0,0])...,:])

    # get vector of wavenumbers and angle averaged structure factor including
    # their uncertainty
    for key in keys(structure_factor_angle_averaged_dict)

        # get wavenumber from wavevector
        push!(wavenumber_vec, reciprocal_lattice_constant*sqrt(key))

        # get angle averaged structure factor
        push!(structure_factor_angle_averaged_vec, 
            Measurements.measurement(Statistics.mean(structure_factor_angle_averaged_dict[key]),
            Statistics.std(structure_factor_angle_averaged_dict[key])))

    end

    # sort wavenumber vector and structure factor vector
    structure_factor_angle_averaged_vec = structure_factor_angle_averaged_vec[sortperm(wavenumber_vec)]
    sort!(wavenumber_vec)

    # create dict to save
    structure_factor_angle_averaged_dict = Dict()

    structure_factor_angle_averaged_dict["wavenumber_vec"] = wavenumber_vec
    structure_factor_angle_averaged_dict["structure_factor_angle_averaged_vec"] = structure_factor_angle_averaged_vec

    # apply gaussian filter if desired
    if gaussian_filter
        filtered_data_x, filtered_data_y = GU.gaussian_filter_1d(wavenumber_vec[2:end], 
            structure_factor_angle_averaged_vec[2:end]; 
            sigma_x=gaussian_filter_sigma_x, 
            filtered_data_x_step_length=gaussian_filter_filtered_data_x_step_length)

        structure_factor_angle_averaged_dict["filtered_wavenumber_vec"] = filtered_data_x
        structure_factor_angle_averaged_dict["filtered_structure_factor_angle_averaged_vec"] = filtered_data_y
        structure_factor_angle_averaged_dict["gaussian_filter_sigma_x"] = gaussian_filter_sigma_x
    end

    # add label to dictionary if label is not nothing
    if label !== nothing
        structure_factor_angle_averaged_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(structure_factor_angle_averaged_dict);
                        save_path=save_path*"_structure_factor_angle_averaged.h5")

    end
                                            
    return structure_factor_angle_averaged_dict
end


"""
Get local number variance for a spherical window of given radius from the
structure factor according to eq 58 in 10.1016/j.physrep.2018.03.001
"""
function get_local_nr_variance(graph_dict::Dict,
    structure_factor_dict::Dict,
    window_radius)

    # check if system is 3d
    if graph_dict["nr_dimensions"] != 3
        @error "Calculation of number variance is only implemented for 3d systems"
    end

    # calculate vertex density
    vertex_density = (graph_dict["nr_vertices"] 
        / graph_dict["supercell_edge_length"]^graph_dict["nr_dimensions"])

    # calculate wavenumber sampling step length
    wavenumber_step_length = (structure_factor_dict["wavenumber_vec"][2] 
                                - structure_factor_dict["wavenumber_vec"][1])

    # calculate local nr variance
    local_nr_variance = (vertex_density * 32 * pi^2 * window_radius^2
        * wavenumber_step_length * sum(
            structure_factor_dict["structure_factor_vec"] 
            .* sin.( structure_factor_dict["wavenumber_vec"] .* window_radius) .^2
            ./ structure_factor_dict["wavenumber_vec"].^2
        )    
    )
    
    return local_nr_variance
end


"""
Measure local nr variance as a function of window radius according
from structure factor
"""
function get_local_nr_variance_by_window_radius_vec(
    graph_dict::Dict;
    structure_factor_dict::Dict = get_structure_factor_isotrope_by_wavenumber_vec(graph_dict),
    window_radius_step_length = 0.2,
    maximal_window_radius = graph_dict["supercell_edge_length"]/2,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\random_networks\sample_name",
    label = nothing)

    # get vector of winow radii
    window_radius_vec = collect(
        window_radius_step_length:window_radius_step_length:maximal_window_radius)

    # initialize local nr variance vector
    local_nr_variance_vec = Vector{Float64}(undef, length(window_radius_vec))

    # get local nr variance vector as a window radius
    for i in eachindex(window_radius_vec)
        local_nr_variance_vec[i] = get_local_nr_variance(graph_dict,
        structure_factor_dict,
        window_radius_vec[i])

    end

    # create dictionary for current plot
    local_nr_variance_dict = Dict("window_radius_vec" => window_radius_vec,
                            "local_nr_variance_vec" => local_nr_variance_vec )

    if label !== nothing
        local_nr_variance_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(local_nr_variance_dict);
                        save_path=save_path*"_local_nr_variance.h5")

    end

    return local_nr_variance_dict
end


"""
Get hyperuniformity metric which is the structure factor
at zero momentum normalized by the height of the first peak in the structure factor
as defined in equation 251 in 10.1016/j.physrep.2018.03.001
"""
function get_hyperuniformity_metric(structure_factor_dict::Dict)

    # locate first peak of structure factor
    pks, vals = Peaks.findmaxima(structure_factor_dict["structure_factor_vec"])

    # cut structure factor data at momentum just above first peak
    structure_factor_cut_vec = structure_factor_dict["structure_factor_vec"][1:pks[2]+1]
    wavenumber_cut_vec = structure_factor_dict["wavenumber_vec"][1:pks[2]+1]

    # set the order of the fitted polynomial
    polynomial_order = 5

    # fit polynomial of given order to cut data
    polynomial_fit = Polynomials.fit(wavenumber_cut_vec, 
                                    structure_factor_cut_vec,
                                    polynomial_order)

    # get extrapolated structure factor at zero momentum
    structure_factor_zero_momentum = polynomial_fit(0)

    # get critical momenta which is roots of first derivative of polynomial
    polynomial_derivative = Polynomials.derivative(polynomial_fit)
    critical_momenta = Polynomials.roots(polynomial_derivative)

    # get real critical momenta
    critical_momenta_real = real.(critical_momenta[imag.(critical_momenta) .== 0])

    # get fitted structure factor at highest peak
    structure_factor_first_peak = maximum( polynomial_fit.(critical_momenta_real) )

    # get hyperuniformity metric
    hyperuniformity_metric = structure_factor_zero_momentum/structure_factor_first_peak

    return [hyperuniformity_metric, polynomial_fit]
end


