"""
these functions can be used to characterize networks
by means of order parameters measuring correlations
"""

"""
Measure structure factor for a given wavenumber,
averaged over angles according to Barlett's isotropic estimator
as described in equation 40 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_isotrope(graph_dict::Dict, wavenumber::Real)

    #check if structure is 3d
    if graph_dict["nr_dimensions"] !== 3
        @error "Structure factor calculation is, so far, only implemented for 3d."
    end

    #initialize double sum
    double_sum = 0

    #perform sum over all combinations of vertices
    for vertex_1 in MetaGraphsNext.labels(graph_dict["spatial_network"])

        #get first vertex position
        vertex_1_pos = graph_dict["spatial_network"][vertex_1]["position"]

        for vertex_2 in vertex_1+1:graph_dict["nr_vertices"]

            #get second vertex position
            vertex_2_pos = graph_dict["spatial_network"][vertex_2]["position"]

            #distance between vertices
            vertex_distance = LinearAlgebra.norm(vertex_2_pos .- vertex_1_pos)

            #calculate structure factor contribution of current vertex combination
            double_sum += sin(wavenumber*vertex_distance) / (
                                    wavenumber*vertex_distance)

        end
    end

    #calculate structure factor
    structure_factor = 1 + 2/graph_dict["nr_vertices"] * double_sum

    return structure_factor
end


"""
Get vector of wavenumbers, for which the structure factor is calculated
"""
function get_wavenumber_vec(graph_dict;
    sampling_distance_step_length::Real = 0.1,
    maximal_sampling_distance = graph_dict["supercell_edge_length"])

    #determine virtual nr of sampling distances
    #(in reality I don't sample in direct space anywhere)
    nr_sampling_distances = floor(maximal_sampling_distance
                                /(sampling_distance_step_length))

    #get nr of wavenumbers which is half the nr of sampling distances
    nr_wavenumbers = Int(floor(nr_sampling_distances/2))

    #get fundamental wavenumber
    fundamental_wavenumber = 2*pi/maximal_sampling_distance

    #get vector of wavenumbers
    wavenumber_vec = collect(1:nr_wavenumbers) * fundamental_wavenumber

    return wavenumber_vec
end


"""
Measure structure factor as a function of wavenumber
averaged over angles according to Barlett's isotropic estimator
as described in equation 40 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict::Dict,
    sampling_distance_step_length::Real = 0.1,
    maximal_sampling_distance = graph_dict["supercell_edge_length"],
    save_result = false,
    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\sample_name",
    label = nothing)

    #get vector of wavenumbers
    wavenumber_vec = get_wavenumber_vec(graph_dict; 
        sampling_distance_step_length = sampling_distance_step_length,
        maximal_sampling_distance = maximal_sampling_distance)

    #initialize structure factor vector
    structure_factor_vec = Vector{Float64}(undef, length(wavenumber_vec))

    #get vector of structure factor as a function of wavenumber
    for i in eachindex(wavenumber_vec)
        structure_factor_vec[i] = get_structure_factor_isotrope(graph_dict, wavenumber_vec[i])

    end

    #create dictionary for current plot
    structure_factor_dict = Dict("wavenumber_vec" => wavenumber_vec,
                            "structure_factor_vec" => structure_factor_vec,
                            "sampling_distance_step_length" => 
                            sampling_distance_step_length,
                            "maximal_sampling_distance" => maximal_sampling_distance,
                            "label" => label )

    #save results if desired
    if save_result
        GU.save_dict_to_h5(copy(structure_factor_dict);
                        save_path=save_path*"_structure_factor_isotrope.h5")

    end

    return [wavenumber_vec, structure_factor_vec]
end


"""
Measure structure factor for a given wavevector by determining
the so called scattering intensity as described in equation 26 
of 10.1007/s11222-023-10219-1
"""
function get_structure_factor(graph_dict::Dict, wavevector::Vector{Real})

    #to be implemented

end


"""
Measure structure factor as a function of wavenumber
by determining the so called scattering intensity
as described in equation 26 of 10.1007/s11222-023-10219-1
"""
function get_structure_factor_by_wavevector_array(graph_dict::Dict)

    #to be implemented

end