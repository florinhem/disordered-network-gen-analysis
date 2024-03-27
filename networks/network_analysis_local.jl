"""
these functions can be used to characterize networks
by means of local order metrics
"""

"""
Measure the standard deviation of bond lengths
"""
function get_bond_length_std(graph_dict::Dict)
    
    # get nr of bonds
    nr_bonds = Int(graph_dict["nr_vertices"] 
                    * graph_dict["coordination_nr"] /2)

    # get vector of bond lengths
    bond_length_vec = Vector{Float64}(undef, nr_bonds)
    bond_count = 1

    for bond in MetaGraphsNext.edge_labels(
                        graph_dict["spatial_network"])

        bond_length_vec[bond_count] = sqrt(
            graph_dict["spatial_network"][bond...]["distance_squared"]
        )

        bond_count += 1
    end

    # determine standard deviation
    bond_length_std = Statistics.std(bond_length_vec)
    
    # get histogram of bond lengths
    # histogram = StatsBase.fit(Histogram, bond_length_vec)
    # bond_length_bin_edges = histogram.edges
    # bond_length_bin_weights = histogram.weights

    return [bond_length_std, bond_length_vec]
end


"""
Measure the standard deviation of bond angles
"""
function get_bond_angle_std(graph_dict::Dict)
    
    # get nr of angles
    nr_angles = (graph_dict["nr_vertices"] 
                    * sum(1:graph_dict["coordination_nr"]-1))

    # initialize vector of bond angles
    bond_angle_vec = Vector{Float64}(undef, nr_angles)
    angle_count = 1

    # loop through all vertices
    for vertex in MetaGraphsNext.labels(
                        graph_dict["spatial_network"])

        # get iterator of bond combinations
        bond_combinations_iter = Combinatorics.combinations(
            collect(MetaGraphsNext.neighbor_labels(
                graph_dict["spatial_network"], vertex)), 2)

        # loop through bond combinations
        for bond_combination in bond_combinations_iter

            # get scalar product of vectors representing the current
            # bond combination
            scalar_product = (sign(bond_combination[1] - vertex)
                * sign(bond_combination[2] - vertex)
                * LinearAlgebra.dot(
graph_dict["spatial_network"][vertex, bond_combination[1]]["vector"], 
graph_dict["spatial_network"][vertex, bond_combination[2]]["vector"] )
              )

            # calculate bond angle
            bond_angle = acos( scalar_product/ sqrt(
graph_dict["spatial_network"][vertex, bond_combination[1]]["distance_squared"]
*graph_dict["spatial_network"][vertex, bond_combination[2]]["distance_squared"]) )

            # save bond_angle
            bond_angle_vec[angle_count] = bond_angle

            angle_count += 1

        end
    end

    # determine standard deviation
    bond_angle_std = Statistics.std(bond_angle_vec)

    return [bond_angle_std, bond_angle_vec]
end


"""
Measure the standard deviation of dihedral angles.
This function might need to be revisited when considering
other coordination numbers than 4 because the dihedral angle
might have a peak around 0 which is not considered here
"""
function get_dihedral_angle_std(graph_dict::Dict)

    # initialize vector of diehedral angles
    dihedral_angle_vec = Vector{Float64}()

    # loop through all bonds
    for bond in MetaGraphsNext.edge_labels(
        graph_dict["spatial_network"])

        # get vector along bond
        bond_vec = graph_dict["spatial_network"][bond...]["vector"]

        # loop through all neighbors of one vertex 
        for first_neighbor in setdiff( MetaGraphsNext.neighbor_labels(
                            graph_dict["spatial_network"], bond[1]), bond[2])

            # get vector from first neighbor to first bond vertex
            first_neighbor_to_bond_vertex_vec = ( sign(bond[1] - first_neighbor)
                *graph_dict["spatial_network"][first_neighbor, bond[1]]["vector"]
            )

            # loop through all neighbors of other vertex
            for second_neighbor in setdiff( MetaGraphsNext.neighbor_labels(
                graph_dict["spatial_network"], bond[2]), bond[1])

                # get vector from second bond vertex to second neighbor
                bond_vertex_to_second_neighbor_vec = ( sign(second_neighbor - bond[2])
                    *graph_dict["spatial_network"][bond[2], second_neighbor]["vector"]
                )

                # calculate dihedral angle according to the equation given in
                # https://en.wikipedia.org/wiki/Dihedral_angle# In_polymer_physics
                dihedral_angle = atan( LinearAlgebra.norm(bond_vec)
                        *LinearAlgebra.dot(first_neighbor_to_bond_vertex_vec,
                            LinearAlgebra.cross(bond_vec,
                                        bond_vertex_to_second_neighbor_vec )
                    ),
                    LinearAlgebra.dot(
                        LinearAlgebra.cross(first_neighbor_to_bond_vertex_vec,
                            bond_vec ),
                        LinearAlgebra.cross(bond_vec,
                              bond_vertex_to_second_neighbor_vec )
                    ))

                # save dihedral angle
                push!(dihedral_angle_vec, dihedral_angle)
            end
        end
    end

    # save one peak of dihedral angle distribution
    lower_limit = 0
    upper_limit =  2 * pi / (graph_dict["coordination_nr"] - 1)

    dihedral_angle_one_peak_vec = dihedral_angle_vec[
        (dihedral_angle_vec .> lower_limit) .& (dihedral_angle_vec .< upper_limit)]

    # determine standard deviation
    dihedral_angle_std = Statistics.std(dihedral_angle_one_peak_vec)


    return [dihedral_angle_std, dihedral_angle_vec]
end


"""
Get dictionary of q_lm averaged over bonds to neighbors for a single vertex
(also called Steinhardt local bond order parameters) 
as in equation 2.1 of 10.1103/PhysRevB.28.784 which is not rotationally
invariant
"""
function get_q_lm_averaged_bonds_to_neighbors_dict_single_vertex(graph_dict::Dict,
    central_vertex::Int64,
    l_max::Int64)

    # initialize vector for spherical harmonics of all neighbors
    y_spherical_harmonic_arr_vec = Vector{SphericalHarmonics.SHVector{
        ComplexF64, 
        Vector{ComplexF64}, 
        Tuple{SphericalHarmonics.ML{SphericalHarmonics.ZeroTo{false}, 
            SphericalHarmonics.FullRange{true}}}
    }}(undef, graph_dict["coordination_nr"])

    # neighbor counter
    neighbor_count = 1

    # loop through bonds to neighbors
    for neighbor in MetaGraphsNext.neighbor_labels(
                        graph_dict["spatial_network"], central_vertex)

        # get vector from vertex to neighbor
        vertex_to_neighbor_vec = (sign(neighbor - central_vertex) 
                * graph_dict["spatial_network"][central_vertex, neighbor]["vector"] )

        # get vector's spherical coordinates
        r_length, theta, phi = convert_cartesian_to_spherical(vertex_to_neighbor_vec)

        # get array of spherical harmonics
        y_spherical_harmonic_arr_vec[neighbor_count] = SphericalHarmonics.computeYlm(
                                                            theta, phi; lmax = l_max)

        neighbor_count += 1
    end

    # average over bonds to neighbors

    # initialize dict of Steinhardt order parameters
    q_lm_averaged_bonds_to_neighbors_dict = Dict{Tuple{Int64, Int64}, Complex{Float64}}()

    # loop through values of l
    for l in 0:l_max

        # loop through values of m
        for m in -l:l

            # initialize current average steinhardt order parameter
            q_lm = 0

            # average over neighbors
            for neighbor_count in 1:graph_dict["coordination_nr"]

                q_lm += (1/graph_dict["coordination_nr"]
                    * 
                    y_spherical_harmonic_arr_vec[neighbor_count][(l,m)]
                )

            end

            # save current average steinhardt order parameter
            q_lm_averaged_bonds_to_neighbors_dict[(l,m)] = q_lm
        end

    end

    return q_lm_averaged_bonds_to_neighbors_dict

end


"""
Get vector of mean values of q_l (rotationally invariant Steinhardt local 
bond order parameters) for a 
single vertex and for all parameters l up to l_max where l is the index of 
the spherical harmonic Y_{lm}. The equations are taken from references
10.1103/PhysRevB.28.784 and 10.1063/1.2977970
"""
function get_q_l_averaged_single_vertex_dict(graph_dict::Dict,
    central_vertex::Int64,
    l_max::Int64)

    # get average steinhardt order parameter dict for current vertex
    # this quantity is nor yet rotationally invariant
    q_lm_averaged_bonds_to_neighbors_dict = (
        get_q_lm_averaged_bonds_to_neighbors_dict_single_vertex(graph_dict,
                                                                central_vertex,
                                                                l_max))

    # initialize dict of Steinhardt order parameters 
    q_l_averaged_single_vertex_dict = Dict{Int64, Float64}()

    # loop through values of l
    for l in 0:l_max

        # initialize sum over m
        q_lm_squared_sum = 0

        # loop through values of m
        for m in -l:l

            # add absolute value squared of q_lm to sum
            q_lm_squared_sum += abs2(q_lm_averaged_bonds_to_neighbors_dict[(l,m)])
        end

        # calculate steinhardt order parameter for current l
        q_l_averaged_single_vertex_dict[l] = sqrt( 4*pi / (2*l + 1) * q_lm_squared_sum )

    end

    return q_l_averaged_single_vertex_dict
end


"""
Get vector of mean values of q_l (rotationally invariant Steinhardt local 
bond order parameters) for the entire network and for all parameters l up 
to l_max where l is the index of the spherical harmonic Y_{lm}.
"""
function get_q_l_total_network_mean_dict(graph_dict::Dict,
    l_max::Int64)

    # initialize dictionary of q_l averaged over entire network with all values
    # set to 0
    q_l_total_network_mean_dict = Dict{Int64, Float64}()

    for l in 0:l_max
        q_l_total_network_mean_dict[l] = 0.0

    end

    # loop through vertices
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        # get vector of steinhardt order parameters for current vertex
        q_l_averaged_single_vertex_dict = (
            get_q_l_averaged_single_vertex_dict(
                graph_dict,
                vertex,
                l_max))

        # for each l, add current vertex' contribution to sum of all vertices
        for l in 0:l_max
            q_l_total_network_mean_dict[l] += (1/graph_dict["nr_vertices"] 
                                        * q_l_averaged_single_vertex_dict[l])
    
        end

    end

    return q_l_total_network_mean_dict
end

