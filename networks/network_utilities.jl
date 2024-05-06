"""
These functions are used in different processes when working with
the network graphs
"""


"""
Calculate vector pointing from position a to position b 
under periodic boundary conditions
"""
function get_distance_vector_pbc(position_a::Vector,
                                position_b::Vector,
                                supercell_edge_length::Real )

    # get vector pointing from position a to b without considering boundary conditions
    distance_vector_without_pbc = (position_b .- position_a)

    # modify the vector according to boundary conditions
    distance_vector = ( ( abs.(distance_vector_without_pbc) 
                    .< (supercell_edge_length/2) ) 
                .* distance_vector_without_pbc
            .+ ( abs.(distance_vector_without_pbc .+ supercell_edge_length) 
                    .< (supercell_edge_length/2) ) 
                .* (distance_vector_without_pbc .+ supercell_edge_length)
            .+ ( abs.(distance_vector_without_pbc .- supercell_edge_length) 
                    .< (supercell_edge_length/2) ) 
                .* (distance_vector_without_pbc .- supercell_edge_length)
            )    

    return distance_vector
    
end


"""
Calculate virtual position of an vertex relative to a central vertex
by placing it outside of the supercell if periodic boundary conditions
have to be taken into account
"""
function get_virtual_position(central_vertex_position::Vector{Float64},
                                other_vertex_position::Vector{Float64},
                                supercell_edge_length::Real )

    # get vector pointing from central vertex to neighbor without considering 
    # boundary conditions
    distance_vector_without_pbc = (other_vertex_position .- central_vertex_position)

    # get the other vertices virtual position according to boundary conditions
    virtual_other_vertex_position = ( ( abs.(distance_vector_without_pbc) 
                    .< (supercell_edge_length/2) ) 
                .* other_vertex_position
            .+ ( abs.(distance_vector_without_pbc .+ supercell_edge_length) 
                    .< (supercell_edge_length/2) ) 
                .* (other_vertex_position .+ supercell_edge_length)
            .+ ( abs.(distance_vector_without_pbc .- supercell_edge_length) 
                    .< (supercell_edge_length/2) ) 
                .* (other_vertex_position .- supercell_edge_length)
            )    

    return virtual_other_vertex_position
    
end


"""
Get matrix with the coordinates of an vertices neighbors
by taking periodic boundary conditions into account
"""
function get_neighbor_positions_mat(graph_dict::Dict, central_vertex::Int64;
                                    exclude_vertices::Vector = [])

    # get central vertex's position
    central_vertex_position = graph_dict["spatial_network"][central_vertex]["position"]

    # create matrix to store neighbors coordinates
    neighbor_positions_mat = Matrix{Float64}(undef, 
                            graph_dict["nr_dimensions"],
                            graph_dict["coordination_nr"]-length(exclude_vertices))
    
    # save coordinates to matrix and array
    current_neighbor = 1

    for neighbor in MetaGraphsNext.neighbor_labels(
        graph_dict["spatial_network"], central_vertex)

        if !in(neighbor, exclude_vertices)

            # get neighbor's virtual coordinates which might be outside of the 
            # supercell if periodic boundary conditions play a role
            neighbor_positions_mat[:,current_neighbor] = get_virtual_position(
                            central_vertex_position,
                            graph_dict["spatial_network"][neighbor]["position"],
                            graph_dict["supercell_edge_length"] )

            current_neighbor += 1
        end

    end

    return neighbor_positions_mat

end


"""
Get array with the coordinates of an vertices next to nearest neighbors
by taking periodic boundary conditions into account
"""
function get_next_neighbor_positions_arr(graph_dict::Dict, central_vertex::Int64)

    # get central vertex's position
    central_vertex_position = graph_dict["spatial_network"][central_vertex]["position"]

    # get central vertices neighbors 
    neighbor_vec = collect(MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], central_vertex))

    # create array to store next to nearest neighbors coordinates
    # The first array index labels the number of the direct neighbor
    next_neighbor_positions_arr = Array{Float64}(undef, 
                                                graph_dict["coordination_nr"],
                                                graph_dict["nr_dimensions"],
                                                graph_dict["coordination_nr"]-1)
    
    # loop through central vertices neighbors
    for i in 1:graph_dict["coordination_nr"]

        current_next_neighbor = 1

        # loop through the current neighbor's neighbors
        for next_neighbor in MetaGraphsNext.neighbor_labels(
                                        graph_dict["spatial_network"], neighbor_vec[i])

            if next_neighbor !== central_vertex

                # get next neighbor's virtual coordinates which might be outside of the 
                # supercell if periodic boundary conditions play a role
                next_neighbor_positions_arr[i,:,current_next_neighbor] = get_virtual_position(
                            central_vertex_position,
                            graph_dict["spatial_network"][next_neighbor]["position"],
                            graph_dict["supercell_edge_length"] )

                current_next_neighbor += 1
            end

        end

    end

    return next_neighbor_positions_arr

end



"""
get all bonds inside and on the edge of cluster
"""
function get_cluster_bonds_vec(graph_dict::Dict,
                                cluster_vertices_to_move_vec::Vector{Int64},
                                cluster_vertices_outer_shell_vec::Vector{Int64})

    # initialize vectors for bonds
    cluster_bonds_inside_vec = []
    cluster_bonds_edge_vec = []

    # get vector of all cluster vertices
    all_cluster_vertices_vec = vcat(cluster_vertices_to_move_vec, 
                                cluster_vertices_outer_shell_vec)

    # loop through all cluster vertices
    for cluster_vertex in all_cluster_vertices_vec

        # loop through all neighbors of current vertex
        for neighbor in MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"],
                                                            cluster_vertex)

            # add current bond to the respective vector if it is not stored yet
            if neighbor in all_cluster_vertices_vec

                if cluster_vertex < neighbor
                    push!(cluster_bonds_inside_vec, (cluster_vertex, neighbor) )
                end

            else
                push!(cluster_bonds_edge_vec, Tuple(sort([cluster_vertex, neighbor])) )
            
            end
    
        end

    end

    return [cluster_bonds_inside_vec, cluster_bonds_edge_vec]
end



"""
Get a dictionary of vectors containing all vertices neighboring
the central vertices up to the given shell
"""
function get_cluster_in_shells_dict(graph_dict::Dict, 
                                    central_vertices::Tuple; 
                                    shell_nr::Int64 = 5,
                                    calculate_cluster_energy::Bool = true)

    # initialize dictionary for all neighbors sorted by shells
    cluster_in_shells_dict = Dict(0 => copy(collect(central_vertices)) )

    # initialize vector for cluster vertices which will be allowed to move
    cluster_vertices_to_move_vec = collect(central_vertices)

    # initialize vector for vertices in the outer shell which will not be
    # allowed to move
    cluster_vertices_outer_shell_vec = Vector{Int64}()

    # loop through neighbor shells
    for current_shell in 1:shell_nr

        # initialize vector of vertices in current shell
        current_shell_vertices_vec = Vector{Int64}()

        # loop through vertices of lower shell
        for lower_shell_vertex in cluster_in_shells_dict[current_shell-1]

            # loop through neighbors of current vertex
            for neighbor in MetaGraphsNext.neighbor_labels(
                            graph_dict["spatial_network"], lower_shell_vertex)

                # if not in outer shell, save as vertex to move
                if current_shell < shell_nr

                    # save current vertex if it was not considered before
                    if !(neighbor in cluster_vertices_to_move_vec)

                        push!(current_shell_vertices_vec, neighbor)
                        push!(cluster_vertices_to_move_vec, neighbor)
                    end

                # if in outer shell, the vertex is not allowed to move
                else
                    if (!(neighbor in cluster_vertices_to_move_vec) 
                            && !(neighbor in cluster_vertices_outer_shell_vec))

                        push!(current_shell_vertices_vec, neighbor)
                        push!(cluster_vertices_outer_shell_vec, neighbor)
                    end
                end

            end

        end

        # save to cluster in shells dict
        cluster_in_shells_dict[current_shell] = current_shell_vertices_vec

    end

    # get all bonds inside and on the edge of cluster
    cluster_bonds_inside_vec, cluster_bonds_edge_vec = get_cluster_bonds_vec(
                                                graph_dict,
                                                cluster_vertices_to_move_vec,
                                                cluster_vertices_outer_shell_vec
                                                )

    # collect cluster information into dictionary
    cluster_dict = Dict("cluster_in_shells_dict" => cluster_in_shells_dict, 
            "cluster_vertices_to_move_vec" => cluster_vertices_to_move_vec, 
            "cluster_vertices_outer_shell_vec" => cluster_vertices_outer_shell_vec,
            "cluster_bonds_inside_vec" => cluster_bonds_inside_vec, 
            "cluster_bonds_edge_vec" => cluster_bonds_edge_vec
            )

    # add cluster energy to dictionary if desired
    if calculate_cluster_energy
        cluster_dict["cluster_energy"] = get_cluster_energy(graph_dict, cluster_dict)
        cluster_dict["cluster_energy_up_to_date"] = true
    else
        cluster_dict["cluster_energy_up_to_date"] = false
    end

    return cluster_dict

end


"""
Check if there are ring up the given member number containing the basis vertex
"""
function has_ring_up_to_member_nr(graph_dict::Dict, basis_vertex::Int64; 
                                max_member_nr::Int64 = 4)

    # get the maximal shell nr where a ring of given member nr could be closed
    maximal_shell_nr = Int(floor(max_member_nr/2))

    # get cluster up to maximal shell nr
    cluster_dict = get_cluster_in_shells_dict(graph_dict, 
                Tuple(basis_vertex); 
                shell_nr = maximal_shell_nr,
                calculate_cluster_energy = false)

    # loop through neighbor shells
    for current_shell in 1:maximal_shell_nr

        # loop through vertices of current shell
        for current_vertex in cluster_dict["cluster_in_shells_dict"][current_shell]

            # get vector of current vertex' neighbors
            current_neighbor_vec = collect(MetaGraphsNext.neighbor_labels(
                graph_dict["spatial_network"],  current_vertex))

            # check if there are two connections to lower shell vertices if current shell 
            # is not the first one
            if current_shell > 1
                if sum(current_neighbor in cluster_dict["cluster_in_shells_dict"][current_shell-1]
                        for current_neighbor in current_neighbor_vec) > 1
                    return true
                end
            end

            # check connection to current shell vertices if max member nr is odd
            if isodd(max_member_nr)
                if sum(current_neighbor in cluster_dict["cluster_in_shells_dict"][current_shell]
                            for current_neighbor in current_neighbor_vec) > 0
                    return true
                end
            end
        end
    end

    # if no ring was found, return false
    return false
end


"""
Check if a proposed bond switch introduces a ring of given member nr
"""
function introduces_ring_up_to_member(graph_dict::Dict, 
    switched_chain::Tuple{Int64, Int64, Int64, Int64}; 
    max_member_nr::Int64 = 4)

    # initialize vector to save current bonds' information 
    current_bonds_info_vec = Vector{Dict{String, Any}}(undef, 2)

    # perform trial bond switch
    for i in 1:2
        # save current bond info
        current_bonds_info_vec[i] = graph_dict["spatial_network"][switched_chain[2*i-1], 
                                                            switched_chain[2*i]]

        # remove current bond
        MetaGraphsNext.rem_edge!(graph_dict["spatial_network"],
            switched_chain[2*i-1], switched_chain[2*i])

        # create new trial bond with meaningless data
        graph_dict["spatial_network"][switched_chain[i], switched_chain[2+i]] = Dict(
            "a" => 1)
    end

    # check if new bonds are part of rings 
    has_ring = (has_ring_up_to_member_nr(graph_dict, switched_chain[1]; 
                            max_member_nr = max_member_nr)
            || has_ring_up_to_member_nr(graph_dict, switched_chain[end]; 
                            max_member_nr = max_member_nr))

    # reverse trial bond switch 
    for i in 1:2

        # remove trial bond
        MetaGraphsNext.rem_edge!(graph_dict["spatial_network"],
        switched_chain[i], switched_chain[2+i])

        # recreate previous bond with its data
        graph_dict["spatial_network"][switched_chain[2*i-1], 
        switched_chain[2*i]] = current_bonds_info_vec[i]
    end

    return has_ring
end


"""
Get all remaining 4-vertex-chains that have not been declined yet
and where the first index label is lower than
the last index label in order to not count chains twice
"""
function get_remaining_chains(graph_dict::Dict,
        declined_chains::Vector;
        min_ring_size::Int64 = 5)

    # initialize vector of remaining chains
    remaining_chains = []

    # loop through all possible chains
    for first_vertex in 1:graph_dict["nr_vertices"]

        for second_vertex in setdiff(
            MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], first_vertex), 
            first_vertex)

            for third_vertex in setdiff(
                MetaGraphsNext.neighbor_labels(
                graph_dict["spatial_network"], second_vertex), 
                first_vertex, second_vertex)

                for fourth_vertex in setdiff(
                    MetaGraphsNext.neighbor_labels(
                    graph_dict["spatial_network"], third_vertex), 
                    first_vertex, second_vertex, third_vertex)

                    # creat current chain
                    current_chain = (
                        first_vertex, second_vertex, third_vertex, fourth_vertex)

                    # save current chain to remaining chains if it has not been
                    # declined yet and if neither 1-3 nor 2-4 are already connected
                    if (fourth_vertex > first_vertex 
                        && !(current_chain in declined_chains)
                        && !(current_chain[1] in collect(MetaGraphsNext.neighbor_labels(
                        graph_dict["spatial_network"], current_chain[3])))
                        && !(current_chain[2] in collect(MetaGraphsNext.neighbor_labels(
                            graph_dict["spatial_network"], current_chain[4])))
                        && !(introduces_ring_up_to_member(graph_dict, 
                        current_chain; 
                        max_member_nr = min_ring_size-1) )
                            )

                        push!(remaining_chains, current_chain)
                    end
                end
            end
        end
    end

    return remaining_chains
    
end



"""
Pick a random chain of four vertices that has not been declined since
the last accepted move and where the first index label is lower than
the last index label in order to not count chains twice
"""
function get_random_chain(graph_dict::Dict; 
        declined_chains::Vector = [], 
        remaining_chains::Vector = [], seed = nothing,
        min_ring_size::Int64 = 5)

    # set seed if desired
    if seed !== nothing
        Random.seed!(seed)
    end

    # if a vector of remaining chains is passed, pick one of those
    if remaining_chains != []
        random_chain = rand(remaining_chains)

    # otherwise get random chain without listing all bonds
    else

        # pick a random vertex 1
        random_chain = [rand(1:graph_dict["nr_vertices"])]

        # pick 3 further vertices that sit along a chain
        for i in 1:3
            push!(random_chain, 
                rand(setdiff(
                MetaGraphsNext.neighbor_labels(
                graph_dict["spatial_network"], random_chain[end]), 
                random_chain)) )

        end

        # sort chain such that first index label is lower than last one
        if random_chain[1] > random_chain[end]
            random_chain = reverse(random_chain)
        end

        # create tuple
        random_chain = Tuple(random_chain)

        # check that new chain has not already been declined and
        # that neither 1 and 3 nor 2 and 4 are already connected 
        # and that no ring of given member nr is introduced
        if (random_chain in declined_chains
            || random_chain[1] in collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], random_chain[3]))
            || random_chain[2] in collect(MetaGraphsNext.neighbor_labels(
                graph_dict["spatial_network"], random_chain[4]))
            || introduces_ring_up_to_member(graph_dict, 
            random_chain; 
            max_member_nr = min_ring_size-1)
            )
            random_chain = get_random_chain(graph_dict; declined_chains = declined_chains)
        else
            random_chain =  random_chain
        end

    end

    return random_chain
end


"""
Check whether all vertices in network have the correct coordination number
"""
function get_incorrectly_coordinated_vertices(graph_dict::Dict)

    incorrectly_coordinated_vertices = []

    # for each vertex, check whether it has the correct coordination nr
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])
        if (length( collect( MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], vertex) ) ) 
            !== graph_dict["coordination_nr"])

            push!(incorrectly_coordinated_vertices, vertex)
        end
    end

    return incorrectly_coordinated_vertices
    
end


"""
This function relaxes a given cluster in two ways, first using the
Newton method which is slower but more accurate, and then more efficiently
but less accurate using the method from 10.1142/S0217984987000065. The two
relaxation methods are compared by plotting the evolution of vertex positions
and cluster energies
"""
function compare_relaxation_methods(original_graph_dict,
    central_cluster_vertices,
    evolution_dict,
    filename;
    nr_max_relaxation_cycles = 25,
    shell_nr::Int64 = 4,
    save_path = raw"..\plots\random_networks\\" )

    # initialize arrays for vertex positions and cluster energy as a
    # function of relaxation cycle
    vertex_position_arr = Array{Float64}(undef, nr_max_relaxation_cycles+1, 3, 2)
    cluster_energy_arr = Array{Float64}(undef, nr_max_relaxation_cycles+1, 2)

    # loop through relaxation methods
    relax_efficiently_vec = [false, true]

    for i in 1:2

        # reset graph dict to original one
        graph_dict = deepcopy(original_graph_dict)

        # get the cluster dict
        cluster_dict = get_cluster_in_shells_dict(
                                            graph_dict, 
                                            central_cluster_vertices; 
                                            shell_nr = shell_nr)

        # store initial vertex position and cluster energy
        vertex_position_arr[1,:,i] = graph_dict["spatial_network"][central_cluster_vertices[1]]["position"]
        cluster_energy_arr[1,i] = cluster_dict["cluster_energy"]

        # perform relaxation cycles
        for relaxation_cycle in 2:nr_max_relaxation_cycles+1

            evolution_dict["relax_efficiently"] = relax_efficiently_vec[i]

            graph_dict, cluster_dict = relax_cluster_one_cycle_keating!(graph_dict, 
            cluster_dict,
            evolution_dict;
            update_cluster_energy = true )

            # keep track of vertex position and cluster energy
            vertex_position_arr[relaxation_cycle,:,i] = graph_dict["spatial_network"][central_cluster_vertices[1]]["position"]
            cluster_energy_arr[relaxation_cycle,i] = cluster_dict["cluster_energy"]

        end
    end

    # plot evolution of vertex position and cluster energy
    Plots.plot(xlabel= Latex.L"x",
    ylabel=Latex.L"y" ,
    legend = true, dpi=250)
    Plots.plot!(vertex_position_arr[:,1,1], vertex_position_arr[:,2,1], label = "Newton", markershape = :auto)
    Plots.plot!(vertex_position_arr[:,1,2], vertex_position_arr[:,2,2], label = "efficient", markershape = :auto)

    Plots.savefig(save_path*filename*"_x_y_pos.png")

    Plots.plot(xlabel= Latex.L"x",
    ylabel=Latex.L"z" ,
    legend = true, dpi=250)
    Plots.plot!(vertex_position_arr[:,1,1], vertex_position_arr[:,3,1], label = "Newton", markershape = :auto)
    Plots.plot!(vertex_position_arr[:,1,2], vertex_position_arr[:,3,2], label = "efficient", markershape = :auto)

    Plots.savefig(save_path*filename*"_x_z_pos.png")

    Plots.plot(xlabel= "relaxation cycle",
    ylabel="cluster energy" ,
    legend = true, dpi=250)
    Plots.plot!(collect(0:nr_max_relaxation_cycles), cluster_energy_arr[:,1], label = "Newton")
    Plots.plot!(collect(0:nr_max_relaxation_cycles), cluster_energy_arr[:,2], label = "efficient")

    Plots.savefig(save_path*filename*"_cluster_energy.png")

    return [vertex_position_arr, cluster_energy_arr]
end



"""
Get extremum of a quadratic function y= a*x^2 + c given by two points
(x_vec[1], y_vec[1]) and (x_vec[2], y_vec[2])
"""
function get_energy_relaxation_coefficients(x_vec::Vector, y_vec::Vector)

    # get coefficients of quadratic function
    a = (y_vec[1] - y_vec[2]) / (x_vec[1]^2 - x_vec[2]^2)
    c = (-x_vec[2]^2 * y_vec[1] + x_vec[1]^2 * y_vec[2])/(x_vec[1]^2 - x_vec[2]^2)

    return [a, c]
end