"""
These functions modify the network spatial networks, e. g. by a bond switch
"""

"""
This function performs a bond switch on a spatial network. The argument
switched_chain is a tuple of four integers which is the edge type of the
MetaGraphsNext package
"""
function switch_chain!(
    spatial_network::MetaGraphsNext.MetaGraph,
    switched_chain::Tuple{Int64, Int64, Int64, Int64} )

    # remove two old bonds and create new ones
    for i in 1:2

        # remove old bond
        MetaGraphsNext.rem_edge!(spatial_network,
            switched_chain[2*i-1], switched_chain[2*i])

        # initialize vector along new bond
        vector_along_new_bond = Vector{Float64}(undef,
            spatial_network[]["nr_dimensions"])

        # determine vector along new bond
        if switched_chain[i] < switched_chain[2+i]
            vector_along_new_bond = get_distance_vector_pbc(
                spatial_network[switched_chain[i]]["position"],
                spatial_network[switched_chain[2+i]]["position"],
                spatial_network[]["supercell_edge_length"] )
        else
            vector_along_new_bond = get_distance_vector_pbc(
                spatial_network[switched_chain[2+i]]["position"],
                spatial_network[switched_chain[i]]["position"],
                spatial_network[]["supercell_edge_length"] )
        end

        # determine length of vector along new bond
        distance_along_new_bond = LinearAlgebra.norm(vector_along_new_bond)

        # create new bond
        spatial_network[switched_chain[i], switched_chain[2+i]] = Dict(
            "vector" => vector_along_new_bond, 
            "distance_squared" => distance_along_new_bond^2 )
    end

    # note, that total energy is not up to date any more
    spatial_network[]["total_energy_up_to_date"] = false

    return spatial_network
end


"""
move an vertex and update its energy and the edges to its neighbors
"""
function move_vertex!(
    spatial_network::MetaGraphsNext.MetaGraph, 
    vertex_to_move::Int64, 
    translation_vector::Vector{Float64};
    update_total_energy::Bool = false,
    large_translation::Bool = false)

    # update vertex position by taking periodic boundary conditions into
    # account
    initial_position = spatial_network[vertex_to_move]["position"]
    spatial_network[vertex_to_move]["position"] = (initial_position
            .+ translation_vector .+ spatial_network[]["supercell_edge_length"]
        ).%spatial_network[]["supercell_edge_length"]

    # update outgoing edges                                                    
    for neighbor in MetaGraphsNext.neighbor_labels(
        spatial_network, vertex_to_move)

        original_distance_vector = spatial_network[vertex_to_move, neighbor][
            "vector"]

        # determine new vector, where the direction of the vector (from vertex
        # with lower label to vertex with higher label) has to be taken into
        # account. If the translation is large, periodic boundary conditions
        # have to be taken into account which is less efficient
        if large_translation
            bond_vector = get_distance_vector_pbc(
                spatial_network[vertex_to_move]["position"],
                spatial_network[neighbor]["position"],
                spatial_network[]["supercell_edge_length"] )
        else
            spatial_network[vertex_to_move, neighbor]["vector"] = (
            original_distance_vector 
            .- float(sign(neighbor - vertex_to_move)) .* translation_vector )
        end
        

        spatial_network[vertex_to_move, neighbor]["distance_squared"] = (
            LinearAlgebra.norm(spatial_network[vertex_to_move, neighbor][
                "vector"])^2 )
    end
    
    # update total energy if desired
    if update_total_energy
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
        spatial_network[]["total_energy_up_to_date"] = true
    else
        spatial_network[]["total_energy_up_to_date"] = false
    end

    return spatial_network
end


"""
calculate translation vector to approximate energy minimum of Keating strain
energy as explained in 10.1142/S0217984987000065
"""
function get_approximate_translation_vector_keating(
    gradient::Vector{Float64}, 
    bond_bending_const::Float64;
    relaxation_overshoot_factor_r = 1.5,
    relaxation_optimization_parameter_l = 1)

    # determine translation vector
    translation_vector = ((- relaxation_overshoot_factor_r/(4 
            + 5*bond_bending_const 
            + relaxation_optimization_parameter_l
                *LinearAlgebra.norm(gradient)))
        .* gradient)

    return translation_vector
end


"""
Approximately relax a single vertex by efficiently determining the approximated
coordinate shift. The corresponding methdo is explained in 
10.1142/S0217984987000065
"""
function relax_single_vertex_keating_efficiently!(
    spatial_network::MetaGraphsNext.MetaGraph,
    vertex_to_relax::Int64;
    relaxation_overshoot_factor_r = 1.5,
    relaxation_optimization_parameter_l = 1,
    update_total_energy::Bool = false)

    # get energy gradient at current vertex position
    gradient = gradient_keating_efficient(spatial_network, vertex_to_relax)

    # calculate translation vector to approximate energy minimum
    translation_vector =  get_approximate_translation_vector_keating(gradient, 
        spatial_network[]["bond_bending_const"];
        relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
        relaxation_optimization_parameter_l 
            = relaxation_optimization_parameter_l)

    # move vertex 
    spatial_network = move_vertex!(spatial_network, vertex_to_relax, 
        translation_vector; update_total_energy = update_total_energy)

    return [spatial_network, gradient]
end


"""
Perform one cycle of relaxing a cluster of vertices 
"""
function relax_cluster_one_cycle_keating!(
    spatial_network::MetaGraphsNext.MetaGraph, 
    cluster_dict::Dict,
    evolution_dict::Dict;
    update_total_energy::Bool = false,
    update_cluster_energy::Bool = true )

    # get initial cluster energy
    if cluster_dict["cluster_energy_up_to_date"]
        initial_cluster_energy = cluster_dict["cluster_energy"]
    else
        initial_cluster_energy = get_cluster_energy(spatial_network, 
            cluster_dict)
    end

    # get total energy if it's not up to date but supposed to be updated later
    if update_total_energy && !spatial_network[]["total_energy_up_to_date"]
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
    end

    # reset absolute value of total cluster force
    cluster_dict["cluster_force"] = 0.0

    # relax each vertex in the given cluster
    for vertex in cluster_dict["cluster_vertices_to_move_vec"]

        # relax efficiently but approximately or exactly but slowly
        spatial_network, gradient = (
            relax_single_vertex_keating_efficiently!(spatial_network,
                vertex;
                relaxation_overshoot_factor_r 
                    =evolution_dict["relaxation_overshoot_factor_r"],
                relaxation_optimization_parameter_l
                    =evolution_dict["relaxation_optimization_parameter_l"],
                update_total_energy = false))

        # add absolute value of force on current vertex to total cluster force
        cluster_dict["cluster_force"] += LinearAlgebra.norm(gradient)
    end

    # update cluster energy if desired
    if update_cluster_energy
        cluster_dict["cluster_energy"] = get_cluster_energy(
            spatial_network, cluster_dict)
        cluster_dict["cluster_energy_up_to_date"] = true
    else
        cluster_dict["cluster_energy_up_to_date"] = false
    end

    # update total energy if desired
    if update_total_energy

        if cluster_dict["cluster_energy_up_to_date"]
            cluster_energy = cluster_dict["cluster_energy"] 
        else
            cluster_energy = get_cluster_energy(spatial_network, cluster_dict)
        end

        spatial_network[]["total_energy"] = (spatial_network[]["total_energy"] 
            + cluster_energy - initial_cluster_energy)

        spatial_network[]["total_energy_up_to_date"] = true
    else
        spatial_network[]["total_energy_up_to_date"] = false
    end

    return [spatial_network, cluster_dict] 
end


"""
Relax the given cluster using an optimization package as opposed to a manual
vertex-by-vertex implementation
"""
function relax_network_keating_optim!(
    spatial_network::MetaGraphsNext.MetaGraph,
    cluster_dict::Dict;
    optimization_method::String = "lbfgs",
    update_total_energy::Bool = false)

    # set the optimization method for the optimization package
    if optimization_method == "lbfgs"
        optimization_fct = Optim.LBFGS()
    elseif optimization_method == "newton"
        optimization_fct = Optim.Newton()
    else
        @error "Optimization method, specified in evolution dict,
        is not known."
    end

    # get a 1d vector of all vertex coordinates in the cluster in the format
    # [x1, y1, z1, x2, ...]
    vertices_to_move_coord_vec = get_vertex_coord_vec(spatial_network, 
        cluster_dict["cluster_vertices_to_move_vec"])
    vertices_outer_shell_coord_vec = get_vertex_coord_vec(spatial_network, 
        cluster_dict["cluster_vertices_outer_shell_vec"])

    # get all vertex labels that are contained in the cluster bonds on the
    # outer edge but are not in the outer shell cluster vertices
    vertices_bonds_edge_vec::Vector{Int64} = []
    for bond in cluster_dict["cluster_bonds_edge_vec"]
        for vertex in bond
            if !(vertex in cluster_dict["cluster_vertices_outer_shell_vec"])
                push!(vertices_bonds_edge_vec, vertex)
            end
        end
    end
    vertices_bonds_edge_coord_vec = get_vertex_coord_vec(spatial_network, 
        vertices_bonds_edge_vec)

    # create a dictionary for the bonds connected to each vertex
    all_cluster_vertices_vec = vcat(
        cluster_dict["cluster_vertices_to_move_vec"], 
        cluster_dict["cluster_vertices_outer_shell_vec"]) 
    all_cluster_bonds_vec = vcat(
        cluster_dict["cluster_bonds_inside_vec"], 
        cluster_dict["cluster_bonds_edge_vec"])
    bonds_connected_to_vertex_dict = Dict{Int64, Vector{Tuple{Int64, Int64}}}()

    for vertex in all_cluster_vertices_vec
        connected_bonds = filter(x -> (x[1] == vertex) || (x[2] == vertex), 
            all_cluster_bonds_vec)

        bonds_connected_to_vertex_dict[vertex] = connected_bonds
    end

    # create another dictionary that, for each bond, stores the vertex labels 
    # of the neighbors for each of the two vertices in the bond. The values
    # of the dictionary are tuples of two lists of integers
    bond_vertex_neighbor_dict = Dict{
        Tuple{Int64, Int64}, Tuple{Vector{Int64}, Vector{Int64}}}()
    for bond in cluster_dict["cluster_bonds_inside_vec"]
        neighbor_label_vec_1::Vector{Int64} = [x for x in 
            MetaGraphsNext.neighbor_labels(spatial_network, bond[1]) 
            if x != bond[2]]
        neighbor_label_vec_2::Vector{Int64} = [x for x in 
            MetaGraphsNext.neighbor_labels(spatial_network, bond[2]) 
            if x != bond[1]]
        bond_vertex_neighbor_dict[bond] = (neighbor_label_vec_1,
            neighbor_label_vec_2)
    end

    # create a dictionary with the l values for all bonds, that is the least 
    # common multiple of the coordination numbers of the two vertices minus one 
    # since the bond itself is not included
    bond_l_dict = Dict{Tuple{Int64, Int64}, Int64}()
    for bond in cluster_dict["cluster_bonds_inside_vec"]
        vertex_1_coordination_nr = spatial_network[bond[1]]["coordination_nr"]
        vertex_2_coordination_nr = spatial_network[bond[2]]["coordination_nr"]
        bond_l_dict[bond] = lcm(
            vertex_1_coordination_nr-1, vertex_2_coordination_nr-1)
    end

    # optimize the Keating energy using a closure
    result = Optim.optimize(
        coord_vec -> cluster_energy_optim(
            coord_vec, 
            vertices_outer_shell_coord_vec, 
            vertices_bonds_edge_coord_vec,
            vertices_bonds_edge_vec,
            bonds_connected_to_vertex_dict,
            bond_vertex_neighbor_dict,
            bond_l_dict,
            cluster_dict,
            spatial_network),
        vertices_to_move_coord_vec, 
        optimization_fct,
        Optim.Options(g_tol=1e-6);
        autodiff = :forward
    )

    relaxed_coord_vec = Optim.minimizer(result)
    
    spatial_network = update_vertex_coords!(
        spatial_network, cluster_dict["cluster_vertices_to_move_vec"], 
        relaxed_coord_vec)

    # update total energy if desired
    if update_total_energy
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
        spatial_network[]["total_energy_up_to_date"] = true
    else
        spatial_network[]["total_energy_up_to_date"] = false
    end

    return spatial_network
end


"""
Relax the given cluster using an optimization package as opposed to a manual
vertex-by-vertex implementation
"""
function relax_network_keating_optim_inefficient!(
    spatial_network::MetaGraphsNext.MetaGraph,
    cluster_dict::Dict;
    optimization_method::String = "lbfgs",
    update_total_energy::Bool = false)

    # set the optimization method for the optimization package
    if optimization_method == "lbfgs"
        optimization_fct = Optim.LBFGS()
    elseif optimization_method == "newton"
        optimization_fct = Optim.Newton()
    else
        @error "Optimization method, specified in evolution dict,
        is not known."
    end

    # get a 1d vector of all vertex coordinates in the cluster in the format
    # [x1, y1, z1, x2, ...]
    vertices_to_move_coord_vec = get_vertex_coord_vec(spatial_network, 
        cluster_dict["cluster_vertices_to_move_vec"])
    vertices_outer_shell_coord_vec = get_vertex_coord_vec(spatial_network, 
        cluster_dict["cluster_vertices_outer_shell_vec"])

    # get all vertex labels that are contained in the cluster bonds on the
    # outer edge but are not in the outer shell cluster vertices
    vertices_bonds_edge_vec::Vector{Int64} = []
    for bond in cluster_dict["cluster_bonds_edge_vec"]
        for vertex in bond
            if !(vertex in cluster_dict["cluster_vertices_outer_shell_vec"])
                push!(vertices_bonds_edge_vec, vertex)
            end
        end
    end
    vertices_bonds_edge_coord_vec = get_vertex_coord_vec(spatial_network, 
        vertices_bonds_edge_vec)

    # optimize the Keating energy using a closure
    result = Optim.optimize(
        coord_vec -> cluster_energy_optim_inefficient(
            coord_vec, 
            vertices_outer_shell_coord_vec, 
            vertices_bonds_edge_coord_vec,
            vertices_bonds_edge_vec,
            cluster_dict,
            spatial_network),
        vertices_to_move_coord_vec, 
        optimization_fct,
        Optim.Options(g_tol=1e-6);
        autodiff = :forward
    )

    relaxed_coord_vec = Optim.minimizer(result)
    
    spatial_network = update_vertex_coords!(
        spatial_network, cluster_dict["cluster_vertices_to_move_vec"], 
        relaxed_coord_vec)

    # update total energy if desired
    if update_total_energy
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
        spatial_network[]["total_energy_up_to_date"] = true
    else
        spatial_network[]["total_energy_up_to_date"] = false
    end

    return spatial_network
end


"""
Fully relax a cluster of vertices. The cluster energy will always be updated
"""
function relax_network_keating!(
    spatial_network::MetaGraphsNext.MetaGraph,
    switched_chain::Tuple{Int64, Int64, Int64, Int64},
    evolution_dict::Dict;
    threshold_total_energy = Inf,
    update_total_energy::Bool = false,
    print_progress::Bool = false)

    # get small cluster of vertices around the given chain
    cluster_dict = get_cluster_in_shells_dict(spatial_network, switched_chain,
        shell_nr = evolution_dict["shell_nr"], calculate_cluster_energy = true)

    # make sure that total energy is up to date
    if !spatial_network[]["total_energy_up_to_date"]
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
    end

    # get threshold cluster energy for initial cluster
    threshold_cluster_energy = (threshold_total_energy
        - spatial_network[]["total_energy"] + cluster_dict["cluster_energy"])
    
    # relax efficiently but approximately
    if evolution_dict["relax_efficiently"]

        # if network is supposed to be relaxed globally, store initial shell nr
        # and set threshold cycle for global relaxation
        if (haskey(evolution_dict, "relax_globally_after_threshold_cycle")
                && evolution_dict["relax_globally_after_threshold_cycle"])
            threshold_cycle_global_relaxation = (
                evolution_dict["reject_during_relaxation_cycle_threshold"]*2+1)
        else
            threshold_cycle_global_relaxation = evolution_dict[
                "nr_max_relaxation_cycles"] + 1
        end

        # create bool that determines whether cluster will be relaxed globally
        # after local relaxation
        continue_with_global_relaxation = true


        # perform the given number of relaxation cycles
        for cycle_nr in 1:threshold_cycle_global_relaxation-1

            # only update cluster energy, if this is needed to get cluster 
            # energy change
            if cycle_nr <= evolution_dict[
                "reject_during_relaxation_cycle_threshold"]-1

                update_cluster_energy = false
            else
                update_cluster_energy = true

                # store previous cluster force and energy
                previous_cluster_force = cluster_dict["cluster_force"]
                previous_cluster_energy = cluster_dict["cluster_energy"]
            end

            # relax cluster for one cycle
            spatial_network, cluster_dict = relax_cluster_one_cycle_keating!(
                spatial_network, 
                cluster_dict,
                evolution_dict;
                update_total_energy = false,
                update_cluster_energy = update_cluster_energy )

            # if cycle nr is above the given threshold, check if the relaxation 
            # can be breaked when it becomes clear that the total energy will 
            # exceed the threshold or because of small relative energy change
            if cycle_nr > evolution_dict[
                "reject_during_relaxation_cycle_threshold"]

                # get vector of last two cluster forces
                cluster_force_vec = [previous_cluster_force,
                    cluster_dict["cluster_force"]]

                # get vector of last two cluster energies
                cluster_energy_vec =[previous_cluster_energy,
                    cluster_dict["cluster_energy"]]

                # estimate relaxed cluster energy
                prefactor_force_squared, relaxed_cluster_energy = (
                    get_energy_relaxation_coefficients(
                        cluster_force_vec, cluster_energy_vec))

                if print_progress
                    println("Prefactor force squared: "*string(
                        prefactor_force_squared))
                    println("Relaxed cluster energy: "*string(
                        relaxed_cluster_energy))
                end

                # break if estimated energy change exceeds the given threshold
                if (prefactor_force_squared < 1
                    && relaxed_cluster_energy > 1.05*threshold_cluster_energy) 

                    if print_progress
                        println("Relaxed energy exceeds threshold: breaking at
                            cycle nr "*string(cycle_nr))
                    end

                    # don't continue with global relaxation
                    continue_with_global_relaxation = false

                    break
                end

                # break if cluster energy changes less than the given threshold
                relative_cluster_energy_change = (
                    abs((previous_cluster_energy 
                        - cluster_dict["cluster_energy"])
                    /cluster_dict["cluster_energy"]))

                if relative_cluster_energy_change < evolution_dict[
                    "break_at_relative_cluster_energy_change"] 
                
                    if print_progress
                        println("Negligeable energy change: breaking at cycle nr "
                            *string(cycle_nr))
                    end
                    break
                end
            end
        end
    
        # relax network globally, if desired
        if (haskey(evolution_dict, "relax_globally_after_threshold_cycle") 
            && evolution_dict["relax_globally_after_threshold_cycle"] 
            && continue_with_global_relaxation)

            if print_progress
                println("Relaxing network globally at cycle nr "
                    *string(threshold_cycle_global_relaxation))
            end

            # get cluster of entire network by setting shell nr to a high value
            cluster_dict = get_cluster_in_shells_dict(spatial_network, 
                switched_chain,
                shell_nr = Int(ceil(log(spatial_network[]["nr_vertices"])))+4,
                calculate_cluster_energy = false)

            for cycle_nr in threshold_cycle_global_relaxation:evolution_dict[
                "nr_max_relaxation_cycles"]

                # only update cluster energy, if this is needed to get cluster
                # energy change
                if cycle_nr > threshold_cycle_global_relaxation

                    # store previous cluster force and energy
                    previous_cluster_force = cluster_dict["cluster_force"]
                    previous_cluster_energy = cluster_dict["cluster_energy"]
                end

                # relax cluster for one cycle
                spatial_network, cluster_dict = (
                    relax_cluster_one_cycle_keating!(
                        spatial_network, 
                        cluster_dict,
                        evolution_dict;
                        update_total_energy = false,
                        update_cluster_energy = true ))

                if cycle_nr > threshold_cycle_global_relaxation

                    # get vector of last two cluster forces
                    cluster_force_vec = [previous_cluster_force,
                        cluster_dict["cluster_force"]]
                
                    # get vector of last two cluster energies
                    cluster_energy_vec =[previous_cluster_energy,
                        cluster_dict["cluster_energy"]]
                
                    # estimate relaxed cluster energy
                    prefactor_force_squared, relaxed_cluster_energy = (
                        get_energy_relaxation_coefficients(
                            cluster_force_vec, cluster_energy_vec))
                
                    if print_progress
                        println("Prefactor force squared: "
                            *string(prefactor_force_squared))
                        println("Relaxed cluster energy: "
                            *string(relaxed_cluster_energy))
                    end
                
                    # break if estimated energy change exceeds given threshold
                    if (prefactor_force_squared < 1 && 
                        relaxed_cluster_energy > 1.05*threshold_total_energy) 

                        if print_progress
                            println("Relaxed energy exceeds threshold: breaking
                                at cycle nr "*string(cycle_nr))
                        end
                        break
                    end
                
                    # break if cluster energy changes less than given threshold
                    relative_cluster_energy_change = (
                        abs((previous_cluster_energy
                            - cluster_dict["cluster_energy"])
                         /cluster_dict["cluster_energy"]))
                
                    if relative_cluster_energy_change < evolution_dict[
                        "break_at_relative_cluster_energy_change"] 
                    
                        if print_progress
                            println("Negligeable energy change: breaking at
                            cycle nr "*string(cycle_nr))
                        end
                        break
                    end
                end
            end
        end

    # otherwise, relax the entire cluster in one step using a gradient based
    # method
    else
        # if global relaxation is desired, print an error that this is not 
        # implemented
        if (haskey(evolution_dict, "relax_globally_after_threshold_cycle") 
            && evolution_dict["relax_globally_after_threshold_cycle"])
            @error "Global relaxation is not implemented for Optim based
            relaxation. Please set relax_globally_after_threshold_cycle to 
            false."
        else
            spatial_network = relax_network_keating_optim!(spatial_network,
                cluster_dict;
                optimization_method = 
                    evolution_dict["inefficient_optimization_method"],
                update_total_energy = false)
        end
        
    end

    # update total energy if desired
    if update_total_energy
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
        spatial_network[]["total_energy_up_to_date"] = true
    else
        spatial_network[]["total_energy_up_to_date"] = false
    end
    
    return spatial_network
end


"""
Introduce thermal fluctuations to cluster by moving
each vertex in the cluster according to 10.1063/1.4867897
"""
function excite_cluster!(
    spatial_network::MetaGraphsNext.MetaGraph, 
    cluster_dict::Dict,
    temperature;
    update_total_energy = false,
    update_cluster_energy = false)

    # if total energy will be updated, store initial cluster energy and make
    # sure that total energy is up to date
    if update_total_energy
        if !spatial_network[]["total_energy_up_to_date"]
            spatial_network[]["total_energy"] = get_total_energy_keating(
                spatial_network)
            spatial_network[]["total_energy_up_to_date"] = true
        end

        if cluster_dict["cluster_energy_up_to_date"]
            initial_cluster_energy  = cluster_dict["cluster_energy"]
        else
            initial_cluster_energy = get_cluster_energy(spatial_network,
                cluster_dict)
        end
    end

    # initialize excitation weight
    cluster_excitation_weight = 1

    # initialize dict where vertex displacements will be stored
    vertex_excitation_dict = Dict()

    # loop through all vertices that are allowed to move
    for vertex in cluster_dict["cluster_vertices_to_move_vec"]

        # get hessian matrix of current vertex at relaxed position
        hessian = hessian_keating_efficient(spatial_network, vertex)

        # get vector of sigmas
        sigma_vec = sqrt.(temperature ./ LinearAlgebra.diag(hessian) ) 

        # draw vector of displacements from Gaussian distribution with standard
        # deviations in the sigma vector
        displacement_vec = sigma_vec .* randn(
            spatial_network[]["nr_dimensions"])

        # get excitation weight for current vertex
        vertex_excitation_weight = prod(displacement_vec)

        # multiply it to cluster excitation weight
        cluster_excitation_weight *= vertex_excitation_weight

        # store displacement vector
        vertex_excitation_dict[vertex] = displacement_vec
    end

    # move vertices according to the given displacements
    for vertex in cluster_dict["cluster_vertices_to_move_vec"]

        # move vertices according to the previously thermal fluctuations
        spatial_network = move_vertex!(spatial_network, vertex, 
            vertex_excitation_dict[vertex]; update_total_energy = false)
    end

    # update cluster energy if desired
    if update_cluster_energy
        cluster_dict["cluster_energy"] = get_cluster_energy(spatial_network, 
            cluster_dict)
        cluster_dict["cluster_energy_up_to_date"] = true
    else
        cluster_dict["cluster_energy_up_to_date"] = false
    end

    # if desired update total energy
    if update_total_energy

        # get excited cluster energy
        if cluster_dict["cluster_energy_up_to_date"]
            excited_cluster_energy = cluster_dict["cluster_energy"]
        else
            excited_cluster_energy = get_cluster_energy(spatial_network,
                cluster_dict)
        end

        # update total energy
        spatial_network[]["total_energy"] = (spatial_network[]["total_energy"] 
                                    + excited_cluster_energy
                                    - initial_cluster_energy)

        spatial_network[]["total_energy_up_to_date"] = true
    else
        spatial_network[]["total_energy_up_to_date"] = false
    end

    return [spatial_network, cluster_dict, cluster_excitation_weight]
end


"""
Perform a Monte Carlo move without thermal fluctuations by switching a bond,
relaxing the network and then accepting the network with Metropolis acceptance
probability
"""
function monte_carlo_move!(
    spatial_network::MetaGraphsNext.MetaGraph, 
    evolution_dict::Dict,
    temperature; 
    switched_chain::Tuple{Int64, Int64, Int64, Int64} 
        = get_random_chain(spatial_network))

     # make sure that initial total energy is up to date
    if !spatial_network[]["total_energy_up_to_date"]
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
        spatial_network[]["total_energy_up_to_date"] = true
    end

    # save original spatial network
    initial_spatial_network = deepcopy(spatial_network)

    # set threshold total energy for Metropolis acceptance probability
    # if there are no thermal fluctuations considered
    if !evolution_dict["thermal_fluctuations"]
        threshold_total_energy = (initial_spatial_network[]["total_energy"]
            - temperature * log(rand()))
    else
        threshold_total_energy = Inf
    end

    # if there are thermal fluctuations, relax the network first and calculate
    # weights of the corresponding shifts
    if evolution_dict["thermal_fluctuations"]

        # get cluster of entire network by setting shell nr to the given value
        cluster_dict = get_cluster_in_shells_dict(spatial_network, 
            switched_chain,
            shell_nr =  Int(ceil(log(spatial_network[]["nr_vertices"])))+4,
            calculate_cluster_energy = false)

        # relax cluster
        spatial_network = relax_network_keating!(spatial_network,
            switched_chain, evolution_dict; update_total_energy = false)

        # get cluster weight corresponding to the relaxation translations
        total_relaxation_weight = get_cluster_fluctuation_weight(
            initial_spatial_network, spatial_network, cluster_dict,
            temperature)
    end

    # switch bonds
    spatial_network = switch_chain!(spatial_network, switched_chain)

    # relax total network and only update total energy if there won't be
    # thermal fluctuations included afterward
    spatial_network = relax_network_keating!(
        spatial_network,
        switched_chain,
        evolution_dict;
        threshold_total_energy = threshold_total_energy,
        update_total_energy = !evolution_dict["thermal_fluctuations"])

    # if desired, include thermal fluctuations by randomly shifting all cluster
    # vertices
    if evolution_dict["thermal_fluctuations"]

        # get cluster of entire network by setting shell nr to the given value
        cluster_dict = get_cluster_in_shells_dict(spatial_network,
            switched_chain,
            shell_nr = Int(ceil(log(spatial_network[]["nr_vertices"])))+4,
            calculate_cluster_energy = false)

        # excite cluster with thermal fluctuations and get the corresponding
        # excitation weight
        spatial_network, cluster_dict, total_excitation_weight = (
            excite_cluster!(spatial_network, cluster_dict, temperature;
                update_total_energy = true, update_cluster_energy = false))

        # set random threshold energy increase for Metropolis acceptance
        # probability
        threshold_total_energy = (initial_spatial_network[]["total_energy"]
            - temperature 
            * (total_excitation_weight/total_relaxation_weight) * log(rand()))

    end

    # accept move if energy increase is below threshold
    if spatial_network[]["total_energy"] <= threshold_total_energy
        move_accepted = true
    else
        move_accepted = false
        spatial_network = initial_spatial_network
    end

    return [spatial_network, move_accepted]
end


"""
Evolve the network with a given number of attempted Monte Carlo moves
"""
function evolve_network!(
    spatial_network::MetaGraphsNext.MetaGraph,
    evolution_dict::Dict,
    nr_attempted_bond_switches::Int64, 
    temperature;
    declined_chains::Vector = [],
    remaining_chains::Vector = [],
    total_energy_vec::Vector{Float64} = Vector{Float64}(undef, 0),
    move_accepted_vec::Vector{Bool} = Vector{Bool}(undef, 0),
    print_progress::Bool = false,
    print_every_nr_attempted_bond_switches::Int64 = 1,
    random_evolution_seed::Int64 = -1,
    quench_counter::Int64 = 0,
    print_lock = Threads.ReentrantLock())

    #print("evolve_network")

    # set seed for random evolution if desired
    if random_evolution_seed != -1
        Random.seed!(random_evolution_seed)
    end

    # determine nr of chains of four vertices
    nr_chains=length(get_all_chains(spatial_network))

    # attempt given number of bond switches
    for i in 1:nr_attempted_bond_switches
        
        #print("i",i)

        # get remaining chains if list of declined chains is long and
        # remaining have not been determined yet
        if (length(declined_chains) > 0.4*nr_chains && remaining_chains == [])

            remaining_chains = get_remaining_chains(
                spatial_network,
                declined_chains;
                min_ring_size=evolution_dict["min_ring_size"])

            # break if network is quenched 
            # (all chains have been attempted without success)
            if remaining_chains == []
                println("Network quenched after "*string(i)*"th attempt.")
                break
            end
        end

        # get random chain that hasn't been declined since the last
        # accepted switch
        switched_chain = get_random_chain(
            spatial_network; 
            declined_chains = declined_chains,
            remaining_chains = remaining_chains,
            min_ring_size = evolution_dict["min_ring_size"])

        # print attempted chain if desired
        if print_progress && (print_every_nr_attempted_bond_switches == 1)
            println("Attempt chain "*string(switched_chain))
        end

        # attempt Monte Carlo move
        spatial_network, move_accepted = monte_carlo_move!(
            spatial_network, 
            evolution_dict,
            temperature; 
            switched_chain = switched_chain)

        # update declined and remaining chains vectors
        if move_accepted
            push!(move_accepted_vec, true)
            declined_chains = []
            remaining_chains = []

            # print progress if desired
            if print_progress
                if print_every_nr_attempted_bond_switches == 1
                    println("Chain "*string(switched_chain)
                        *" accepted. Energy: "
                        *string(spatial_network[]["total_energy"]))

                elseif i%print_every_nr_attempted_bond_switches == 0
                    lock(print_lock) do
                        if (haskey(evolution_dict, 
                                "mean_nr_monte_carlo_steps_for_quenching")
                            && temperature == 0 
                            && i == 50*nr_chains)

                            Format.printfmtln(
                                "Thread {1}: attempted bond switch nr {2} at 
                                    T={3} accepted. T={3} is {4:.3f} % done. 
                                    Finished quenches: {5}",
                                Threads.threadid(), i, temperature, 
                                i/(evolution_dict[
                                    "mean_nr_monte_carlo_steps_for_quenching"]
                                *nr_chains)*100,
                                quench_counter )
                        else
                            Format.printfmtln("Thread {1}: attempted bond
                                switch nr {2} at T={3} accepted. T={3} is 
                                {4:.3f} % done. Finished quenches: {5}",
                                Threads.threadid(), i, temperature, 
                                i/nr_attempted_bond_switches*100,
                                quench_counter )
                        end
                    end
                end
            end

        else
            push!(move_accepted_vec, false)
            push!(declined_chains, switched_chain)

            # print progress if desired
            if print_progress
                if print_every_nr_attempted_bond_switches == 1
                    println("Chain "*string(switched_chain)*" declined.")

                elseif i%print_every_nr_attempted_bond_switches == 0
                    lock(print_lock) do
                        if (haskey(evolution_dict,
                            "mean_nr_monte_carlo_steps_for_quenching")
                            && temperature == 0 
                            && i == 50*nr_chains)

                            Format.printfmtln("Thread {1}: attempted bond 
                                switch nr {2} at T={3} declined. T={3} is {4:.3f} 
                                % done. Finished quenches: {5}",
                                Threads.threadid(), i, temperature, 
                                i/(evolution_dict[
                                    "mean_nr_monte_carlo_steps_for_quenching"]
                                    *nr_chains)*100,
                                quench_counter )
                        else
                            Format.printfmtln("Thread {1}: attempted bond 
                                switch nr {2} at T={3} declined. T={3} is 
                                {4:.3f} % done. Finished quenches: {5}",
                                Threads.threadid(), i, temperature, 
                                i/nr_attempted_bond_switches*100,
                                quench_counter )
                        end
                    end
                end
            end

            # break if network is quenched 
            # (all chains have been attempted without success)
            if length(remaining_chains) == 1
                lock(print_lock) do
                    println("Network quenched after "*string(i)*"th attempt.")
                end
                break

            # otherwise update remaining chains
            else
                deleteat!(remaining_chains,
                    findall(x->x==switched_chain,remaining_chains))
            end
            
        end

        # update total energy
        push!(total_energy_vec, spatial_network[]["total_energy"])

    end

    return [spatial_network, total_energy_vec, move_accepted_vec]
end


"""
Evolve network according to a given order of temperatures and nr of Monte Carlo
steps (which I define as nr_bonds attempted Monte Carlo moves) per temperature
"""
function evolve_network_temperature_sequence!(
    spatial_network::MetaGraphsNext.MetaGraph,
    evolution_dict::Dict;
    total_energy_vec::Vector{Float64} = Vector{Float64}(undef, 0),
    move_accepted_vec::Vector{Bool} = Vector{Bool}(undef, 0),
    print_progress::Bool = false,
    print_every_nr_attempted_bond_switches::Int64 = 100,
    random_evolution_seed::Int64 = -1,
    save_network_after_each_temperature::Bool = false,
    filename::String = "sample_name",
    save_path::String = "../../data/structures/",
    print_lock = Threads.ReentrantLock())

    # set seed for random evolution if desired
    if random_evolution_seed != -1
        Random.seed!(random_evolution_seed)
    end

    # determine nr of chains of four vertices
    nr_chains=length(get_all_chains(spatial_network))

    # count finished quenches
    quench_counter = 0

    # evolve network according to given temperature sequence 
    for i in eachindex(evolution_dict["temperature_vec"])

        nr_attempted_bond_switches = Int(round(nr_chains
        * evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][i]))

        spatial_network, total_energy_vec_new, move_accepted_vec_new = (
            evolve_network!(spatial_network,
                evolution_dict,
                nr_attempted_bond_switches, 
                evolution_dict["temperature_vec"][i];
                print_progress = print_progress,
                print_every_nr_attempted_bond_switches 
                    = print_every_nr_attempted_bond_switches,
                quench_counter = quench_counter,
                print_lock = print_lock))

        # update quench_counter
        if (evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][i] 
            == 50 && evolution_dict["temperature_vec"][i] == 0)
            quench_counter += 1
        end

        # concatenate new vectors to previous ones 
        total_energy_vec = vcat(total_energy_vec, total_energy_vec_new)
        move_accepted_vec = vcat(move_accepted_vec, move_accepted_vec_new)

        # print progress if desired
        lock(print_lock) do
            if print_progress
                println("T="
                    *string(evolution_dict["temperature_vec"][i])*" done")
            end
        end

        # if desired, save network
        if save_network_after_each_temperature

            # save evolution data to evolution dict
            evolution_dict["total_energy_vec"] = total_energy_vec
            evolution_dict["move_accepted_vec"] = move_accepted_vec

            # save network and evolution_dict
            save_spatial_network_to_gml(
                spatial_network,
                filename*"_"*string(i);
                evolution_dict = evolution_dict,
                save_path = save_path)

        end

    end

    return [spatial_network, total_energy_vec, move_accepted_vec]
end


"""
Thermally excite entire network
"""
function excite_entire_network!(
    spatial_network::MetaGraphsNext.MetaGraph,
    evolution_dict::Dict;
    relax_first::Bool = false,
    update_total_energy::Bool = true)

    # create tuple containing all vertices
    all_vertices = Tuple(collect(1:spatial_network[]["nr_vertices"]))

    # get cluster for entire network
    cluster_dict = get_cluster_in_shells_dict(spatial_network, 
                                    all_vertices; 
                                    shell_nr = 0)
    
    # if total energy is supposed to be updated, make sure that it is
    # initially up to date
    if update_total_energy && !spatial_network[]["total_energy_up_to_date"]
        spatial_network[]["total_energy"] = cluster_dict["cluster_energy"]
        spatial_network[]["total_energy_up_to_date"] = true
    end

    # if desired, relax network first
    if relax_first
        spatial_network = relax_network_keating!(
            spatial_network,
            get_random_chain(spatial_network),
            evolution_dict;
            threshold_total_energy = Inf,
            update_total_energy = false)
    end

    # excite network
    spatial_network, cluster_dict, cluster_excitation_weight = excite_cluster!(
        spatial_network,
        cluster_dict,
        evolution_dict["temperature_vec"][1];
        update_cluster_energy = update_total_energy,
        update_total_energy = update_total_energy)

    return spatial_network
end


"""
Randomly displace all vertices in a network by choosing a random direction and
sampling a displacement length from a Gaussian distribution with given sigma
"""
function randomly_displace_all_vertices!(
    spatial_network::MetaGraphsNext.MetaGraph;
    sigma::Float64 = 0.1,
    update_total_energy::Bool = true)

    # loop through all vertices
    for vertex in MetaGraphsNext.labels(spatial_network)

        # get random direction using sphere point picking
        theta = 2 * π * rand()
        phi = acos(2 * rand() - 1)
        random_direction = [
            sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi)]

        # get random displacement length
        displacement_length = sigma * randn()

        # get translation vector
        translation_vector = displacement_length .* random_direction

        # move vertex
        spatial_network = move_vertex!(spatial_network, vertex, 
            translation_vector; update_total_energy = false)
    end

    # update total energy if desired
    if update_total_energy
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
        spatial_network[]["total_energy_up_to_date"] = true
    else
        spatial_network[]["total_energy_up_to_date"] = false
    end

    return spatial_network
end
