"""
These functions modify the network graphs, 
e. g. by a bond switch
"""



"""
This function performs a bond switch on a graph.
The argument switched_bond is a tuple of two integers
which is the edge type of the MetaGraphsNext package
"""
function switch_bond!(graph_dict::Dict, switched_bond::Tuple{Int64, Int64} )

    #find the other vertex's neighbors that are the closest to the current vertex
    new_bond_vertex_vec = Vector{Int64}(undef, 2)
    vector_to_new_bond_vertex_vec = Vector{Vector{Float64}}(undef, 2)
    distance_to_new_bond_vertex_vec = Vector{Float64}(undef, 2)

    #get vectors of original neighbors
    original_neighbors_vec_vec = [collect(MetaGraphsNext.neighbor_labels(
        graph_dict["spatial_network"], switched_bond[1]) ),
        collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[2]) )]

    for i in 1:2

        #get the vertex position of bond vertex
        vertex_position_vec = graph_dict["spatial_network"][switched_bond[i]]["position"]

        #get the other bond vertex's neighbors excluding 
        #the bond vertex and the bond vertex's neighbors
        considered_new_bond_vertices_vec = setdiff(original_neighbors_vec_vec[3-i], 
                                            switched_bond[i], 
                                            original_neighbors_vec_vec[i])

        #break if there are no possible new bond vertices
        if considered_new_bond_vertices_vec == []
            new_bond_vec = []
            return [graph_dict, new_bond_vec]
        
        #otherwise, pick a random new bond vertex
        else
            new_bond_vertex_vec[i] = rand(considered_new_bond_vertices_vec)
        end

        #determine vector to new bond vertex
        if switched_bond[i] < new_bond_vertex_vec[i]
            vector_to_new_bond_vertex_vec[i] = get_distance_vector_pbc(
                    vertex_position_vec,
                    graph_dict["spatial_network"][new_bond_vertex_vec[i]]["position"],
                    graph_dict["supercell_edge_length"] )
        else
            vector_to_new_bond_vertex_vec[i] = get_distance_vector_pbc(
                    graph_dict["spatial_network"][new_bond_vertex_vec[i]]["position"],
                    vertex_position_vec,
                    graph_dict["supercell_edge_length"] )
        end

        #determine length of vector to new bond vertex
        distance_to_new_bond_vertex_vec[i] = LinearAlgebra.norm(vector_to_new_bond_vertex_vec[i])

    end

    #create vector to save new bonds
    new_bond_vec = Vector{Tuple{Int64, Int64}}(undef, 2)

    #for each bond vertex, break bond to one neighbor and reconnect to
    #random neighbor of the other vertex
    for i in 1:2

        MetaGraphsNext.rem_edge!(graph_dict["spatial_network"],
            switched_bond[i], new_bond_vertex_vec[3-i])

        new_bond_vec[i] = (switched_bond[i], new_bond_vertex_vec[i])

        graph_dict["spatial_network"][new_bond_vec[i]...] = Dict(
            "vector" => vector_to_new_bond_vertex_vec[i], 
            "distance_squared" => distance_to_new_bond_vertex_vec[i]^2 )
    end

    #note, that total energy is not up to date any more
    graph_dict["total_energy_up_to_date"] = false

    return [graph_dict, new_bond_vec]

end


"""
move an vertex and update its energy and the edges to its neighbors
"""
function move_vertex!(graph_dict::Dict, 
                    vertex_to_move::Int64, 
                    translation_vector::Vector{Float64};
                    update_total_energy::Bool = false,
                    total_energy_fct = get_total_energy_keating)

    #update vertex position by taking periodic boundary conditions into account
    initial_position = graph_dict["spatial_network"][vertex_to_move]["position"]
    graph_dict["spatial_network"][vertex_to_move]["position"] = (
                                                    initial_position .+ translation_vector
                                                            ).%graph_dict["supercell_edge_length"]

    #update outgoing edges                                                    
    for neighbor in MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex_to_move)

        original_distance_vector = graph_dict["spatial_network"][vertex_to_move, neighbor]["vector"]

        #determine new vector, where the direction of the vector (from vertex with lower label
        #to vertex with higher label) has to be taken into account
        graph_dict["spatial_network"][vertex_to_move, neighbor]["vector"] = (
                    original_distance_vector .- sign(neighbor - vertex_to_move)*translation_vector )

        graph_dict["spatial_network"][vertex_to_move, neighbor]["distance_squared"] = (
            LinearAlgebra.norm(graph_dict["spatial_network"][vertex_to_move, neighbor]["vector"])^2
                                                                                    )
    end
    
    #update total energy if desired
    if update_total_energy
        graph_dict["total_energy"] = total_energy_fct(graph_dict)
        graph_dict["total_energy_up_to_date"] = true
    else
        graph_dict["total_energy_up_to_date"] = false
    end

    return graph_dict
end


"""
Relax a single vertex by moving it to the energy minimum 
while fixing its neighbors' positions
"""
function relax_single_vertex_keating!(graph_dict::Dict, vertex_to_relax::Int64;
    optimization_method = Optim.Newton(),
    update_total_energy::Bool = false)
    
    #get initial position of vertex to relax 
    initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

    #get matrix of the vertex's neighbors' positions 
    neighbor_positions_mat = get_neighbor_positions_mat(graph_dict, vertex_to_relax)

    #get next to nearest neighbors' positions
    next_neighbor_positions_arr = get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

    #set energy, gradient and hessian for energy minimization
    energy(x) = energy_from_position_keating(x, graph_dict,
                                                neighbor_positions_mat,
                                                next_neighbor_positions_arr )
                                                
    gradient!(gradient, x) = gradient_keating!(gradient, x, graph_dict, 
                                                neighbor_positions_mat,
                                                next_neighbor_positions_arr)
    
    hessian!(hessian, x) = hessian_keating!(hessian, x, graph_dict,
                                                neighbor_positions_mat,
                                                next_neighbor_positions_arr)
    #find energy minimum
    minimizer_result = Optim.optimize(
                                energy, 
                                gradient!, 
                                hessian!,
                                initial_position, 
                                optimization_method)

    #if minimization converged update dictionary
    if Optim.converged(minimizer_result)

        #get relaxed position and local keating energy
        relaxed_position = Optim.minimizer(minimizer_result)

        #calculate translation vector for relaxed vertex
        translation_vector = relaxed_position .- initial_position

        #move vertex 
        graph_dict = move_vertex!(graph_dict, 
                                vertex_to_relax, 
                                translation_vector;
                                update_total_energy = update_total_energy)

    else
        @warn "Using gradient descent for vertex "*string(vertex_to_relax)
        
        #find energy minimum
        minimizer_result = Optim.optimize(
                                energy, 
                                gradient!, 
                                hessian!,
                                initial_position, 
                                Optim.GradientDescent())

        #if minimization converged update dictionary
        if Optim.converged(minimizer_result)

            #get relaxed position and local keating energy
            relaxed_position = Optim.minimizer(minimizer_result)

            #calculate translation vector for relaxed vertex
            translation_vector = relaxed_position .- initial_position

            #move vertex 
            graph_dict = move_vertex!(graph_dict, 
                                    vertex_to_relax, 
                                    translation_vector;
                                    update_total_energy = update_total_energy)
        else
            @error "No minimum found for vertex "*string(vertex_to_relax)
        end

    end

    return graph_dict

end


"""
calculate translation vector to approximate energy minimum of
Keating strain energy as explained in 10.1142/S0217984987000065
"""
function get_approximate_translation_vector_keating(gradient::Vector{Float64}, 
                                        bond_bending_const::Float64;
                                        relaxation_overshoot_factor_r::Real = 1.5,
                                        relaxation_optimization_parameter_l::Real = 1)

    #determine translation vector
    translation_vector = ((- relaxation_overshoot_factor_r/(4 + 5*bond_bending_const 
                    + relaxation_optimization_parameter_l*LinearAlgebra.norm(gradient)))
            .* gradient
        )

    return translation_vector
end


"""
Approximately relax a single vertex by efficiently determining
the approximated coordinate shift. The corresponding methdo is explained in
10.1142/S0217984987000065
"""
function relax_single_vertex_keating_efficiently!(graph_dict::Dict,
    vertex_to_relax::Int64;
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    update_total_energy::Bool = false)

    #get energy gradient at current vertex position
    gradient = gradient_keating_efficient(graph_dict, vertex_to_relax)

    #calculate translation vector to approximate energy minimum
    translation_vector =  get_approximate_translation_vector_keating(gradient, 
        graph_dict["bond_bending_const"];
        relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
        relaxation_optimization_parameter_l = relaxation_optimization_parameter_l)

    #move vertex 
    graph_dict = move_vertex!(graph_dict, 
                            vertex_to_relax, 
                            translation_vector;
                            update_total_energy = update_total_energy)


    return graph_dict

end


"""
Perform one cycle of relaxing a cluster of vertices 
"""
function relax_cluster_one_cycle_keating!(graph_dict::Dict, 
    cluster_dict::Dict;
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    inefficient_optimization_method = Optim.Newton(),
    update_total_energy::Bool = false,
    update_cluster_energy::Bool = true )

    #get initial cluster energy
    if cluster_dict["cluster_energy_up_to_date"]
        initial_cluster_energy = cluster_dict["cluster_energy"]
    else
        initial_cluster_energy = get_cluster_energy(graph_dict, cluster_dict)
    end

    #get total energy if its not up to date but supposed to be updated later
    if update_total_energy && !graph_dict["total_energy_up_to_date"]
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
    end

    #relax each vertex in the given cluster
    for vertex in cluster_dict["cluster_vertices_to_move_vec"]

        #relax efficiently but approximately or exactly but slowly
        if relax_efficiently
            graph_dict = relax_single_vertex_keating_efficiently!(graph_dict,
    vertex;
    relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
    relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
    update_total_energy = false)
        else
            graph_dict = relax_single_vertex_keating!(graph_dict, vertex;
                            optimization_method=inefficient_optimization_method,
                            update_total_energy = false)
        end

    end

    #update cluster energy if desired
    if update_cluster_energy
        cluster_dict["cluster_energy"] = get_cluster_energy(graph_dict, cluster_dict)
        cluster_dict["cluster_energy_up_to_date"] = true
    else
        cluster_dict["cluster_energy_up_to_date"] = false
    end

    #update total energy if desired
    if update_total_energy

        if cluster_dict["cluster_energy_up_to_date"]
            cluster_energy = cluster_dict["cluster_energy"] 
        else
            cluster_energy = get_cluster_energy(graph_dict, cluster_dict)
        end

        graph_dict["total_energy"] = (graph_dict["total_energy"] 
                                + cluster_energy
                                - initial_cluster_energy)

        graph_dict["total_energy_up_to_date"] = true
    else
        graph_dict["total_energy_up_to_date"] = false
    end

    return [graph_dict, cluster_dict] 
    
end


"""
Fully relax a cluster of vertices. The cluster energy will always be updated
"""
function relax_cluster_keating!(graph_dict::Dict,
    cluster_dict::Dict; 
    nr_max_relaxation_cycles::Int64 = 25,
    break_at_relative_cluster_energy_change::Float64 = 0.001,
    reject_during_relaxation_cycle_threshold::Int64 = 5,
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    update_total_energy::Bool = false,
    print_progress::Bool = false)

    #make sure that cluster energy is up to date
    if !cluster_dict["cluster_energy_up_to_date"]
        cluster_dict["cluster_energy"] = get_cluster_energy(graph_dict, cluster_dict)
    end

    #store initial cluster energy if total energy will be updated
    if update_total_energy
        initial_cluster_energy  = cluster_dict["cluster_energy"]
    end

    #make sure that total energy is up to date if it will be updated later
    if !graph_dict["total_energy_up_to_date"] && update_total_energy
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
    end

    #perform the given number of relaxation cycles
    for cycle_nr in 1:nr_max_relaxation_cycles

        #store previous cluster energy
        previous_cluster_energy = cluster_dict["cluster_energy"]

        graph_dict, cluster_dict = relax_cluster_one_cycle_keating!(graph_dict, 
        cluster_dict;
        relax_efficiently = relax_efficiently,
        relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
        relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
        update_total_energy = false,
        update_cluster_energy = true )

        #break if cluster energy changes less than the given threshold
        relative_cluster_energy_change = (
            abs((previous_cluster_energy - cluster_dict["cluster_energy"])
                    /cluster_dict["cluster_energy"]))

        if (relative_cluster_energy_change < break_at_relative_cluster_energy_change 
                && cycle_nr > reject_during_relaxation_cycle_threshold)
            if print_progress
                println("Breaking at cycle nr "*string(cycle_nr))
            end
            break
        end

        #if cycle nr is above the given threshold, check if the relaxation can 
        #be rejected before full relaxation by estimating the final energy
        if cycle_nr > reject_during_relaxation_cycle_threshold

            #to be implemented

        end

    end

    #update total energy if desired
    if update_total_energy
        graph_dict["total_energy"] = (graph_dict["total_energy"] 
                                    + cluster_dict["cluster_energy"]
                                    - initial_cluster_energy)

        graph_dict["total_energy_up_to_date"] = true
    else
        graph_dict["total_energy_up_to_date"] = false
    end
    
    return [graph_dict, cluster_dict]
end


"""
Move each vertex in the cluster according to 
"""
function excite_cluster!(graph_dict::Dict, cluster_dict::Dict,
                        temperature::Real;
                        update_total_energy = false,
                        update_cluster_energy = false)

    #if total energy will be updated, store initial cluster energy and make sure
    #that total energy is up to date
    if update_total_energy
        if !graph_dict["total_energy_up_to_date"]
            graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
        end

        if cluster_dict["cluster_energy_up_to_date"]
            initial_cluster_energy  = cluster_dict["cluster_energy"]
        else
            initial_cluster_energy = get_cluster_energy(graph_dict, cluster_dict)
        end
    end

    #initialize excitation weight
    cluster_excitation_weight = 1

    #initialize dict where vertex displacements will be stored
    vertex_excitation_dict = Dict()

    #loop through all vertices that are allowed to move
    for vertex in cluster_dict["cluster_vertices_to_move_vec"]

        #get hessian matrix of current vertex at relaxed position
        hessian = hessian_keating_efficient(graph_dict, vertex)

        #get vector of sigmas
        sigma_vec = sqrt.(temperature ./ LinearAlgebra.diag(hessian) ) 

        #draw vector of displacements from Gaussian distribution with standard
        #deviations in the sigma vector
        displacement_vec = sigma_vec .* randn(graph_dict["nr_dimensions"])

        #get excitation weight for current vertex
        vertex_excitation_weight = prod(displacement_vec)

        #multiply it to cluster excitation weight
        cluster_excitation_weight *= vertex_excitation_weight

        #store displacement vector
        vertex_excitation_dict[vertex] = displacement_vec
    end

    #move vertices according to the given displacements
    for vertex in cluster_dict["cluster_vertices_to_move_vec"]

        #move vertices according to the previously thermal fluctuations
        graph_dict = move_vertex!(graph_dict, 
                            vertex, 
                            vertex_excitation_dict[vertex];
                            update_total_energy = false)
    end

    #update cluster energy if desired
    if update_cluster_energy
        cluster_dict["cluster_energy"] = get_cluster_energy(graph_dict, cluster_dict)
        cluster_dict["cluster_energy_up_to_date"] = true
    else
        cluster_dict["cluster_energy_up_to_date"] = false
    end

    #if desired update total energy
    if update_total_energy

        #get excited cluster energy
        if cluster_dict["cluster_energy_up_to_date"]
            excited_cluster_energy = cluster_dict["cluster_energy"]
        else
            excited_cluster_energy = get_cluster_energy(graph_dict, cluster_dict)
        end

        #update total energy
        graph_dict["total_energy"] = (graph_dict["total_energy"] 
                                    + excited_cluster_energy
                                    - initial_cluster_energy)

        graph_dict["total_energy_up_to_date"] = true
    else
        graph_dict["total_energy_up_to_date"] = false
    end

    return [graph_dict, cluster_dict, cluster_excitation_weight]
end


"""
Perform a Monte Carlo move without thermal fluctuations
by switching a bond, relaxing the network and then accepting the network
with Metropolis acceptance probability
"""
function monte_carlo_move!(graph_dict::Dict, 
    temperature::Real; 
    switched_bond::Tuple{Int64, Int64} = get_random_bond(graph_dict),
    nr_max_relaxation_cycles::Int64 = 25,
    reject_during_relaxation_cycle_threshold::Int64 = 5,
    break_at_relative_cluster_energy_change::Float64 = 0.001,
    shell_nr::Int64 = 5,
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    thermal_fluctuations::Bool = true,
    print_progress::Bool = false)

    #save original graph dict 
    initial_graph_dict = deepcopy(graph_dict)

    #get initial cluster before bond switch
    initial_cluster_dict = get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_bond; 
                                    shell_nr = shell_nr)

    #make sure that total energy is up to date
    if !graph_dict["total_energy_up_to_date"]
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
        graph_dict["total_energy_up_to_date"] = true
    end

    #initialize cluster weights which will contribute to the acceptance probability
    #when thermal fluctuations are included
    cluster_relaxation_weight = 1
    cluster_excitation_weight = 1

    #if there are thermal fluctuations, relax cluster first and calculate
    #weights of the corresponding shifts
    if thermal_fluctuations

        #deep copy initial cluster, such that it does not get modified
        cluster_dict = deepcopy(initial_cluster_dict)

        #relax cluster
        graph_dict, cluster_dict = relax_cluster_keating!(graph_dict,
            cluster_dict; 
            nr_max_relaxation_cycles = nr_max_relaxation_cycles,
            break_at_relative_cluster_energy_change = break_at_relative_cluster_energy_change,
            reject_during_relaxation_cycle_threshold = reject_during_relaxation_cycle_threshold,
            relax_efficiently = relax_efficiently,
            relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
            relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
            update_total_energy = false,
            print_progress = print_progress)

        #get cluster weight corresponding to the relaxation translations
        cluster_relaxation_weight = get_cluster_fluctuation_weight(initial_graph_dict, 
                                                            graph_dict, 
                                                            cluster_dict,
                                                            temperature)
    end

    #switch bond
    graph_dict, new_bond_vec = switch_bond!(graph_dict, switched_bond )

    #return immediately if bond switch was not possible
    if new_bond_vec == []
        move_accepted = false
        return [graph_dict, move_accepted, new_bond_vec]
    end

    #get cluster after bond switch
    cluster_dict = get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_bond; 
                                    shell_nr = shell_nr)

    #relax cluster around switched bond and only update energy when there won't be
    #thermal fluctuations included afterward
    graph_dict, cluster_dict = relax_cluster_keating!(graph_dict,
        cluster_dict; 
        nr_max_relaxation_cycles = nr_max_relaxation_cycles,
        break_at_relative_cluster_energy_change = break_at_relative_cluster_energy_change,
        reject_during_relaxation_cycle_threshold = reject_during_relaxation_cycle_threshold,
        relax_efficiently = relax_efficiently,
        relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
        relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
        update_total_energy = false)


    #if desired, include thermal fluctuations by randomly shifting all cluster vertices
    if thermal_fluctuations

        #excite cluster with thermal fluctuations and get the corresponding excitation
        #weight
        graph_dict, cluster_dict, cluster_excitation_weight = excite_cluster!(graph_dict,
                                                            cluster_dict,
                                                            temperature;
                                                            update_total_energy = false,
                                                            update_cluster_energy = true)
    end

    #update total energy
    graph_dict["total_energy"] = (graph_dict["total_energy"] 
                            + cluster_dict["cluster_energy"]
                            - initial_cluster_dict["cluster_energy"])

    graph_dict["total_energy_up_to_date"] = true
    
    #set threshold energy for Metropolis acceptance probability
    total_threshold_energy = (initial_graph_dict["total_energy"]
                - temperature 
                * (cluster_excitation_weight/cluster_relaxation_weight) * log(rand()) )

    #if print_progress
    #    println("Initial energy: "*string(initial_graph_dict["total_energy"]))
    #    println("Cluster relaxation weight: "*string(cluster_relaxation_weight))
    #    println("Cluster excitation weight: "*string(cluster_excitation_weight))
    #    println("Threshold energy: "*string(total_threshold_energy))
    #    println("New energy: "*string(graph_dict["total_energy"]))
    #end

    #accept move if total energy is below threshold
    move_accepted = false

    if graph_dict["total_energy"] <= total_threshold_energy
        move_accepted = true
    else
        graph_dict = initial_graph_dict
        new_bond_vec = []
    end

    return [graph_dict, move_accepted, new_bond_vec]
end


"""
Evolve the network with a given number of attempted Monte Carlo moves
"""
function evolve_network(graph_dict::Dict,
    nr_attempted_bond_switches, 
    temperature::Real; 
    nr_max_relaxation_cycles::Int64 = 25,
    reject_during_relaxation_cycle_threshold::Int64 = 5,
    break_at_relative_cluster_energy_change::Float64 = 0.001,
    shell_nr::Int64 = 5,
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    print_progress::Bool = false,
    random_evolution_seed = Nothing,
    thermal_fluctuations::Bool = true)

    #initialize vectors to keep track of network evolution
    declined_bonds = []
    total_energy_vec = []
    move_accepted_vec = []

    #set seed for random evolution if desired
    if random_evolution_seed !== Nothing
        Random.seed!(random_evolution_seed)
    end

    #attempt given number of bond switches
    for i in 1:nr_attempted_bond_switches

        #get random bond that hasn't been declined since the last
        #accepted switch
        switched_bond = get_random_bond(graph_dict; 
                                        declined_bonds = declined_bonds)

        if print_progress
            println("Attempt bond "*string(switched_bond))
        end

        #attempt Monte Carlo move
        graph_dict, move_accepted, new_bond_vec = monte_carlo_move!(
        graph_dict, 
        temperature; 
        switched_bond = switched_bond,
        nr_max_relaxation_cycles = nr_max_relaxation_cycles,
        break_at_relative_cluster_energy_change = break_at_relative_cluster_energy_change,
        reject_during_relaxation_cycle_threshold = reject_during_relaxation_cycle_threshold,
        shell_nr = shell_nr,
        relax_efficiently = relax_efficiently,
        relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
        relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
        thermal_fluctuations = thermal_fluctuations,
        print_progress = print_progress)

        #update declined bond vec
        if move_accepted
            push!(move_accepted_vec, true)
            declined_bonds = []

            if print_progress
                println("Bond "*string(switched_bond)*" accepted. Energy: "
                *string(graph_dict["total_energy"]))
            end

        else
            push!(move_accepted_vec, false)
            push!(declined_bonds, switched_bond)

            if print_progress
                println("Bond "*string(switched_bond)*" declined.")
            end
        end

        #update total energy
        push!(total_energy_vec, graph_dict["total_energy"])

    end

    return [graph_dict, total_energy_vec, move_accepted_vec]

end


"""
Thermally excite entire network
"""
function excite_entire_network!(graph_dict::Dict, temperature::Real;
    relax_first::Bool = false,
    nr_max_relaxation_cycles::Int64 = 25,
    reject_during_relaxation_cycle_threshold::Int64 = 5,
    break_at_relative_cluster_energy_change::Float64 = 0.001,
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    update_total_energy::Bool = true)

    #create tuple containing all vertices
    all_vertices = Tuple(collect(1:graph_dict["nr_vertices"]))

    #get cluster for entire network
    cluster_dict = get_cluster_in_shells_dict(graph_dict, 
                                    all_vertices; 
                                    shell_nr = 0)
    
    #if total energy is supposed to be updated, make sure that it is
    #initially up to date
    if update_total_energy && !graph_dict["total_energy_up_to_date"]
        graph_dict["total_energy"] = cluster_dict["cluster_energy"]
        graph_dict["total_energy_up_to_date"] = true
    end

    #if desired, relax network first
    if relax_first
        graph_dict, cluster_dict = relax_cluster_keating!(graph_dict,
        cluster_dict; 
        nr_max_relaxation_cycles = nr_max_relaxation_cycles,
        break_at_relative_cluster_energy_change = break_at_relative_cluster_energy_change,
        reject_during_relaxation_cycle_threshold = reject_during_relaxation_cycle_threshold,
        relax_efficiently = relax_efficiently,
        relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
        relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
        update_total_energy = false)
    end

    #excite network
    graph_dict, cluster_dict, cluster_excitation_weight = excite_cluster!(graph_dict,
                                cluster_dict,
                                temperature;
                                update_cluster_energy = update_total_energy,
                                update_total_energy = update_total_energy)

    return graph_dict
end