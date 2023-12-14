"""
Functions for graph plotting
"""


"""
cut all bonds that reach out of the supercell and replace 
them by half way bonds
"""
function cut_bonds_out_of_supercell!(plot_dict::Dict)
    
    #get vector of all bonds in network
    bond_vec = collect(MetaGraphsNext.edge_labels(plot_dict["spatial_network"]))

    #count current vortex
    vortex_count = copy(plot_dict["nr_atoms"])

    #loop through all bonds
    for bond in bond_vec

        #check if bond crosses supercell edge due to periodic boundary
        #conditions
        if (LinearAlgebra.norm(plot_dict["spatial_network"][bond[1]]["position"] 
                .- plot_dict["spatial_network"][bond[2]]["position"]) 
            > plot_dict["supercell_edge_length"]/2)
            
            #determine half way vector
            new_vector = (1/2) .* plot_dict["spatial_network"][bond...]["vector"]

            #add two new vertices and bonds half way of original bond
            for i in 1:2
                plot_dict["spatial_network"][vortex_count + i] = (
                    Dict("position" => (plot_dict["spatial_network"][bond[i]]["position"] 
                                    .+ (-1)^(i+1) .* new_vector )) )

                plot_dict["spatial_network"][bond[i], vortex_count + i] = (
                    Dict("vector" => (-1)^(i+1) .* new_vector, 
                        "distance_squared" => (
                    (1/4) .* plot_dict["spatial_network"][bond...]["distance_squared"] )) )

            end

            #update votex count
            vortex_count += 2

            #cut original bond
            MetaGraphsNext.rem_edge!(plot_dict["spatial_network"],
            bond...)
        end
    end

    plot_dict["nr_atoms"] = vortex_count

    return plot_dict
end


"""
get all atom positions in a vector of tuples
"""
function get_atom_position_vec(plot_dict::Dict)

    #initialize vector of atomic positions
    atom_position_vec = Vector{Tuple}(undef, plot_dict["nr_atoms"])

    #loop through all atoms
    for atom in MetaGraphsNext.labels(plot_dict["spatial_network"])

        #save atomic position as a tuple
        atom_position_vec[atom] = Tuple(
                    plot_dict["spatial_network"][atom]["position"])

    end

    return atom_position_vec
end


"""
Get color vectors for nodes and edges
"""
function get_node_edge_color_vecs(plot_dict::Dict,
    highlight_nodes::Tuple = (),
    highlight_edges::Vector = [] )

    #get node color vector
    node_color_vec = [GLMakie.Colors.colorant"black" for i in 1:plot_dict["nr_atoms"]]

    for highlight_node in highlight_nodes
        node_color_vec[highlight_node] = GLMakie.Colors.colorant"red"
    end

    #get edge color vector

    #get all edges
    edge_vec = collect(MetaGraphsNext.edge_labels(plot_dict["spatial_network"]))

    #initialize edge color vector
    edge_color_vec = Vector{GLMakie.Colors.RGB{GLMakie.Colors.FixedPointNumbers.N0f8}}(
                                                undef, length(edge_vec)) 

    #count current edge
    edge_count = 1

    for edge in edge_vec

        if edge in highlight_edges
            edge_color_vec[edge_count] = GLMakie.Colors.colorant"red"
        else
            edge_color_vec[edge_count] = GLMakie.Colors.colorant"black"
        end

        edge_count += 1
    end

    return [node_color_vec, edge_color_vec]
end



"""
Plot a network in 2d or 3d
"""
function plot_network(graph_dict::Dict;
    highlight_nodes::Tuple = (),
    highlight_edges::Vector = [])

    #get original nr of atoms
    nr_atoms = graph_dict["nr_atoms"]

    #create graph dict to plot
    plot_dict = deepcopy(graph_dict)
    
    #cut all bonds that reach out of supercell and replace
    #them by half way bonds
    plot_dict = cut_bonds_out_of_supercell!(plot_dict)

    #get nr of new virtual vertices
    nr_virtual_vertices = plot_dict["nr_atoms"] - nr_atoms

    #set node size to zero for all virtual vertices
    node_size_vec = vcat(10 .* ones(Int64, nr_atoms), 
                        zeros(Int64, nr_virtual_vertices)) 

    #Get color vectors for nodes and edges
    node_color_vec, edge_color_vec = get_node_edge_color_vecs(
                                    plot_dict,
                                    highlight_nodes,
                                    highlight_edges )

    #get all atom positions in a vector of tuples
    atom_position_vec = get_atom_position_vec(plot_dict)

    #generate plot
    f, ax, p = GraphMakie.graphplot(plot_dict["spatial_network"],
    layout=GraphMakie.NetworkLayout.Spring(dim=3),
    node_size = node_size_vec,
    node_color = node_color_vec,
    edge_color = edge_color_vec)     #ilabels=repr.(1:graph_dict["nr_atoms"]),
    
    #adjust atomic positions
    fixed_layout(_) = atom_position_vec
    p.layout = fixed_layout

    #hide axes and grid
    #GLMakie.hidedecorations!(GLMakie.Axis(f)) 
    #GLMakie.hidespines!(GLMakie.Axis(f))

    return f
end


function get_rainbow_color_vecs(plot_dict )

    #get node color vector
    node_color_vec = [GLMakie.Colors.HSV(rand(1:360), rand(1:360), rand(1:360)) for i in 1:plot_dict["nr_atoms"]]


    #get edge color vector

    #get all edges
    edge_vec = collect(MetaGraphsNext.edge_labels(plot_dict["spatial_network"]))

    #initialize edge color vector
    edge_color_vec = Vector{Tuple{GLMakie.Colors.RGB{GLMakie.Colors.FixedPointNumbers.N0f8}, Float64}}(
                                                undef, length(edge_vec)) 

    #count current edge
    edge_count = 1

    for edge in edge_vec
        edge_color_vec[edge_count] = (GLMakie.Colors.HSV(rand(1:360), rand(1:360), rand(1:360)), 0.6)
        
        edge_count += 1
    end


    return [node_color_vec, edge_color_vec]
end


"""
Plot a network in 2d or 3d
"""
function plot_network_rainbow(graph_dict::Dict)

    #get original nr of atoms
    nr_atoms = graph_dict["nr_atoms"]

    #create graph dict to plot
    plot_dict = deepcopy(graph_dict)
    
    #cut all bonds that reach out of supercell and replace
    #them by half way bonds
    plot_dict = cut_bonds_out_of_supercell!(plot_dict)

    #get nr of new virtual vertices
    nr_virtual_vertices = plot_dict["nr_atoms"] - nr_atoms

    #set node size to zero for all virtual vertices
    node_size_vec = vcat(10 .* ones(Int64, nr_atoms), 
                        zeros(Int64, nr_virtual_vertices)) 

    #get all atom positions in a vector of tuples
    atom_position_vec = get_atom_position_vec(plot_dict)

    #Get color vectors for nodes and edges
    node_color_vec, edge_color_vec = get_rainbow_color_vecs(
                                    plot_dict )

    #get all atom positions in a vector of tuples
    atom_position_vec = get_atom_position_vec(plot_dict)

    GLMakie.set_theme!(GLMakie.theme_black())

    #generate plot
    f, ax, p = GraphMakie.graphplot(plot_dict["spatial_network"],
    layout=GraphMakie.NetworkLayout.Spring(dim=3),
    node_size = 0,
    edge_width = 4,
    node_color = node_color_vec,
    edge_color = edge_color_vec)  

    #adjust atomic positions
    fixed_layout(_) = atom_position_vec
    p.layout = fixed_layout

    #hide axes and grid
    #GLMakie.hidedecorations!(GLMakie.Axis(f)) 
    #GLMakie.hidespines!(GLMakie.Axis(f))

    return f
end