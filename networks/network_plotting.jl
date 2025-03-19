"""
Functions for spatial network plotting
"""

"""
cut all bonds that reach out of the supercell and replace  them by half way
bonds
"""
function cut_bonds_out_of_supercell!(
    spatial_network::MetaGraphsNext.MetaGraph; 
    vector_out_of_supercell_length = 1/2)
    
    # get vector of all bonds in network
    bond_vec = collect(MetaGraphsNext.edge_labels(spatial_network))

    # count current vortex
    vortex_count = copy(spatial_network[]["nr_vertices"])

    # loop through all bonds
    for bond in bond_vec

        # check if bond crosses supercell edge due to periodic boundary
        # conditions
        if (LinearAlgebra.norm(spatial_network[bond[1]]["position"] 
                .- spatial_network[bond[2]]["position"]) 
            > spatial_network[]["supercell_edge_length"]/2)
            
            # determine half way vector
            new_vector = (vector_out_of_supercell_length 
                        .* spatial_network[bond...]["vector"])

            # add two new vertices and bonds half way of original bond
            for i in 1:2
                spatial_network[vortex_count + i] = (
                    Dict("position" => (spatial_network[bond[i]]["position"] 
                        .+ (-1)^(i+1) .* new_vector )) )

                spatial_network[bond[i], vortex_count + i] = (
                    Dict("vector" => (-1)^(i+1) .* new_vector, 
                        "distance_squared" => (vector_out_of_supercell_length^2 
                            .* spatial_network[bond...]["distance_squared"] )))

            end

            vortex_count += 2

            # cut original bond
            MetaGraphsNext.rem_edge!(spatial_network, bond...)
        end
    end

    spatial_network[]["nr_vertices"] = vortex_count

    return spatial_network
end


"""
For each bond in the network, if the both its vertices are on the same side of
the supercell but at least one of them lies close to the supercell edge,
duplicate the bond on the other side of the supercell just outside the edge.
This is required when cylinders are assigned to the bonds and it is plotted or
used in an optical simulation.
"""
function duplicate_bonds_close_to_supercell_edge!(
    spatial_network::MetaGraphsNext.MetaGraph;
    bond_radius::Float64 = 0.35)

    # count current vortex
    vortex_count = copy(spatial_network[]["nr_vertices"])

    # loop through bonds
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # get bond's start and target positions and its direction vector
        start_pos = spatial_network[bond[1]]["position"]
        target_pos = spatial_network[bond[2]]["position"]

        # if one of the two vertices is close to the supercell edge but the
        # vertices are not on opposite sides of the supercell, save another
        # cylinder just outside the supercell on the other side
        if ((any(start_pos .< bond_radius ) 
            || any(target_pos .< bond_radius ) 
            || any((spatial_network[]["supercell_edge_length"] .- start_pos) 
                .< bond_radius )
            || any((spatial_network[]["supercell_edge_length"] .- target_pos) 
                .< bond_radius ) )
            && LinearAlgebra.norm(start_pos .- target_pos) 
                < spatial_network[]["supercell_edge_length"]/2
            && all(start_pos .< spatial_network[]["supercell_edge_length"] )
            && all(target_pos .< spatial_network[]["supercell_edge_length"] )
            && all(start_pos .> 0.0 )
            && all(target_pos .> 0.0 ) )

            # check on which side of supercell the additional bond should be
            # added and calculate new start and target positions
            if (any(start_pos .< bond_radius ) 
                || any(target_pos .< bond_radius ) )
                
                new_start_pos = (
                    start_pos .+ spatial_network[]["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) 
                    .|| (target_pos .< bond_radius ) ))

                new_target_pos = (target_pos 
                    .+ spatial_network[]["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) 
                    .|| (target_pos .< bond_radius ) ))
            else
                new_start_pos = (start_pos 
                    .- spatial_network[]["supercell_edge_length"]
                    .* (((spatial_network[]["supercell_edge_length"] 
                        .- start_pos) .< bond_radius )
                    .|| ((spatial_network[]["supercell_edge_length"] 
                        .- target_pos) .< bond_radius )))

                new_target_pos = (target_pos 
                    .- spatial_network[]["supercell_edge_length"]
                    .* (((spatial_network[]["supercell_edge_length"] 
                        .- start_pos) .< bond_radius )
                    .|| ((spatial_network[]["supercell_edge_length"] 
                        .- target_pos) .< bond_radius )))
            end

            # add two new vertices and the bond between them to the spatial
            # network
            spatial_network[vortex_count + 1] = (
                    Dict("position" => new_start_pos) )
            spatial_network[vortex_count + 2] = (
                        Dict("position" => new_target_pos) )

            spatial_network[vortex_count + 1, vortex_count + 2] = (
                Dict("vector" => (new_target_pos .- new_start_pos), 
                    "distance_squared" => (
                LinearAlgebra.norm(new_target_pos .- new_start_pos)^2 )) )

            vortex_count += 2
        end
    end

    spatial_network[]["nr_vertices"] = vortex_count

    return spatial_network
end


"""
get all vertex positions in a vector of tuples
"""
function get_vertex_position_vec(spatial_network::MetaGraphsNext.MetaGraph)

    # initialize vector of vertex positions
    vertex_position_vec = Vector{Tuple}(undef, 
        spatial_network[]["nr_vertices"])

    # loop through all vertices
    for vertex in MetaGraphsNext.labels(spatial_network)

        # save vertex position as a tuple
        vertex_position_vec[vertex] = Tuple(
                    spatial_network[vertex]["position"])

    end

    return vertex_position_vec
end


"""
Get color vectors for nodes and edges
"""
function get_node_edge_color_vecs(
    spatial_network::MetaGraphsNext.MetaGraph,
    highlight_nodes::Tuple = (),
    highlight_edges::Vector = [] )

    # get node color vector
    node_color_vec = [GLMakie.Colors.colorant"black" for i 
        in 1:spatial_network[]["nr_vertices"]]

    for highlight_node in highlight_nodes
        node_color_vec[highlight_node] = GLMakie.Colors.colorant"red"
    end

    # get all edges
    edge_vec = collect(MetaGraphsNext.edge_labels(spatial_network))

    # initialize edge color vector
    edge_color_vec = Vector{GLMakie.Colors.RGB{
        GLMakie.Colors.FixedPointNumbers.N0f8}}(undef, length(edge_vec)) 

    # count current edge
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
Plot a spatial network in 3d
"""
function plot_spatial_network(
    spatial_network::MetaGraphsNext.MetaGraph;
    highlight_nodes::Tuple = (),
    highlight_edges::Vector = [])

    # get original nr of vertices
    nr_vertices = spatial_network[]["nr_vertices"]
    
    # cut all bonds that reach out of supercell and replace
    # them by half way bonds
    spatial_network = cut_bonds_out_of_supercell!(spatial_network)

    # get nr of new virtual vertices
    nr_virtual_vertices = spatial_network[]["nr_vertices"] - nr_vertices

    # set node size to zero for all virtual vertices
    node_size_vec = vcat(10 .* ones(Int64, nr_vertices), 
        zeros(Int64, nr_virtual_vertices)) 

    # Get color vectors for nodes and edges
    node_color_vec, edge_color_vec = get_node_edge_color_vecs(
        spatial_network, highlight_nodes, highlight_edges )

    # get all vertex positions in a vector of tuples
    vertex_position_vec = get_vertex_position_vec(spatial_network)

    # generate plot
    f, ax, p = GraphMakie.graphplot(spatial_network,
        layout=NetworkLayout.Spring(dim=3),
        node_size = node_size_vec,
        node_color = node_color_vec,
        edge_color = edge_color_vec) 
    
    # adjust vertex positions
    fixed_layout(_) = vertex_position_vec
    p.layout = fixed_layout

    # hide axes and grid
    GLMakie.hidedecorations!(GLMakie.Axis(f)) 
    GLMakie.hidespines!(GLMakie.Axis(f))

    return f
end


"""
Plot a spatial network in 3d without the need of GraphMakie.graphplot which
does not always work
"""
function plot_spatial_network_2(
        spatial_network::MetaGraphsNext.MetaGraph,
        color::String = "black")

    # create spatial network to plot
    spatial_network = deepcopy(spatial_network)
    
    # cut all bonds that reach out of supercell and replace
    # them by half way bonds
    spatial_network = cut_bonds_out_of_supercell!(spatial_network)

    # change Plots backend to plotlyjs
    Plots.plotlyjs()

    # create a line plot for each bond
    figure = Plots.plot()
    for bond in MetaGraphsNext.edge_labels(spatial_network)
        pos_1 = spatial_network[bond[1]]["position"]
        pos_2 = spatial_network[bond[2]]["position"]

        len=2.0
        if LinearAlgebra.norm(pos_1 .- pos_2)<len

            Plots.plot!([pos_1[1], pos_2[1]], 
                        [pos_1[2], pos_2[2]], 
                        [pos_1[3], pos_2[3]], 
            type="scatter3d", mode="lines", color=color, showlegend=false)
        else
            println("We ignored the line in the plot from $pos_1 to $pos_2")
        end
        #=
        markersize_1 = 1 #rand(1:2)
        markersize_2 = 2 #rand(1:2)
        color_1 = :black #Colors.RGB(rand(), rand(), rand())
        color_2 = :black #Colors.RGB(rand(), rand(), rand())
        Plots.scatter!([pos_1[1]], [pos_1[2]], [pos_1[3]], 
                       markersize=markersize_1, marker=:circle, color=color_1, showlegend=false)
        Plots.scatter!([pos_2[1]], [pos_2[2]], [pos_2[3]], 
                       markersize=markersize_2, marker=:circle, color=color_2, showlegend=false)
        =#
    end

    Plots.gr()
    return figure
end


function get_rainbow_color_vecs(spatial_network )

    # get node color vector
    node_color_vec = [GLMakie.Colors.HSV(rand(1:360), rand(1:360), rand(1:360)) 
        for i in 1:spatial_network[]["nr_vertices"]]

    # get all edges
    edge_vec = collect(MetaGraphsNext.edge_labels(spatial_network))

    # initialize edge color vector
    edge_color_vec = Vector{Tuple{GLMakie.Colors.RGB{
        GLMakie.Colors.FixedPointNumbers.N0f8}, Float64}}(undef,
        length(edge_vec)) 

    # count current edge
    edge_count = 1

    # TODO edge was never used:
    for edge in edge_vec
        edge_color_vec[edge_count] = (
            GLMakie.Colors.HSV(rand(1:360), rand(1:360), rand(1:360)), 0.6)
        
        edge_count += 1
    end

    return [node_color_vec, edge_color_vec]
end


"""
Plot a network in 2d or 3d
"""
function plot_network_rainbow(spatial_network::MetaGraphsNext.MetaGraph)

    # get original nr of vertices
    nr_vertices = spatial_network[]["nr_vertices"]
    
    # cut all bonds that reach out of supercell and replace
    # them by half way bonds
    spatial_network = cut_bonds_out_of_supercell!(spatial_network)

    # get nr of new virtual vertices
    nr_virtual_vertices = spatial_network[]["nr_vertices"] - nr_vertices

    # set node size to zero for all virtual vertices
    node_size_vec = vcat(10 .* ones(Int64, nr_vertices), 
                        zeros(Int64, nr_virtual_vertices)) 

    # get all vertex positions in a vector of tuples
    vertex_position_vec = get_vertex_position_vec(spatial_network)

    # Get color vectors for nodes and edges
    node_color_vec, edge_color_vec = get_rainbow_color_vecs(
                                    spatial_network )

    # get all vertex positions in a vector of tuples
    vertex_position_vec = get_vertex_position_vec(spatial_network)

    GLMakie.set_theme!(GLMakie.theme_black())

    # generate plot
    f, ax, p = GraphMakie.graphplot(spatial_network,
    layout=GraphMakie.NetworkLayout.Spring(dim=3),
    node_size = 0,
    edge_width = 4,
    node_color = node_color_vec,
    edge_color = edge_color_vec)  

    # adjust vertex positions
    fixed_layout(_) = vertex_position_vec
    p.layout = fixed_layout

    # hide axes and grid
    GLMakie.hidedecorations!(GLMakie.Axis(f)) 
    GLMakie.hidespines!(GLMakie.Axis(f))

    return f
end
