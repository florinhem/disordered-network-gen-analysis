"""
Functions for spatial network plotting
"""

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

    # get edge color vector

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

    # create spatial network to plot
    spatial_network = deepcopy(spatial_network)
    
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
function plot_spatial_network_2(spatial_network::MetaGraphsNext.MetaGraph)

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
        Plots.plot!([pos_1[1], pos_2[1]], 
        [pos_1[2], pos_2[2]], 
        [pos_1[3], pos_2[3]], 
        type="scatter3d", mode="lines", color="black", showlegend=false)
    end

    Plots.gr()
    return figure
end


function get_rainbow_color_vecs(spatial_network )

    # get node color vector
    node_color_vec = [GLMakie.Colors.HSV(rand(1:360), rand(1:360), rand(1:360)) 
        for i in 1:spatial_network[]["nr_vertices"]]


    # get edge color vector

    # get all edges
    edge_vec = collect(MetaGraphsNext.edge_labels(spatial_network))

    # initialize edge color vector
    edge_color_vec = Vector{Tuple{GLMakie.Colors.RGB{
        GLMakie.Colors.FixedPointNumbers.N0f8}, Float64}}(undef,
        length(edge_vec)) 

    # count current edge
    edge_count = 1

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

    # create spatial network to plot
    spatial_network = deepcopy(spatial_network)
    
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
