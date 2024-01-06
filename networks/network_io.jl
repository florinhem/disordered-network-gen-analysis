"""
Functions for loading and saving graph data
"""

"""
Save the coordinates of start and end of all edges in a graph to a CSV file
"""
function save_graph_to_csv(graph_dict::Dict,
    filename::String;
    save_path::String 
        = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\\")

    edges = collect(MetaGraphsNext.edge_labels(graph_dict["spatial_network"]))

    #create an empty array with the following entries
    data_arr = Array{Float64}(undef, length(edges), 2*graph_dict["nr_dimensions"])

    #save start and end coordinates of edges to array
    edge_count = 1

    for edge in edges
        data_arr[edge_count, 1:graph_dict["nr_dimensions"]] = (
            graph_dict["spatial_network"][edge[1]]["position"])
        data_arr[edge_count, graph_dict["nr_dimensions"]+1:2*graph_dict["nr_dimensions"]] = (
            graph_dict["spatial_network"][edge[2]]["position"])

        edge_count += 1

    end

    #save data
    FileIO.save(save_path*filename, DataFrames.DataFrame(data_arr, :auto) )

    return

end


"""
Save spatial network to an MGformat file and the rest of graph_dict and
evolution_dict to an h5 file
"""
function save_graph_to_h5_and_MGformat(graph_dict::Dict,
    filename::String;
    save_path::String 
        = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\\")

    #save graph to MGformat
    MetaGraphsNext.savegraph(save_path*filename*".mg", graph_dict["spatial_network"])

    #remove spatial_network from graph_dict
    delete!(graph_dict, "spatial_network")

    #save graph dict
    GU.save_dict_to_h5(graph_dict; save_path=save_path*filename*".h5")

    return
end


"""
Load graph to an MGformat file and its properties to an h5 dictionary
"""
function load_graph_from_h5_and_MGformat(dict_path_without_format::String)

    #load spatial network in MGformat
    spatial_network = MetaGraphsNext.loadgraph(
            dict_path_without_format*".mg", MetaGraphsNext.MGFormat())

    #load rest of graph dict
    graph_dict = GU.load_h5_dict(dict_path_without_format*".h5")

    #add spatial network key to graph dict
    graph_dict["spatial_network"] = spatial_network

    return graph_dict
end


"""
Get mesh from network
"""
function save_mesh_from_network(graph_dict::Dict, filename::String;
    bond_radius::Real = 0.05,
    path::String = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\\")

    #loop through bonds
    for bond in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        #get bond's start and target positions and its direction vector
        start_pos = graph_dict["spatial_network"][bond[1]]["position"]
        target_pos = graph_dict["spatial_network"][bond[2]]["position"]
        #direction_vec = graph_dict["spatial_network"][bond...]["vector"]

        #create cylinder object
        current_cylinder = GeometryBasics.Cylinder(
            GeometryBasics.Point( start_pos...),
            GeometryBasics.Point( target_pos...),
            bond_radius)
        
        #mesh cylinder object
        current_cylinder_mesh = GeometryBasics.mesh(current_cylinder)

        #save mesh
        total_path = path*filename*"\\"*string(bond[1])*"_"*string(bond[2])*".obj"

        FileIO.save(total_path, current_cylinder_mesh)

    end

    return
end