"""
Functions for loading and saving graph data
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