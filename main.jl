
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import MetaGraphsNext
import Graphs

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007

# julia --threads 20

# load some network
dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_4\216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched.gml"

# create an empty network graph where vertexic positions and edge vectors will be stored
spatial_network = MetaGraphsNext.MetaGraph(Graphs.Graph(); 
label_type = Int64,
vertex_data_type = Dict{String, Any},
edge_data_type = Dict{String, Any} )

# load gml file to string
gml_string = read(dict_path, String)

# extract node and edge strings
nodes_string = gml_string[findfirst("node", gml_string)[end]+1:findfirst("edge", gml_string)[1]-1]
edges_string = gml_string[findfirst("edge", gml_string)[end]+1:findlast("]", gml_string)[end]-1]

# split strings into individual nodes and edges
nodes_string_list = split(nodes_string, "node")
edges_string_list = split(edges_string, "edge")

for node_string in nodes_string_list

    # get vertex and position
    # Regular expressions to extract integer and float values
    id_regex = r"id (\d+)"
    position_regex = r"x ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?) y ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?) z ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"

    # Extracting id
    id_match = match(id_regex, node_string)
    vertex = parse(Int, id_match.captures[1])

    # Extracting position
    position_match = match(position_regex, node_string)
    x_value = parse(Float64, position_match.captures[1])
    y_value = parse(Float64, position_match.captures[2])
    z_value = parse(Float64, position_match.captures[3])

    # add vertex to graph
    spatial_network[vertex] = Dict("position" =>  [x_value, y_value, z_value])

end

for edge_string in edges_string_list

    # get source, target, vector and distance squared
    # Regular expressions to extract integer and float values
    source_regex = r"source (\d+)"
    target_regex = r"target (\d+)"
    vector_regex = r"x ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?) y ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?) z ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
    distance_squared_regex = r"distance_squared (\d+\.\d+)"

    # Extracting source
    source_match = match(source_regex, edge_string)
    source_value = parse(Int, source_match.captures[1])

    # Extracting target
    target_match = match(target_regex, edge_string)
    target_value = parse(Int, target_match.captures[1])

    # Extracting vector
    vector_match = match(vector_regex, edge_string)
    x_value = parse(Float64, vector_match.captures[1])
    y_value = parse(Float64, vector_match.captures[2])
    z_value = parse(Float64, vector_match.captures[3])

    # Extracting distance squared
    distance_squared_match = match(distance_squared_regex, edge_string)
    distance_squared_value = parse(Float64, distance_squared_match.captures[1])

    println("source: ", source_value, " target: ", target_value)

    # add edge to graph
    spatial_network[source_value, target_value] = Dict("vector" => [x_value, y_value, z_value],
        "distance_squared" => distance_squared_value)

end