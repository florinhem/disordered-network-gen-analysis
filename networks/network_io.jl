"""
Functions for loading and saving spatial network data
"""

"""
Save the coordinates of start and target position of all edges in a spatial
network to a CSV file
"""
function save_spatial_network_to_csv(
    spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    save_path::String = raw"..\structures\random_networks\\")

    edges = collect(MetaGraphsNext.edge_labels(spatial_network))

    # create an empty array with the following entries
    data_arr = Array{Float64}(
        undef, length(edges), 2*spatial_network[]["nr_dimensions"])

    # save start and end coordinates of edges to array
    edge_count = 1

    for edge in edges
        data_arr[edge_count, 1:spatial_network[]["nr_dimensions"]] = (
            spatial_network[edge[1]]["position"])
        data_arr[edge_count, (spatial_network[]["nr_dimensions"]
                +1):2*spatial_network[]["nr_dimensions"]] = (
            spatial_network[edge[2]]["position"])

        edge_count += 1

    end

    FileIO.save(save_path*filename, DataFrames.DataFrame(data_arr, :auto) )

    return
end


"""
Save spatial network to a GML format file 
"""
function save_spatial_network_to_gml(
    spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    evolution_dict = nothing,
    save_path::String = raw"..\structures\random_networks\\")

    # save evolution dict if passed
    if evolution_dict !== nothing
        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
    end

    # open new file
    open(save_path*filename*".gml", "w") do opened_file
        # write header
        write(opened_file, "graph [ \n")

        # write network properties
        for (key, value) in spatial_network[]
            if value isa Bool
                value = Int(value)
            end

            write(opened_file, Format.format(
                "  {1} {2}\n",
                key,
                value))
        end

        # loop through vertices
        for vertex in MetaGraphsNext.labels(spatial_network)

            # write vertex
            write(opened_file, Format.format(
                "  node [\n    id {1}\n    label \"{2}\"\n    coordination_nr {3}\n    position [ x {4} y {5} z {6} ]\n  ]\n",
                vertex,
                vertex,
                spatial_network[vertex]["coordination_nr"],
                spatial_network[vertex]["position"][1],
                spatial_network[vertex]["position"][2],
                spatial_network[vertex]["position"][3]))
        end

        # loop through edges
        for edge in MetaGraphsNext.edge_labels(spatial_network)

            # write edge
            write(opened_file, Format.format(
                "  edge [\n    label \"{1}\"\n    source {2}\n    target {3}\n    vector [ x {4} y {5} z {6} ]\n    distance_squared {7}\n    ]\n",
            string(edge[1])*" "*string(edge[2]),
            edge[1],
            edge[2], 
            spatial_network[edge...]["vector"][1],
            spatial_network[edge...]["vector"][2],
            spatial_network[edge...]["vector"][3],
            spatial_network[edge...]["distance_squared"]))

        end

        # write footer
        write(opened_file, "]\n")
    end
    return
end 


"""
Load spatial network from a GML format file 
"""
function load_spatial_network_from_gml(spatial_network_path::String)

    # create an empty spatial network where vertex positions and edge vectors
    # will be stored
    spatial_network = MetaGraphsNext.MetaGraph(
        Graphs.Graph(); 
        label_type = Int64,
        vertex_data_type = Dict{String, Any},
        edge_data_type = Dict{String, Any},
        graph_data = Dict{String, Any}() )

    # load gml file to string
    gml_string = read(spatial_network_path, String)

    # extract network data, node and edge strings
    network_data_string = gml_string[1:findfirst("node", gml_string)[end]]
    nodes_string = gml_string[findfirst("node [", gml_string)[end]+1:findfirst(
        "edge [", gml_string)[1]-1]
    edges_string = gml_string[findfirst("edge [", gml_string)[end]+1:findlast(
        "]", gml_string)[end]-1]

    # Regular expression to match network data keys and values
    pattern = r"(\w+)\s+([\w\.e\-\+]+)"

    # Function to parse values as Int, Float64, or leave as string
    function parse_value(value_str)
        if value_str == "0" 
            return false
        elseif value_str == "1"  
            return true
        elseif !isnothing(tryparse(Int, value_str))  
            return parse(Int, value_str)
        elseif !isnothing(tryparse(Float64, value_str))  
            return parse(Float64, value_str)
        else 
            return value_str
        end
    end

    # Extract the matches using the regex and save them to the network data
    # dictionary
    for m in eachmatch(pattern, network_data_string)
        key = m.captures[1]
        value = parse_value(m.captures[2])
        spatial_network[][key] = value
    end

    # split strings into individual nodes and edges
    nodes_string_list = split(nodes_string, "node")
    edges_string_list = split(edges_string, "edge")

    for node_string in nodes_string_list

        # get vertex and position
        # Regular expressions to extract integer and float values
        id_regex = r"id (\d+)"
        coordination_nr_regex = r"coordination_nr (\d+)"
        position_regex = r"x ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?) y ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?) z ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"

        # Extracting id
        id_match = match(id_regex, node_string)
        vertex = parse(Int, id_match.captures[1])

        # Extracting coordination_nr
        coordination_nr_match = match(coordination_nr_regex, node_string)
        coordination_nr = parse(Int, coordination_nr_match.captures[1])

        # Extracting position
        position_match = match(position_regex, node_string)
        x_value = parse(Float64, position_match.captures[1])
        y_value = parse(Float64, position_match.captures[2])
        z_value = parse(Float64, position_match.captures[3])

        # add vertex to spatial network
        spatial_network[vertex] = Dict(
            "coordination_nr" => coordination_nr,
            "position" =>  [x_value, y_value, z_value],
            )

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
        distance_squared_value = parse(
            Float64, distance_squared_match.captures[1])

        # add edge to spatial network
        spatial_network[source_value, target_value] = Dict(
            "vector" => [x_value, y_value, z_value],
            "distance_squared" => distance_squared_value)

    end
    
    return spatial_network
end


"""
Load an excel file containing node and segment data and extract the relevant
information to create a spatial network.
"""
function import_excel_data(filepath::String)
    # Open the Excel file
    xf = XLSX.openxlsx(filepath)

    # Get tables from sheets
    nodes_table = XLSX.gettable(xf["Nodes"]; infer_eltypes=true)
    segments_table = XLSX.gettable(xf["Segments"]; infer_eltypes=true)

    # Convert to DataFrames
    nodes_df = DataFrames.DataFrame(nodes_table)
    segments_df = DataFrames.DataFrame(segments_table)

    # Extract columns from Nodes sheet
    node_ids = nodes_df[!, "Node ID"] .+ 1  # Adjust for 1-based indexing
    x_coords = nodes_df[!, "X Coord"]
    y_coords = nodes_df[!, "Y Coord"]
    z_coords = nodes_df[!, "Z Coord"]

    # Extract columns from Segments sheet
    node_id_1 = segments_df[!, "Node ID #1"] .+ 1  # Adjust for 1-based indexing
    node_id_2 = segments_df[!, "Node ID #2"] .+ 1  # Adjust for 1-based indexing

    return (
        node_ids,
        x_coords,
        y_coords,
        z_coords,
        node_id_1,
        node_id_2
    )
end


"""
Convert the list of node ids, coordinates and coordination numbers as well as 
the list of segments into a spatial network in the original units as imported
"""
function get_unscaled_spatial_network(node_ids,
    x_coords,
    y_coords,
    z_coords,
    node_id_1,
    node_id_2)

    nr_vertices = length(node_ids)
    nr_dimensions = 3
    min_coords = [minimum(x_coords), minimum(y_coords), minimum(z_coords)]
    max_coords = [maximum(x_coords), maximum(y_coords), maximum(z_coords)]
    edge_lengths = max_coords .- min_coords
    supercell_edge_length = maximum(edge_lengths) * 2

    # get the vector to translate all coordinates such that the network is
    # centered in the supercell
    translation_vec = ((supercell_edge_length .- edge_lengths) ./ 2 
        .- min_coords)

    # create an empty network graph where vertex positions and edge vectors
    # will be stored
    spatial_network_unscaled = MetaGraphsNext.MetaGraph(
        Graphs.Graph(); 
        label_type = Int64,
        vertex_data_type = Dict{String, Any},
        edge_data_type = Dict{String, 
            Union{Float64, Vector{Union{Float64}}}},
        graph_data = Dict{String, Any}(
            "nr_vertices" => nr_vertices,
            "nr_dimensions" => nr_dimensions,
            "supercell_edge_length" => supercell_edge_length
        )
    )

    # label each vertex by its code integer and assign it its position vector
    for vertex in node_ids

        position_vector = [
                x_coords[vertex],
                y_coords[vertex],
                z_coords[vertex]
            ] .+ translation_vec

        spatial_network_unscaled[vertex] = Dict{String, Any}(
            "position" => position_vector
        )
    end

    # add edges to the graph
    for (n1, n2) in zip(node_id_1, node_id_2)
        if n1 == n2
            println("Warning: Self-loop detected at node $n1. Skipping this edge.")
            continue
        end
        if n1 < n2
            edge_vector = (spatial_network_unscaled[n2]["position"] 
            .- spatial_network_unscaled[n1]["position"])
        else
            edge_vector = (spatial_network_unscaled[n1]["position"] 
            .- spatial_network_unscaled[n2]["position"])
        end

        distance_squared = LinearAlgebra.dot(edge_vector, edge_vector)

        spatial_network_unscaled[n1, n2] =  Dict(
            "vector" => edge_vector, 
            "distance_squared" => distance_squared
        )
    end

    # loop through vertices again to assign coordination numbers
    for vertex in node_ids
        spatial_network_unscaled[vertex]["coordination_nr"] = (
            length(MetaGraphsNext.neighbor_labels(
                spatial_network_unscaled, vertex)))
    end

    return spatial_network_unscaled
end 


"""
Convert a spatial network in original units into a spatial network in units of
the equilibrium bond length
"""
function get_scaled_spatial_network(
    spatial_network_unscaled::MetaGraphsNext.MetaGraph)

    # calculate the mean bond length to scale everything in units of the
    # equilibrium bond length
    bond_length_vec = sqrt.([spatial_network_unscaled[edge...][
        "distance_squared"] for edge in MetaGraphsNext.edge_labels(
            spatial_network_unscaled)])
    mean_bond_length = Statistics.mean(bond_length_vec)

    # create a new spatial network where vertex positions and edge vectors
    # will be stored
    spatial_network = MetaGraphsNext.MetaGraph(
        Graphs.Graph(); 
        label_type = Int64,
        vertex_data_type = Dict{String, Any},
        edge_data_type = Dict{String, Union{Float64, Vector{Float64}}},
        graph_data = Dict{String, Any}(
            "nr_vertices" => spatial_network_unscaled[]["nr_vertices"],
            "nr_dimensions" => spatial_network_unscaled[]["nr_dimensions"],
            "supercell_edge_length" => (
                spatial_network_unscaled[]["supercell_edge_length"] 
                / mean_bond_length),
            "mean_bond_length_nm" => mean_bond_length
        )
    )

    # label each vertex by its code integer and assign it its position vector
    for vertex in MetaGraphsNext.labels(spatial_network_unscaled)

        position_vector = (spatial_network_unscaled[vertex]["position"] 
            ./ mean_bond_length)

        spatial_network[vertex] = Dict{String, Any}(
            "position" => position_vector,
            "coordination_nr" => spatial_network_unscaled[vertex][
                "coordination_nr"]
        )
    end

    # add edges to the graph
    for edge in MetaGraphsNext.edge_labels(spatial_network_unscaled)

        edge_vector = (spatial_network_unscaled[edge...]["vector"] 
            ./ mean_bond_length)

        distance_squared = (spatial_network_unscaled[edge...]["distance_squared"] 
            / mean_bond_length^2)

        spatial_network[edge...] =  Dict(
            "vector" => edge_vector,
            "distance_squared" => distance_squared
        )
    end

    return spatial_network
end


"""
Load an excel file containing node and segment data and extract the relevant
information to create a spatial network and save it to a gml file. The
equilibrium bond length can be specified if the coordinates in the excel file
are not given in units of the equilibrium bond length.
"""
function load_spatial_network_from_excel(
    load_path::String,
    load_filename::String;
    save_network = true,
    save_path::String = load_path,
    save_filename::String = load_filename)

    node_ids, x_coords, y_coords, z_coords, node_id_1, node_id_2 = (
        import_excel_data(load_path*load_filename*".xlsx"))

    spatial_network_unscaled = get_unscaled_spatial_network(
        node_ids,
        x_coords,
        y_coords,
        z_coords,
        node_id_1,
        node_id_2
    )

    spatial_network = get_scaled_spatial_network(spatial_network_unscaled)

    if save_network
        # save spatial network to gml file
        save_spatial_network_to_gml(
            spatial_network,
            save_filename;
            save_path = save_path)
    end

    return spatial_network
end


"""
For each evolution dict in a list of filenames in one directory generate a new
spatial network with same evolution in another directory
"""
function generate_spatial_networks_from_evolution_dicts_single_thread(
    filenames,
    evolution_dicts_directory_path::String,
    save_path::String;
    print_every_nr_attempted_bond_switches::Int64 = 100,
    print_progress::Bool = false,
    random_evolution_seed::Int64 = -1,
    save_network_after_each_temperature::Bool = false,
    further_evolve_previous_networks::Bool = false,
    print_lock = Threads.ReentrantLock())

    # loop through files
    for filename in filenames

        # check that file is evolution dict and that the corresponding gml file
        # does not exist yet
        if (endswith(filename, "_evolution.h5") 
            && !isfile(save_path*filename[1:end-13]*".gml"))

            # load evolution dict
            evolution_dict = GU.load_h5_dict(
                evolution_dicts_directory_path*filename)

            # remove move_accepted_vec and total_energy_vec
            delete!(evolution_dict, "move_accepted_vec")
            delete!(evolution_dict, "total_energy_vec")

            # print current thread id and filename if desired
            if print_progress
                lock(print_lock) do
                    Format.printfmtln("Thread {1} is evolving file {2}",
                        Threads.threadid(), filename)
                end
            end
            
            # load previous network if desired 
            if further_evolve_previous_networks

                # get all filenames in spatial network path that contain the
                # current filename
                all_previous_spatial_network_filenames = readdir(save_path)
                current_previous_spatial_network_filenames = [
                    current_filename for current_filename 
                        in all_previous_spatial_network_filenames 
                        if (occursin(filename[1:end-13], current_filename) 
                            && endswith(current_filename, ".gml") )]
                
                # get the filename
                regex = r"(\d+)\.gml$"
                spatial_network_numbers = Vector{Int64}()
                for current_previous_spatial_network_filename in current_previous_spatial_network_filenames
                    number_match = match(regex, 
                        current_previous_spatial_network_filename)
                    extracted_number = parse(Int, number_match.captures[1])[1]
                    push!(spatial_network_numbers, extracted_number)
                end
                maximum_spatial_network_number = maximum(
                    spatial_network_numbers)
                current_previous_spatial_network_filename = (filename[1:end-13]
                    *"_"*string(maximum_spatial_network_number))

                # load spatial network and evolution dict
                spatial_network = load_spatial_network_from_gml(
                    save_path*current_previous_spatial_network_filename*".gml")
                original_evolution_dict = deepcopy(evolution_dict)
                evolution_dict = GU.load_h5_dict(
                    save_path*current_previous_spatial_network_filename
                    *"_evolution.h5")

                # remove those entries from evolution dict, that have already
                # been evolved
                evolution_dict["temperature_vec"] = evolution_dict[
                    "temperature_vec"][maximum_spatial_network_number+1:end]
                evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = (
                    evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][
                        maximum_spatial_network_number+1:end])

                total_energy_vec = evolution_dict["total_energy_vec"]
                move_accepted_vec = evolution_dict["move_accepted_vec"]

            # otherwise generate initial spatial network
            else
                spatial_network = get_periodic_network(evolution_dict)

                total_energy_vec = Vector{Float64}(undef, 0)
                move_accepted_vec = Vector{Bool}(undef, 0)
            end

            # evolve spatial network
            spatial_network, total_energy_vec, move_accepted_vec = (
                evolve_network_temperature_sequence!(
                    spatial_network, evolution_dict; 
                    print_progress = print_progress,
                    total_energy_vec = total_energy_vec,
                    move_accepted_vec = move_accepted_vec,
                    print_every_nr_attempted_bond_switches 
                        = print_every_nr_attempted_bond_switches,
                    random_evolution_seed = random_evolution_seed,
                    save_network_after_each_temperature 
                        = save_network_after_each_temperature,
                    filename = filename[1:end-13],
                    save_path = save_path,
                    print_lock = print_lock))

            # save move_accepted_vec and total_energy_vec
            evolution_dict["total_energy_vec"] = total_energy_vec
            evolution_dict["move_accepted_vec"] = move_accepted_vec

            # restore original temperature and nr monte carlo steps per
            # temperature
            if further_evolve_previous_networks
                evolution_dict["temperature_vec"] = (
                    original_evolution_dict["temperature_vec"])
                evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = (
                    original_evolution_dict[
                        "nr_monte_carlo_steps_per_temperature_vec"])
            end

            # save spatial network
            save_spatial_network_to_gml(
                spatial_network,
                filename[1:end-13];
                evolution_dict = evolution_dict,
                save_path = save_path)
        end

    end

    return
end


"""
Get all evolution dicts in one directory and for each generate a new spatial
network with same evolution in another directory. This is done in a 
multi-threaded (parallel) fashion by splitting all filenames into chunks that
are run on different threads
"""
function generate_spatial_networks_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path::String,
    save_path::String;
    print_every_nr_attempted_bond_switches::Int64 = 100,
    print_progress::Bool = false,
    random_evolution_seed::Int64 = -1,
    save_network_after_each_temperature::Bool = false,
    further_evolve_previous_networks::Bool = false,
    print_lock = Threads.ReentrantLock())

    # get all files in directory
    filenames = readdir(evolution_dicts_directory_path)

    # get filenames of all evolution dicts
    filenames_evolution_dicts = filter(filename
        -> endswith(filename, "_evolution.h5"), filenames)

    # split filenames into chunks for multi-threading
    filename_chunks = Iterators.partition(filenames_evolution_dicts, 
        length(filenames_evolution_dicts) ÷ Threads.nthreads())

    # run all filename chunk in parallel in different threads
    map(filename_chunks) do filename_chunk

        Threads.@spawn (
            generate_spatial_networks_from_evolution_dicts_single_thread(
                filename_chunk,
                evolution_dicts_directory_path,
                save_path;
                print_every_nr_attempted_bond_switches 
                    = print_every_nr_attempted_bond_switches,
                print_progress = print_progress,
                random_evolution_seed = random_evolution_seed,
                save_network_after_each_temperature 
                    = save_network_after_each_temperature,
                further_evolve_previous_networks 
                    = further_evolve_previous_networks,
                print_lock = print_lock))
    end
    
    return
end


"""
This is just an intermediate function to evolve multiple networks from
evolution dicts in multiple runs
"""
function (
    generate_spatial_networks_from_evolution_dicts_single_thread_multiple_runs(
        save_path_filename_tuple_chunks,
        evolution_dicts_directory_path::String;
        print_every_nr_attempted_bond_switches::Int64 = 100,
        print_progress::Bool = false,
        random_evolution_seed::Int64 = -1,
        save_network_after_each_temperature::Bool = false,
        further_evolve_previous_networks::Bool = false,
        print_lock = Threads.ReentrantLock()))

    # loop through the vector of save paths and evolution dict filenames and
    # evolve each of them separately
    for save_path_filename_tuple in save_path_filename_tuple_chunks

        generate_spatial_networks_from_evolution_dicts_single_thread(
            [save_path_filename_tuple[2]],
            evolution_dicts_directory_path,
            save_path_filename_tuple[1];
            print_every_nr_attempted_bond_switches 
                = print_every_nr_attempted_bond_switches,
            print_progress = print_progress,
            random_evolution_seed = random_evolution_seed,
            save_network_after_each_temperature 
                = save_network_after_each_temperature,
            further_evolve_previous_networks 
                = further_evolve_previous_networks,
            print_lock = print_lock)
    end

    return
end


"""
Get all evolution dicts in one directory and for each generate the given number
of spatial networks in separate files. This is done in a multi-threaded
(parallel) fashion by splitting all filenames into chunks that are run on
different threads
"""
function (
    generate_spatial_networks_from_evolution_dicts_in_directory_multiple_runs(
        evolution_dicts_directory_path::String,
        save_path::String;
        print_every_nr_attempted_bond_switches::Int64 = 100,
        print_progress::Bool = false,
        random_evolution_seed::Int64 = -1,
        save_network_after_each_temperature::Bool = false,
        further_evolve_previous_networks::Bool = false,
        runs_vec = collect(1:5),
        print_lock = Threads.ReentrantLock()))

    # get all filenames by reading the evolution dict directory
    filenames = readdir(evolution_dicts_directory_path)
    filenames_evolution_dicts = filter(
        filename -> endswith(filename, "_evolution.h5"), filenames)

    # get vector of filenames and save paths for multiple runs
    save_path_all_runs_vec = Vector{String}(undef, 0)
    filenames_evolution_dicts_all_runs_vec = Vector{String}(undef, 0)

    for run in runs_vec
        # only consider those files, that have not been generated yet
        filenames_evolution_dicts_current_run = filter( filename -> !isfile(
            save_path .* "run_" * string(run) * "/" * filename[1:end-13] 
            * ".gml"),
            filenames_evolution_dicts)

        append!(save_path_all_runs_vec, (save_path .* "run_" 
            .* string.( Int.( ones(length(
                filenames_evolution_dicts_current_run)) .* run ) ) .* "/" ) )
        append!(filenames_evolution_dicts_all_runs_vec, 
            filenames_evolution_dicts_current_run)
    end

    # create tuples out of the elements of both vectors
    save_path_filename_tuple_vec = Vector{Tuple{String, String}}(undef, 
        length(save_path_all_runs_vec))

    for i in eachindex(save_path_filename_tuple_vec)
        save_path_filename_tuple_vec[i] = (save_path_all_runs_vec[i],
            filenames_evolution_dicts_all_runs_vec[i] )
    end

    # split filenames and save_paths into chunks for multi-threading
    save_path_filename_tuple_chunks = Iterators.partition(
        save_path_filename_tuple_vec, 
        length(save_path_filename_tuple_vec) ÷ Threads.nthreads())

    # run all filename chunks in parallel in different threads
    map(save_path_filename_tuple_chunks) do save_path_filename_tuple_chunk

        Threads.@spawn generate_spatial_networks_from_evolution_dicts_single_thread_multiple_runs(
            save_path_filename_tuple_chunk,
            evolution_dicts_directory_path;
            print_every_nr_attempted_bond_switches 
                = print_every_nr_attempted_bond_switches,
            print_progress = print_progress,
            random_evolution_seed = random_evolution_seed,
            save_network_after_each_temperature 
                = save_network_after_each_temperature,
            further_evolve_previous_networks 
                = further_evolve_previous_networks,
            print_lock = print_lock)
    end

    return
end    


"""
cut all bonds that reach out of the supercell and replace them by one bond on
each side with the specified length
"""
function cut_bonds_out_of_supercell!(
    spatial_network::MetaGraphsNext.MetaGraph; 
    vector_out_of_supercell_length = 1)
    
    # get vector of all bonds in network
    bond_vec = collect(MetaGraphsNext.edge_labels(spatial_network))

    # count current vertex
    vertex_count = copy(spatial_network[]["nr_vertices"])

    # get highest vertex index
    current_highest_index = maximum(MetaGraphsNext.labels(spatial_network))

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
                spatial_network[current_highest_index + i] = (
                    Dict(
                        "position" => (spatial_network[bond[i]]["position"] 
                        .+ (-1)^(i+1) .* new_vector ),
                        "coordination_nr" => 1
                        ) )

                spatial_network[bond[i], current_highest_index + i] = (
                    Dict("vector" => (-1)^(i+1) .* new_vector, 
                        "distance_squared" => (vector_out_of_supercell_length^2 
                            .* spatial_network[bond...]["distance_squared"] )))

            end

            vertex_count += 2
            current_highest_index += 2

            # cut original bond
            MetaGraphsNext.rem_edge!(spatial_network,
            bond...)
        end
    end

    spatial_network[]["nr_vertices"] = vertex_count

    return spatial_network
end


"""
Create a padding of 0.6 bond length outside of the supercell for all bonds that
are close to the supercell edge. I choose 0.6 times the equilibrium bond length
to account for bonds that are slightly longer than the equilibrium bond length
"""
function duplicate_bonds_close_to_supercell_edge!(
    spatial_network::MetaGraphsNext.MetaGraph;
    padding = 0.6)

    # count current vertex
    vertex_count = copy(spatial_network[]["nr_vertices"])

    # get highest vertex index
    current_highest_index = maximum(MetaGraphsNext.labels(spatial_network))

    initial_bonds = collect(MetaGraphsNext.edge_labels(spatial_network))

    # loop through bonds
    for bond in initial_bonds

        # get bond's start and target positions and its direction vector
        start_pos = spatial_network[bond[1]]["position"]
        target_pos = spatial_network[bond[2]]["position"]

        # Loop over shifts: -1, 0, 1 in each direction
        for dx in -1:1, dy in -1:1, dz in -1:1
            # skip the case where all shifts are zero
            if dx == 0 && dy == 0 && dz == 0
                continue
            end

            # Shift vector
            shift = [dx, dy, dz] .* spatial_network[]["supercell_edge_length"]

            # Check if point is within l of a face corresponding to the shift
            xcond = ((dx == 1 && start_pos[1] < padding) 
                || (dx == -1 && start_pos[1] 
                    > spatial_network[]["supercell_edge_length"] - padding) 
                || (dx == 1 && target_pos[1] < padding) 
                || (dx == -1 && target_pos[1] 
                    > spatial_network[]["supercell_edge_length"] - padding) 
                ||(dx == 0))
            ycond = ((dy == 1 && start_pos[2] < padding) 
                || (dy == -1 && start_pos[2] 
                    > spatial_network[]["supercell_edge_length"] - padding)
                || (dy == 1 && target_pos[2] < padding)
                || (dy == -1 && target_pos[2] 
                    > spatial_network[]["supercell_edge_length"] - padding)
                || (dy == 0))
            zcond = ((dz == 1 && start_pos[3] < padding)
                || (dz == -1 && start_pos[3] 
                    > spatial_network[]["supercell_edge_length"] - padding)
                || (dz == 1 && target_pos[3] < padding)
                || (dz == -1 && target_pos[3] 
                    > spatial_network[]["supercell_edge_length"] - padding)
                || (dz == 0))

            # If the point is within l of a face, create a new vertex
            if xcond && ycond && zcond
                # Calculate new positions for the start and target vertices
                new_start_pos = start_pos .+ shift
                new_target_pos = target_pos .+ shift

                # check if the new positions are already in the network
                not_in_network = true
                for bond in MetaGraphsNext.edge_labels(spatial_network)
                    if (LinearAlgebra.norm(
                            spatial_network[bond[1]]["position"] 
                            .- new_start_pos) < 1e-4
                        && LinearAlgebra.norm(
                            spatial_network[bond[2]]["position"] 
                            .- new_target_pos) < 1e-4)
                        # if so, do not add the new vertices and bond
                        not_in_network = false
                        break
                    end
                end

                # if the new positions are not in the network, we can add them
                # to the network
                if not_in_network
                    spatial_network[current_highest_index + 1] = (
                        Dict("position" => new_start_pos, 
                            "coordination_nr" => 1 ))
                    spatial_network[current_highest_index + 2] = (
                        Dict("position" => new_target_pos, 
                        "coordination_nr" => 1))

                    spatial_network[current_highest_index + 1, 
                        current_highest_index + 2] = (
                        Dict("vector" => (new_target_pos .- new_start_pos), 
                            "distance_squared" => (LinearAlgebra.norm(
                                new_target_pos .- new_start_pos)^2 )))

                    vertex_count += 2
                    current_highest_index += 2
                end
            end
        end
    end

    spatial_network[]["nr_vertices"] = vertex_count

    return spatial_network
end



"""
Modify the spatial network to prepare it for plotting or optical simulations
by cutting all bonds that reach out of the supercell and by duplicating those
bonds that are close to the supercell edge on the other side of the supercell
just outside the supercell edge
"""
function get_spatial_network_for_simulation!(
    spatial_network::MetaGraphsNext.MetaGraph;
    vector_out_of_supercell_length = 1,
    duplicate_bonds_close_to_supercell_edge::Bool = true,
    save_result::Bool = false,
    filename::String = "some_network",
    save_path::String = raw"..\structures\random_networks\\")
    
    # cut all bonds that reach out of supercell and replace
    # them by new bonds to duplicated vertices outside of the supercell
    spatial_network = cut_bonds_out_of_supercell!(spatial_network; 
        vector_out_of_supercell_length = vector_out_of_supercell_length)

    # duplicate bonds close to supercell edge if desired
    if duplicate_bonds_close_to_supercell_edge
        spatial_network = duplicate_bonds_close_to_supercell_edge!(
            spatial_network)
    end

    # save spatial network to gml format
    if save_result
        save_spatial_network_to_gml(spatial_network, filename*"_for_sim"; 
            save_path=save_path)
    end

    return spatial_network
end


"""
Get mesh from network
"""
function save_mesh_from_spatial_network(
    spatial_network::MetaGraphsNext.MetaGraph, 
    filename::String;
    bond_radius::Float64 = 0.3131,
    vector_out_of_supercell_length = 1,
    save_path::String = raw"..\structures\random_networks\\",
    duplicate_bonds_close_to_supercell_edge::Bool = true,
    format::String = "obj")

    # cut all bonds that reach out of the supercell and by duplicate those 
    # bonds that are close to the supercell edge on the other side of the
    # supercell just outside the supercell edge
    spatial_network = cut_bonds_out_of_supercell!(
        spatial_network; 
        vector_out_of_supercell_length = 1)

    # duplicate bonds close to supercell edge if desired
    if duplicate_bonds_close_to_supercell_edge
        spatial_network = duplicate_bonds_close_to_supercell_edge!(
            spatial_network)
    end

    # loop through bonds
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # get bond's start and target positions and its direction vector
        start_pos = spatial_network[bond[1]]["position"]
        target_pos = spatial_network[bond[2]]["position"]
        # direction_vec = spatial_network[bond...]["vector"]

        # create cylinder object
        current_cylinder = GeometryBasics.Cylinder(
            GeometryBasics.Point( start_pos...),
            GeometryBasics.Point( target_pos...),
            bond_radius)
        
        # mesh cylinder object
        current_cylinder_mesh = GeometryBasics.mesh(current_cylinder)

        # save mesh
        total_path = (save_path*filename*"\\"*string(bond[1])*"_"
            *string(bond[2])*"."*format)

        FileIO.save(total_path, current_cylinder_mesh)

    end

    return
end
