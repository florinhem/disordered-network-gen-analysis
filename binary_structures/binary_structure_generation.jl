"""
These functions generate 3d binary structures
"""


"""
get nodal equation for given surface taken from
10.1016/S0009-2614(00)01418-4
diamond
"""
function nodal_eqn_d(x,y,z)

    return (cos(x)*cos(y)*cos(z) 
            - sin(x)*sin(y)*sin(z) )

end


"""
get nodal equation for given surface taken from
10.1016/S0009-2614(00)01418-4
gyroid
"""
function nodal_eqn_g(x,y,z)

    return (sin(x)*cos(y)
            + sin(y)*cos(z)
            + sin(z)*cos(x) )

end


"""
get nodal equation for given surface taken from
10.1016/S0009-2614(00)01418-4
simple cubic
"""
function nodal_eqn_p(x,y,z)

    return (cos(x) + cos(y) + cos(z) )

end


"""
get nodal equation for given surface taken from
10.1016/S0009-2614(00)01418-4
I-WP surface which has BCC symmetry
"""
function nodal_eqn_i_wp(x,y,z)

    return ( 2 * (cos(x)*cos(y)
                + cos(y)*cos(z)
                + cos(z)*cos(x) )
            - (cos(2*x) + cos(2*y) + cos(2*z) ) )

end


"""
get desired nodal equation.
the following surfact types are supported:
- D 
- G 
- P 
- I-WP
"""
function get_nodal_eqn(surface_type::String)
    
    if surface_type == "D"
        nodal_eqn_fct = nodal_eqn_d

    elseif surface_type  == "G"
        nodal_eqn_fct = nodal_eqn_g
    
    elseif surface_type  == "P"
        nodal_eqn_fct = nodal_eqn_p

    elseif surface_type  == "I-WP"
        nodal_eqn_fct = nodal_eqn_i_wp

    else
        @error "surface type "*surface_type*" is not supported"
    end

    return nodal_eqn_fct

end


"""
obtain binary 3d data array for a given surface.
The unit cell size and voxel_edge_length are given in nanometers.
"""
function get_binary_data_from_nodal_eqn(unit_cell_length=500, 
    nr_unit_cells=10,
    surface_type::String="D";
    voxel_edge_length=10,
    volume_fraction_parameter = 0,
    save_result::Bool=false, 
    save_path=raw"..\structures\nodal_surfaces\\"*surface_type*"_surface_structure.h5")

    # get desired nodal equation
    nodal_eqn = get_nodal_eqn(surface_type)

    # determine the unit cell length in units of voxels
    unit_cell_length_in_voxels = unit_cell_length / voxel_edge_length

    # determine the data edge length in units of voxels
    data_edge_length = Int( round( unit_cell_length_in_voxels * nr_unit_cells  ) )

    # generate array of zeros where data will be stored in
    data_binary = zeros(Bool, data_edge_length, data_edge_length, data_edge_length)

    # determine the wavenumber
    wavenumber = 2*pi/unit_cell_length_in_voxels

    # loop through data array and check nodal equation
    for i in 1:data_edge_length
        for j in 1:data_edge_length
            for k in 1:data_edge_length

                if ( nodal_eqn(wavenumber*i,
                                wavenumber*j,
                                wavenumber*k) < volume_fraction_parameter )

                    data_binary[i,j,k] = 1
                end

            end
        end
    end

    # get essential information about the structure data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = get_data_essentials(
        data_binary)

    # save everything in dictionary
    structure_dict = Dict("data_binary" => data_binary, 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => size_data, 
                            "mean_edge_length_data" => mean_edge_length_data, 
                            "nr_dimensions_data" => nr_dimensions_data,
                            "voxel_edge_length" => voxel_edge_length ,
                            "label" => surface_type,
                            "unit_cell_length" => unit_cell_length,
                            "nr_unit_cells" => nr_unit_cells,
                            "volume_fraction_parameter" => volume_fraction_parameter )

    # if desired, save corrected data
    if save_result
        GU.save_dict_to_h5(structure_dict, save_path)

    end

    return structure_dict

end


"""
Get binary data for a spatial network where the bonds are represented by
a single line of voxels without a 'finite' bond radius
"""
function get_binary_data_from_spatial_network_bonds_only(graph_dict::Dict;
    voxel_edge_length::Float64 = 0.1)

    # generate array of zeros where data will be stored in
    data_binary = zeros(Bool, Int(round(graph_dict["supercell_edge_length"] / voxel_edge_length)), 
                            Int(round(graph_dict["supercell_edge_length"] / voxel_edge_length)), 
                            Int(round(graph_dict["supercell_edge_length"] / voxel_edge_length)))

    # loop through bonds in spatial network and set those voxels to 1 that lie
    # closer to the bond than the bond radius
    for bond in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        # get the bond vector
        direction_vec = graph_dict["spatial_network"][bond...]["vector"]

        # get the bond length
        bond_length = sqrt(graph_dict["spatial_network"][bond...]["distance_squared"])

        # get the number of voxels that are needed to represent the bond
        nr_voxels = Int( round( bond_length / voxel_edge_length ) )

        # get the voxel vector
        voxel_vector = direction_vec ./ nr_voxels

        # loop through voxels and set them to 1
        for i in 1:nr_voxels

            # get the position of the voxel by accounting for periodic boundary conditions
            voxel_position = (graph_dict["spatial_network"][bond[1]]["position"] 
                            .+ i * voxel_vector 
                            .+ graph_dict["supercell_edge_length"]) .% graph_dict["supercell_edge_length"]

            # get the index of the voxel
            voxel_index = (Int.( round.( (voxel_position .+ (voxel_edge_length/2) )/ voxel_edge_length ) )
                            .- 1) .% size(data_binary) .+ 1

            # set the voxel to 1
            data_binary[voxel_index[1], voxel_index[2], voxel_index[3]] = 1

        end

    end

    return data_binary

end


function add_volume_to_bonds(data_binary_bonds_only::Array{Bool,3};
    bond_radius::Float64 = 0.35 ,
    voxel_edge_length::Float64 = 0.1)

    # give volume to the bonds by setting all voxels that lie within the bond radius to 1
    data_binary = zeros(Bool, size(data_binary_bonds_only)...)

    for i in 1:size(data_binary_bonds_only)[1]
        for j in 1:size(data_binary_bonds_only)[2]
            for k in 1:size(data_binary_bonds_only)[3]

                if data_binary_bonds_only[i,j,k] == 1

                    # check window around current voxel
                    for l in (-Int( round( bond_radius / voxel_edge_length ) )
                                :Int( round( bond_radius / voxel_edge_length ) ))
                        for m in (-Int( round( bond_radius / voxel_edge_length ) )
                                    :Int( round( bond_radius / voxel_edge_length ) ))
                            for n in (-Int( round( bond_radius / voxel_edge_length ) )
                                        :Int( round( bond_radius / voxel_edge_length ) ))

                                # check if voxel is within bond radius
                                if  sqrt(l^2 + m^2 + n^2) * voxel_edge_length <= bond_radius

                                    # get index of current voxel by accounting for periodic boundary conditions
                                    cartesian_index_pbc = ( [i+l,j+m,k+n] .+ size(data_binary) .- 1 
                                                            ) .% size(data_binary) .+ 1

                                    # set voxel to 1
                                    data_binary[cartesian_index_pbc...] = 1

                                end

                            end
                        end
                    end

                end

            end
        end
    end

    return data_binary

end


"""
Get binary data for a spatial network with a given bond radius
"""
function get_binary_data_from_spatial_network(graph_dict::Dict;
    bond_radius::Float64 = 0.35,
    voxel_edge_length::Float64 = 0.1,
    save_path::String = raw"..\structures\random_networks\binary_structures\\",
    filename::String = "some_structure",
    save_result::Bool=false)

    # get binary data for for only the bonds without a finite bond radius
    data_binary_bonds_only = get_binary_data_from_spatial_network_bonds_only(graph_dict; 
        voxel_edge_length = voxel_edge_length)

    # give volume to the bonds by setting all voxels that lie within the bond radius to 1
    data_binary = add_volume_to_bonds(data_binary_bonds_only; 
        bond_radius = bond_radius, 
        voxel_edge_length = voxel_edge_length)

    # get essential information about the structure data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = get_data_essentials(
        data_binary)

    # save everything in structure dictionary
    structure_dict = Dict("data_binary" => data_binary, 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => size_data, 
                            "mean_edge_length_data" => graph_dict["supercell_edge_length"], 
                            "nr_dimensions_data" => graph_dict["nr_dimensions"],
                            "voxel_edge_length" => voxel_edge_length ,
                            "label" => filename,
                            "coordination_nr" => graph_dict["coordination_nr"],
                            "nr_vertices" => graph_dict["nr_vertices"],
                            "bond_radius" => bond_radius )

    # if desired, save corrected data
    if save_result
        GU.save_dict_to_h5(structure_dict, save_path*filename*"_structure.h5")

    end

    return structure_dict

end