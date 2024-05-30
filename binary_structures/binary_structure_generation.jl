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
