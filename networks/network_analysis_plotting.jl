"""
Functions to plot network analysis results
"""


"""
Plot heat map of real part, imaginary part and absolute value of spectral 
density by keeping one component of the wavevector fixed. The fixed wavevector 
value is given in units of (1/nm)
"""
function plot_structure_factor_heatmap(
    structure_factor_dict::Dict,
    save_path::String;
    title="Structure factor",
    save_plot = false,
    clims = (0, 0.1 ),
    x_y_lims = nothing,
    wavevector_component_to_fix::Int64 = 3,
    wavevector_value_fixed = 0)

    #discriminate between different wavevector components that are fixed
    if wavevector_component_to_fix == 1

        #set vectors of x and y axes
        wavenumber_vec_x = structure_factor_dict["wavenumber_vec_vec"][2]
        wavenumber_vec_y = structure_factor_dict["wavenumber_vec_vec"][3]

        #find index of fixed wavenumber value
        wavevector_fixed_index = argmin( abs.( 
            structure_factor_dict["wavenumber_vec_vec"][1]
                                    .- wavevector_value_fixed ) )

        structure_factor_2d_array = structure_factor_dict[
            "structure_factor_array"][wavevector_fixed_index,:,:] 
        
        #set labels and title for the plot
        xlabel = Latex.L"k_y / d^{-1}" 
        ylabel = Latex.L"k_z /  d^{-1}"
        title = (title*Format.format(Latex.L", k_x = {:.2f}", 
            structure_factor_dict["wavenumber_vec_vec"][1][
                wavevector_fixed_index] )*" "*Latex.L" d^{-1}")


    elseif wavevector_component_to_fix == 2
        
        #set vectors of x and y axes
        wavenumber_vec_x = structure_factor_dict["wavenumber_vec_vec"][1]
        wavenumber_vec_y = structure_factor_dict["wavenumber_vec_vec"][3]

        #find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
            structure_factor_dict["wavenumber_vec_vec"][2]
            .- wavevector_value_fixed ) )

        structure_factor_2d_array = structure_factor_dict[
            "structure_factor_array"][:,wavevector_fixed_index,:] 

        
        #set labels and title for the plot
        xlabel = Latex.L"k_x / d^{-1}" 
        ylabel = Latex.L"k_z /  d^{-1}"
        title = (title
        *Format.format(Latex.L", k_y = {:.2f}", structure_factor_dict[
            "wavenumber_vec_vec"][2][wavevector_fixed_index] )
            *" "*Latex.L" d^{-1}")

    elseif wavevector_component_to_fix == 3
        
        #set vectors of x and y axes
        wavenumber_vec_x = structure_factor_dict["wavenumber_vec_vec"][1]
        wavenumber_vec_y = structure_factor_dict["wavenumber_vec_vec"][2]

        #find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
            structure_factor_dict["wavenumber_vec_vec"][3] 
            .- wavevector_value_fixed ) )

        structure_factor_2d_array = structure_factor_dict[
            "structure_factor_array"][:,:,wavevector_fixed_index] 

        xlabel = Latex.L"k_x / d^{-1}" 
        ylabel = Latex.L"k_y /  d^{-1}"
        title = (title
        *Format.format(Latex.L", k_z = {:.2f} ", structure_factor_dict[
            "wavenumber_vec_vec"][3][wavevector_fixed_index] )*" "
            *Latex.L" d^{-1}")

    else
        @error ("Wavevector component to fix must 
                be 1, 2 or 3, but is "*string(wavevector_component_to_fix))
    end

    #permute dimensions of spectral density array, such that they match the
    # axes
    structure_factor_2d_permuted_array = permutedims(structure_factor_2d_array)

    #create plot
    structure_factor_plot = Plots.heatmap(wavenumber_vec_x,
        wavenumber_vec_y,
        structure_factor_2d_permuted_array,
        xlabel=xlabel,
        ylabel=ylabel,
        colorbar_title = "\n"*Latex.L" S (\vec{k} ) " ,
        right_margin = 8Plots.mm,
        legend = true, title=title,
        c = :bluesreds,
        clims = clims,
        aspect_ratio = :equal,
        dpi=250, 
        size = 500 .* (1.2 , 1 )) 

    #set x and y lims if desired
    if x_y_lims !== nothing
        Plots.heatmap!(structure_factor_plot, xlims = x_y_lims, 
            ylims = x_y_lims, size = 500 .* (1 , 1 ))
    end

    #if specified by the argument, save the plot
    if  save_plot
        Plots.savefig(structure_factor_plot, 
            save_path*"_structure_factor_2d.png")

    #otherwise display the plot
    else
        Plots.display(structure_factor_plot)
    end

    return
end


"""
Plot heat map of autocovariance fct by keeping one component of the 
wavevector fixed. The sampling_vector_value fixed needs to be given in units of
nm
"""
function plot_autocovariance_fct_heatmap(
    plot_dict::Dict,
    save_path::String;
    save_plot = false,
    clims = nothing,
    x_y_lims = nothing,
    sampling_vector_component_to_fix::Int64 = 3,
    sampling_vector_value_fixed = 0)

    # discriminate between different wavevector components that are fixed
    if sampling_vector_component_to_fix == 1

        # set vectors of x and y axes
        sampling_distance_vec_x = plot_dict["sampling_distance_vec_vec"][2]
        sampling_distance_vec_y = plot_dict["sampling_distance_vec_vec"][3]

        # find index of fixed sampling vector value
        sampling_vector_fixed_index = argmin( abs.( 
                                plot_dict["sampling_distance_vec_vec"][1] 
                            .- sampling_vector_value_fixed ) )

        autocovariance_fct_2d_array = (plot_dict["autocovariance_fct_array"])[
            sampling_vector_fixed_index,:,:] 
        
        # set labels and title for the plot
        xlabel = Latex.L"r_y / d " 
        ylabel = Latex.L"r_z / d"
        title = (Latex.L"r_x = "
        *Format.format(Latex.L"{1:.1f}", plot_dict[
            "sampling_distance_vec_vec"][1][sampling_vector_fixed_index]  ) 
            *" "*Latex.L" d" )


    elseif sampling_vector_component_to_fix == 2
        
        # set vectors of x and y axes
        sampling_distance_vec_x = plot_dict["sampling_distance_vec_vec"][1]
        sampling_distance_vec_y = plot_dict["sampling_distance_vec_vec"][3]

        # find index of fixed sampling vector value
        sampling_vector_fixed_value, sampling_vector_fixed_index = findmin( 
            abs.( plot_dict["sampling_distance_vec_vec"][2] 
            .- sampling_vector_value_fixed ) ) 

        autocovariance_fct_2d_array = (plot_dict["autocovariance_fct_array"]
            )[:,sampling_vector_fixed_index,:] 

        
        # set labels and title for the plot
        xlabel = Latex.L"r_x / d " 
        ylabel = Latex.L"r_z / d"
        title = (Latex.L"r_y = "
        *Format.format(Latex.L"{1:.1f}", plot_dict[
            "sampling_distance_vec_vec"][2][sampling_vector_fixed_index] ) 
            *" "*Latex.L" d" )


    elseif sampling_vector_component_to_fix == 3
        
        # set vectors of x and y axes
        sampling_distance_vec_x = plot_dict["sampling_distance_vec_vec"][1]
        sampling_distance_vec_y = plot_dict["sampling_distance_vec_vec"][2]

        # find index of fixed sampling vector value
        sampling_vector_fixed_value, sampling_vector_fixed_index = findmin( 
            abs.( plot_dict["sampling_distance_vec_vec"][3] 
            .- sampling_vector_value_fixed  ) )

        autocovariance_fct_2d_array = (plot_dict["autocovariance_fct_array"]
            )[:,:,sampling_vector_fixed_index] 
        
        # set labels and title for the plot
        xlabel = Latex.L"r_x / d " 
        ylabel = Latex.L"r_y / d"
        title = (Latex.L"r_z = "  
        *Format.format(Latex.L"{1:.1f}", plot_dict[
            "sampling_distance_vec_vec"][3][sampling_vector_fixed_index] ) 
            *" "*Latex.L" d" )
        
    else
        @error ("Sampling vector component to fix must be 1, 2 or 3, but is "
            *string(sampling_vector_component_to_fix))
    end

    # create plots
    autocovariance_fct_plot = Plots.heatmap(sampling_distance_vec_x ,
        sampling_distance_vec_y ,
        permutedims(autocovariance_fct_2d_array),
        xlabel=xlabel,
        ylabel=ylabel,
        colorbar_title = "\n"*Latex.L"\chi (\vec{r}) ) " ,
        right_margin = 8Plots.mm,
        legend = true, 
        title=title,
        c = :bluesreds,
        aspect_ratio = :equal,
        dpi=250, 
        size = 500 .* (1.2 , 1 ) ) 

    # set clims if desired
    if clims !== nothing
        Plots.heatmap!(autocovariance_fct_plot, clims = clims)
    end

    # set x and y lims if desired
    if x_y_lims !== nothing
        Plots.heatmap!(autocovariance_fct_plot, xlims = x_y_lims,
            ylims = x_y_lims, size = 500 .* (1 , 1 ))
    end

    # if specified by the argument, save the plot
    if  save_plot
        Plots.savefig(autocovariance_fct_plot, save_path
            *"_autocovariance_fct_direction.png")

    # otherwise display the plot
    else
        Plots.display(autocovariance_fct_plot)
    end

    return
end


"""
Plot heat map of real part, imaginary part and absolute value of spectral
density by keeping one component of the wavevector fixed. The fixed wavevector
value is given in units of (1/nm)
"""
function plot_spectral_density_heatmap(
    plot_dict::Dict,
    save_path::String;
    title="Spectral density",
    save_plot = false,
    clims = (0, 0.1 ),
    x_y_lims = nothing,
    wavevector_component_to_fix::Int64 = 3,
    wavevector_value_fixed = 0,
    plot_im_re::Bool = false)

    # discriminate between different wavevector components that are fixed
    if wavevector_component_to_fix == 1

        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["wavenumber_vec_vec"][2]
        wavenumber_vec_y = plot_dict["wavenumber_vec_vec"][3]

        # find index of fixed wavenumber value
        wavevector_fixed_index = argmin( abs.( 
                                        plot_dict["wavenumber_vec_vec"][1]
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][
            wavevector_fixed_index,:,:] 
        
        # set labels and title for the plot
        xlabel = Latex.L"k_y / d^{-1}" 
        ylabel = Latex.L"k_z / d^{-1}"
        title = (title*", "
        * Latex.L"k_x = " 
        *Format.format(Latex.L"{:.2f}", plot_dict["wavenumber_vec_vec"][1][
            wavevector_fixed_index])*" "*Latex.L"d^{-1}")


    elseif wavevector_component_to_fix == 2
        
        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["wavenumber_vec_vec"][1]
        wavenumber_vec_y = plot_dict["wavenumber_vec_vec"][3]

        # find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
            plot_dict["wavenumber_vec_vec"][2] .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict[
            "spectral_density_array"][:,wavevector_fixed_index,:] 

        
        # set labels and title for the plot
        xlabel = Latex.L"k_x / d^{-1}" 
        ylabel = Latex.L"k_z / d^{-1}"
        title = (title*", "
        * Latex.L"k_y = " 
        *Format.format(Latex.L"{:.2f}", plot_dict["wavenumber_vec_vec"][2][
            wavevector_fixed_index] )*" "*Latex.L"d^{-1}")

    elseif wavevector_component_to_fix == 3
        
        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["wavenumber_vec_vec"][1]
        wavenumber_vec_y = plot_dict["wavenumber_vec_vec"][2]

        # find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
            plot_dict["wavenumber_vec_vec"][3] .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict[
            "spectral_density_array"][:,:,wavevector_fixed_index] 

        # set labels and title for the plot
        xlabel = Latex.L"k_x / d^{-1}" 
        ylabel = Latex.L"k_y / d^{-1}"
        title = (title*", "
        * Latex.L"k_z = " *Format.format(Latex.L"{:.2f}", plot_dict[
            "wavenumber_vec_vec"][3][wavevector_fixed_index]  )
            *" "*Latex.L"d^{-1}")

    else
        @error ("Wavevector component to fix must 
                be 1, 2 or 3, but is "*string(wavevector_component_to_fix))
    end

    # permute dimensions of spectral density array, such that they match the
    # axes
    spectral_density_2d_permuted_array = permutedims(spectral_density_2d_array)

    # normalize spectral density array to its maximum aboslute value
    spectral_density_2d_normalized_array = (spectral_density_2d_permuted_array 
        ./ maximum( abs.(spectral_density_2d_permuted_array) ) )

    # create plots
    abs_plot = Plots.heatmap(wavenumber_vec_x,
        wavenumber_vec_y,
        abs.(spectral_density_2d_normalized_array),
        xlabel=xlabel,
        ylabel=ylabel,
        colorbar_title = "\n"*Latex.L" \mathrm{Abs}( \tilde{\chi} (\vec{k}) ) 
            / \mathrm{Abs}( \tilde{\chi} )_\mathrm{max}  " ,
        right_margin = 8Plots.mm,
        legend = true, title=title,
        c = :bluesreds,
        clims = clims,
        aspect_ratio = :equal,
        dpi=250, 
        size = 500 .* (1.2 , 1 )) 

    # create plots for real and imaginary part if desired
    if plot_im_re                            
        re_plot = Plots.heatmap(wavenumber_vec_x,
            wavenumber_vec_y,
            real.(spectral_density_2d_normalized_array),
            xlabel=xlabel,
            ylabel=ylabel,
            colorbar_title = "\n"*Latex.L"\mathrm{Re}( \tilde{\chi} 
                (\vec{k}) ) / \mathrm{Abs}( \tilde{\chi} )_\mathrm{max}" ,
            right_margin = 8Plots.mm,
            legend = true, title=title,
            c = :bluesreds,
            clims = clims,
            aspect_ratio = :equal,
            dpi=250, 
            size = 500 .* (1.2 , 1 ))

        im_plot = Plots.heatmap(wavenumber_vec_x,
            wavenumber_vec_y,
            imag.(spectral_density_2d_normalized_array),
            xlabel=xlabel,
            ylabel=ylabel,
            colorbar_title = "\n"*Latex.L"\mathrm{Im}( \tilde{\chi} 
                (\vec{k}) ) / \mathrm{Abs}( \tilde{\chi} )_\mathrm{max}" ,
            right_margin = 8Plots.mm,
            legend = true, title=title,
            c = :bluesreds,
            aspect_ratio = :equal,
            dpi=250, 
            size = 500 .* (1.2 , 1 ))
    end

    # set x and y lims if desired
    if x_y_lims !== nothing
        Plots.heatmap!(abs_plot, xlims = x_y_lims, ylims = x_y_lims, 
            size = 500 .* (1 , 1 ))

        if plot_im_re
            Plots.heatmap!(re_plot, xlims = x_y_lims, ylims = x_y_lims, 
                size = 500 .* (1 , 1 ))
            Plots.heatmap!(im_plot, xlims = x_y_lims, ylims = x_y_lims, 
                size = 500 .* (1 , 1 ))
        end
    end

    # if specified by the argument, save the plot
    if  save_plot
        Plots.savefig(abs_plot, save_path*"_spectral_density_abs.png")

        if plot_im_re
            Plots.savefig(re_plot, save_path*"_spectral_density_re.png")
            Plots.savefig(im_plot, save_path*"_spectral_density_im.png")
        end

    # otherwise display the plot
    else
        Plots.display(abs_plot)
    end

    return
end


"""
Plot binary structure as a 3D scatter plot
"""
function plot_binary_structure(data_binary::Array{Bool, 3})

    # create a list of cartesian coordinates for the ones
    coords = findall(x -> x, data_binary)

    # convert list of coordinates to a matrix
    coords_vecs = [collect(index.I) for index in coords]
    coordinates_matrix = hcat(coords_vecs...)

    # create a 3D scatter plot with voxels
    scene = GLMakie.meshscatter(coordinates_matrix[1,:],
        coordinates_matrix[2,:],
        coordinates_matrix[3,:],
        marker=GLMakie.FRect3D(GLMakie.Point3f0(-0.5, -0.5, -0.5), 
            GLMakie.Vec3f0(1, 1, 1)), markersize = 1)

    # display the plot
    GLMakie.display(scene)

    return
end
