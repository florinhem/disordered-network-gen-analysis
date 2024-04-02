"""
These functions are used to plot different statistical measures of 3d binary
structures

This code requires the following Packages:
import Measurements    #for handling data with uncertainty and error propagation
import Plots    #for plotting
import LaTeXStrings as Latex #to display latex symbols in plot labels
"""


"""
Plot window volume times local volume fraction variance as a function of window edge length
to determine whether structure is hyperuniform.
If the plot tends to zero for edge length -> infintity, then the structure is hyperuniform.
The argument plot_dict_vec contains dictionaries with the following keys:
- window_edge_length_vec
- local_volume_fract_variance_vec
- label
"""
function plot_volume_fraction_variance_times_window_volume(nr_dimensions_data::Int64,
                                    plot_dict_vec::Vector,
                                    save_path::String;
                                    title="Hyperuniformity test",
                                    save_plot = false)

    #create empty plot
    hyperuniform_plot = Plots.plot(xlabel="window edge length "*Latex.L"L"*" / voxel edge length",
                                    ylabel=Fmt.format(Latex.L"\sigma_\tau^2(L) \cdot L^{:d} / ", 
                                                        nr_dimensions_data)
                                                    *"voxel edge length"
                                                    *Fmt.format(Latex.L"^{:d}", nr_dimensions_data) ,
                                    legend = true, dpi=250, title=title)

    #loop through plot dictionaries to add plots
    for plot_dict in plot_dict_vec

        #get vector of window volumes
        window_volume_vec = plot_dict["window_edge_length_vec"] .^ nr_dimensions_data

        #plot window volume times local volume fraction variance against window edge length
        Plots.plot!(plot_dict["window_edge_length_vec"], 
                    Measurements.value.(plot_dict["local_volume_fract_variance_vec"]) .* window_volume_vec, 
                    ribbon = Measurements.uncertainty.(plot_dict["local_volume_fract_variance_vec"]) .* window_volume_vec,
                    label = plot_dict["label"])

    end

    #if specified by the argument, save the plot
    if save_plot
        Plots.savefig(save_path*"_volume_fraction_variance_times_window_volume.png")

    else
        Plots.display(hyperuniform_plot)
    end

    return

end



"""
Plot local volume fraction variance as a function of window edge length
to determine whether structure is hyperuniform.
If the plot tends to zero for edge length -> infintity, then the structure is hyperuniform.
The argument plot_dict_vec contains dictionaries with the following keys:
- window_edge_length_vec
- local_volume_fract_variance_vec
- voxel_edge_length
- label
"""
function plot_volume_fraction_variance(plot_dict_vec::Vector,
                                    save_path::String;
                                    title="Local volume fraction variance",
                                    save_plot = false
                                    )

    yticks = 10. .^ collect(-6:-1)

    #create empty plot
    volume_fraction_plot = Plots.plot(xlabel= "measuring window diameter " * Latex.L"L / \mathrm{nm}",
                                    ylabel=Latex.L"\sigma_\tau^2(L) " ,
                                    yaxis=:log,  yticks=yticks,
                                    legend = true, dpi=250, title=title)


    #loop through plot dictionaries to add plots
    for plot_dict in plot_dict_vec

        #plot window volume times local volume fraction variance against window edge length
        Plots.plot!(plot_dict["window_edge_length_vec"]*plot_dict["voxel_edge_length"], 
                    Measurements.value.(plot_dict["local_volume_fract_variance_vec"]) , 
                    ribbon = Measurements.uncertainty.(plot_dict["local_volume_fract_variance_vec"]) ,
                    label = plot_dict["label"])

    end

    #if specified by the argument, save the plot
    if save_plot
        Plots.savefig(save_path*"_volume_fraction_variance.png")

    else
        Plots.display(volume_fraction_plot)
    end

    return

end



"""
Plot autocovariance function as a function of the distance r between two points inside
a structure.
The argument plot_dict_vec contains dictionaries with the following keys:
- sampling_distance_vec
- autocovariance_fct_vec
- voxel_edge_length
- label
"""
function plot_autocovariance_fct(plot_dict_vec::Vector,
                                    save_path::String;
                                    title="Autocovariance function",
                                    save_plot = false,
                                    ylims=nothing)

    #create empty plot
    autocovariance_fct_plot = Plots.plot(xlabel= "sampling distance " * Latex.L"r / \mathrm{nm}",
                                    ylabel=Latex.L"\chi(r) " ,
                                    legend = true, dpi=250, title=title)
    
    #set ylims if specified
    if ylims !== nothing
        Plots.plot!(ylims = ylims)
    end

    #loop through plot dictionaries to add plots
    for plot_dict in plot_dict_vec

        #plot window volume times local volume fraction variance against window edge length
        Plots.plot!(plot_dict["sampling_distance_vec"]*plot_dict["voxel_edge_length"], 
                    Measurements.value.(plot_dict["autocovariance_fct_vec"]), 
                    ribbon = Measurements.uncertainty.(plot_dict["autocovariance_fct_vec"]),
                    label = plot_dict["label"])

    end

    #if specified by the argument, save the plot
    if save_plot && ylims == nothing
        Plots.savefig(save_path*"_autocovariance_fct.png")

    elseif save_plot
        Plots.savefig(save_path*"_autocovariance_fct_zoom.png")

    else
        Plots.display(autocovariance_fct_plot)
    end

    return

end



"""
Plot absolute value, real and imaginary part of spectral density to 
determine whether structure is hyperuniform.
If the plot tends to zero for k -> 0, then the structure is hyperuniform.
The argument plot_dict_vec contains dictionaries with the following keys:
- wavenumber_vec
- spectral_density_vec
- voxel_edge_length
- label
"""
function plot_spectral_density(plot_dict_vec::Vector,
                                save_path::String;
                                title="Spectral density",
                                save_plot = false,
                                xlims = nothing)

    #create empty plot
    spectral_density_plot = Plots.plot(xlabel="wavenumber " * Latex.L"k / ( 1/ \mathrm{nm} )",
                                    ylabel=Latex.L"\tilde{\chi} (k) ",
                                    legend = true, dpi=250, title=title)

    #if desired, zoom in along x axes
    if xlims !== nothing
        Plots.plot!(spectral_density_plot, xlims = xlims)
    end

    #loop through plot dictionaries to add plots
    for plot_dict in plot_dict_vec

        #plot spectral density
        Plots.plot!(plot_dict["wavenumber_vec"] ./ plot_dict["voxel_edge_length"], 
                    real.( Measurements.value.(plot_dict["spectral_density_vec"]) ), 
                    ribbon = real.( Measurements.uncertainty.(plot_dict["spectral_density_vec"]) ),
                    label = plot_dict["label"])

    end

    #if specified by the argument, save the plot
    if save_plot && xlims!== nothing
        Plots.savefig(save_path*"_spectral_density_zoom.png")

    elseif save_plot
        Plots.savefig(save_path*"_spectral_density.png")

    #otherwise display the plot
    else
        Plots.display(spectral_density_plot)

    end

    return

end


"""
Determine a reasonable nr of measurements per distance by performing a convergence
analysis of the autocovariance function for r=0 where its supposed to yield
(1-volume_fract_tot)*volume_fract_tot.
The plots suggest to choose nr_measurements_per_distance >~ 10000
"""
function convergence_analysis_autocovariance_fct_nr_measurements_per_distance(
                    structure_dict::Dict;
                    nr_measurements_per_distance_vec = Int.( round.( 10 .^ collect(1:0.1:4) ) ),
                    title="Convergence analysis autocovariance function",
                    save_plot = false,
                    save_path=raw"..\plots\\",
                    save_filename="convergence_analysis_nr_measurements_autocovariance_fct.png"  )

    #in this array the autocovariance fct at distance 0 will be stored as a function of
    #the nr of measurements
    autocovariance_fct_vec = Vector{Measurements.Measurement}(undef, length(nr_measurements_per_distance_vec))

    #loop through nr of measurements
    for i in eachindex(nr_measurements_per_distance_vec)

    autocovariance_fct_vec[i] = get_autocovariance_fct_isotrope(0,
        structure_dict["size_data"],
        structure_dict["volume_fract_tot"],
        structure_dict["data_binary"];
        nr_measurements_per_distance = nr_measurements_per_distance_vec[i])

    end
    
    #initialize plot and set labels and characteristics
    convergence_analysis_plot = Plots.plot(xlabel="# measurements",
                ylabel=Latex.L"\chi(r=0) " ,
                xaxis=:log,
                legend = true, dpi=250, title=title)

    #plot measured autocovariance functions
    Plots.plot!(nr_measurements_per_distance_vec, 
                Measurements.value.(autocovariance_fct_vec),
                ribbon = Measurements.uncertainty.( autocovariance_fct_vec),
                label = "from measurement")

    #determine and plot theoretical values
    theoretical_value = structure_dict["volume_fract_tot"]*(1 - structure_dict["volume_fract_tot"])

    Plots.plot!( [minimum(nr_measurements_per_distance_vec), 
                maximum(nr_measurements_per_distance_vec)], 
                [theoretical_value, theoretical_value],
                label = "from theory")          

    #if specified by the argument, save the plot
    if save_plot
        Plots.savefig(save_path*save_filename)

    else
        Plots.display(convergence_analysis_plot)
    end

    return

end



"""
Plot heat map of autocovariance fct by keeping one component of the 
wavevector fixed.
The sampling_vector_value fixed needs to be given in units of nm
"""
function plot_autocovariance_fct_heatmap(plot_dict::Dict,
    save_path::String;
    save_plot = false,
    clims = nothing,
    x_y_lims = nothing,
    sampling_vector_component_to_fix::Int64 = 3,
    sampling_vector_value_fixed = 0)

    #discriminate between different wavevector components that are fixed
    if sampling_vector_component_to_fix == 1

        #set vectors of x and y axes
        sampling_distance_vec_x = plot_dict["sampling_distance_vec_vec"][2]
        sampling_distance_vec_y = plot_dict["sampling_distance_vec_vec"][3]

        #find index of fixed sampling vector value
        sampling_vector_fixed_index = argmin( abs.( 
                                plot_dict["sampling_distance_vec_vec"][1] 
                            .- (sampling_vector_value_fixed / plot_dict["voxel_edge_length"] ) ) )

        autocovariance_fct_2d_array = Measurements.value.(
                                            plot_dict["autocovariance_fct_array"]
                                        )[sampling_vector_fixed_index,:,:] 

        uncertainty_2d_array = Measurements.uncertainty.(
                                            plot_dict["autocovariance_fct_array"]
                                        )[sampling_vector_fixed_index,:,:] 
        
        #set labels and title for the plot
        xlabel = Latex.L"r_y / \mathrm{nm} " 
        ylabel = Latex.L"r_z / \mathrm{nm}"
        title = (Latex.L"r_x = "
        *Fmt.format(Latex.L"{1:.1f}",
                        (plot_dict["sampling_distance_vec_vec"][1][sampling_vector_fixed_index]
                                        * plot_dict["voxel_edge_length"] ) ) 
        *" "*Latex.L" \mathrm{nm},  \Delta \chi (\vec{r}) ) = "
        *Fmt.format(Latex.L"{1:.3f}",
        Statistics.mean( uncertainty_2d_array) ) )


    elseif sampling_vector_component_to_fix == 2
        
        #set vectors of x and y axes
        sampling_distance_vec_x = plot_dict["sampling_distance_vec_vec"][1]
        sampling_distance_vec_y = plot_dict["sampling_distance_vec_vec"][3]

        #find index of fixed sampling vector value
        sampling_vector_fixed_value, sampling_vector_fixed_index = findmin( abs.( 
                                plot_dict["sampling_distance_vec_vec"][2] 
                            .- (sampling_vector_value_fixed / plot_dict["voxel_edge_length"] ) ) )

        autocovariance_fct_2d_array = Measurements.value.(
                                            plot_dict["autocovariance_fct_array"]
                                        )[:,sampling_vector_fixed_index,:] 

        uncertainty_2d_array = Measurements.uncertainty.(
                                            plot_dict["autocovariance_fct_array"]
                                        )[:,sampling_vector_fixed_index,:] 
        
        #set labels and title for the plot
        xlabel = Latex.L"r_x / \mathrm{nm} " 
        ylabel = Latex.L"r_z / \mathrm{nm}"
        title = (Latex.L"r_y = "
        *Fmt.format(Latex.L"{1:.1f}",
                        (plot_dict["sampling_distance_vec_vec"][2][sampling_vector_fixed_index]
                                        * plot_dict["voxel_edge_length"] ) ) 
        *" "*Latex.L" \mathrm{nm},  \Delta \chi (\vec{r}) ) = "
        *Fmt.format(Latex.L"{1:.3f}",
        Statistics.mean( uncertainty_2d_array) ) )


    elseif sampling_vector_component_to_fix == 3
        
        #set vectors of x and y axes
        sampling_distance_vec_x = plot_dict["sampling_distance_vec_vec"][1]
        sampling_distance_vec_y = plot_dict["sampling_distance_vec_vec"][2]

        #find index of fixed sampling vector value
        sampling_vector_fixed_value, sampling_vector_fixed_index = findmin( abs.( 
                                plot_dict["sampling_distance_vec_vec"][3] 
                            .- (sampling_vector_value_fixed / plot_dict["voxel_edge_length"] ) ) )

        autocovariance_fct_2d_array = Measurements.value.(
                                            plot_dict["autocovariance_fct_array"]
                                        )[:,:,sampling_vector_fixed_index] 

        uncertainty_2d_array = Measurements.uncertainty.(
                                            plot_dict["autocovariance_fct_array"]
                                        )[:,:,sampling_vector_fixed_index] 
        
        #set labels and title for the plot
        xlabel = Latex.L"r_x / \mathrm{nm} " 
        ylabel = Latex.L"r_y / \mathrm{nm}"
        title = (Latex.L"r_z = "  
        *Fmt.format(Latex.L"{1:.1f}",
                        (plot_dict["sampling_distance_vec_vec"][3][sampling_vector_fixed_index]
                                        * plot_dict["voxel_edge_length"] ) ) 
        *" "*Latex.L" \mathrm{nm},  \Delta \chi (\vec{r}) ) = "
        *Fmt.format(Latex.L"{1:.3f}",
        Statistics.mean( uncertainty_2d_array) ) )
        
    else
        @error ("Sampling vector component to fix must 
                be 1, 2 or 3, but is "*string(sampling_vector_component_to_fix))
    end

    #create plots
    autocovariance_fct_plot = Plots.heatmap(sampling_distance_vec_x .* plot_dict["voxel_edge_length"],
                                sampling_distance_vec_y .* plot_dict["voxel_edge_length"],
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
                                size = 500 .* (1.2 , 1 )  ) #600 .* (1, length(sampling_distance_vec_y)/length(sampling_distance_vec_x) )
                                #dpi=250, 

    #set clims if desired
    if clims !== nothing
        Plots.heatmap!(autocovariance_fct_plot, clims = clims)
    end

    #set x and y lims if desired
    if x_y_lims !== nothing
        Plots.heatmap!(autocovariance_fct_plot, xlims = x_y_lims, ylims = x_y_lims, size = 500 .* (1 , 1 ))
    end

    #if specified by the argument, save the plot
    if  save_plot
        Plots.savefig(autocovariance_fct_plot, save_path*"_autocovariance_fct_direction.png")


    #otherwise display the plot
    else
        Plots.display(autocovariance_fct_plot)

    end

    return
end



"""
Plot heat map of real part, imaginary part and absolute value of spectral density
by keeping one component of the wavevector fixed.
The fixed wavevector value is given in units of (1/nm)
"""
function plot_spectral_density_heatmap(plot_dict::Dict,
    save_path::String;
    title="Spectral density",
    save_plot = false,
    clims = (0, 0.1 ),
    x_y_lims = nothing,
    wavevector_component_to_fix::Int64 = 3,
    wavevector_value_fixed = 0)

    #discriminate between different wavevector components that are fixed
    if wavevector_component_to_fix == 1

        #set vectors of x and y axes
        wavenumber_vec_x = plot_dict["wavenumber_vec_vec"][2]
        wavenumber_vec_y = plot_dict["wavenumber_vec_vec"][3]

        #find index of fixed wavenumber value
        wavevector_fixed_index = argmin( abs.( 
                                        plot_dict["wavenumber_vec_vec"][1]
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][wavevector_fixed_index,:,:] 
        
        #set labels and title for the plot
        xlabel = Latex.L"k_y / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_z / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_x = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["wavenumber_vec_vec"][1][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")


    elseif wavevector_component_to_fix == 2
        
        #set vectors of x and y axes
        wavenumber_vec_x = plot_dict["wavenumber_vec_vec"][1]
        wavenumber_vec_y = plot_dict["wavenumber_vec_vec"][3]

        #find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
                                    plot_dict["wavenumber_vec_vec"][2]
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][:,wavevector_fixed_index,:] 

        
        #set labels and title for the plot
        xlabel = Latex.L"k_x / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_z / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_y = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["wavenumber_vec_vec"][2][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")

    elseif wavevector_component_to_fix == 3
        
        #set vectors of x and y axes
        wavenumber_vec_x = plot_dict["wavenumber_vec_vec"][1]
        wavenumber_vec_y = plot_dict["wavenumber_vec_vec"][2]

        #find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
                                    plot_dict["wavenumber_vec_vec"][3] 
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][:,:,wavevector_fixed_index] 

        #set labels and title for the plot
        xlabel = Latex.L"k_x / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_y / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_z = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["wavenumber_vec_vec"][3][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")

    else
        @error ("Wavevector component to fix must 
                be 1, 2 or 3, but is "*string(wavevector_component_to_fix))
    end

    #scale x and y axes in units of 1/nm
    x_axis = wavenumber_vec_x  ./ plot_dict["voxel_edge_length"]
    y_axis = wavenumber_vec_y ./ plot_dict["voxel_edge_length"]

    #permute dimensions of spectral density array, such that they match the axes
    spectral_density_2d_permuted_array = permutedims(spectral_density_2d_array)

    #normalize spectral density array to its maximum aboslute value
    spectral_density_2d_normalized_array = (spectral_density_2d_permuted_array 
                                ./ maximum( abs.(spectral_density_2d_permuted_array) ) )

    #create plots
    abs_plot = Plots.heatmap(x_axis,
                                y_axis,
                                abs.(spectral_density_2d_normalized_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = "\n"*Latex.L" \mathrm{Abs}( \tilde{\chi} (\vec{k}) ) / \mathrm{Abs}( \tilde{\chi} )_\mathrm{max}  " ,
                                right_margin = 8Plots.mm,
                                legend = true, title=title,
                                c = :bluesreds,
                                clims = clims,
                                aspect_ratio = :equal,
                                dpi=250, 
                                size = 500 .* (1.2 , 1 )) #2/3 * length(y_axis)/length(x_axis)

    re_plot = Plots.heatmap(x_axis,
                                y_axis,
                                real.(spectral_density_2d_normalized_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = "\n"*Latex.L"\mathrm{Re}( \tilde{\chi} (\vec{k}) ) / \mathrm{Abs}( \tilde{\chi} )_\mathrm{max}" ,
                                right_margin = 8Plots.mm,
                                legend = true, title=title,
                                c = :bluesreds,
                                clims = clims,
                                aspect_ratio = :equal,
                                dpi=250, 
                                size = 500 .* (1.2 , 1 ))

    im_plot = Plots.heatmap(x_axis,
                                y_axis,
                                imag.(spectral_density_2d_normalized_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = "\n"*Latex.L"\mathrm{Im}( \tilde{\chi} (\vec{k}) ) / \mathrm{Abs}( \tilde{\chi} )_\mathrm{max}" ,
                                right_margin = 8Plots.mm,
                                legend = true, title=title,
                                c = :bluesreds,
                                aspect_ratio = :equal,
                                dpi=250, 
                                size = 500 .* (1.2 , 1 ))

    #set x and y lims if desired
    if x_y_lims !== nothing
        Plots.heatmap!(abs_plot, xlims = x_y_lims, ylims = x_y_lims, size = 500 .* (1 , 1 ))
        Plots.heatmap!(re_plot, xlims = x_y_lims, ylims = x_y_lims, size = 500 .* (1 , 1 ))
        Plots.heatmap!(im_plot, xlims = x_y_lims, ylims = x_y_lims, size = 500 .* (1 , 1 ))
    end

    #if specified by the argument, save the plot
    if  save_plot
        Plots.savefig(abs_plot, save_path*"_spectral_density_abs.png")
        Plots.savefig(re_plot, save_path*"_spectral_density_re.png")
        Plots.savefig(im_plot, save_path*"_spectral_density_im.png")

    #otherwise display the plot
    else
        Plots.display(abs_plot)

    end

    return
end



"""
This function plots several statistical measures, that is
- local volume fraction variance
- autocovariance function
- spectral density
"""
function plot_statistical_measures(data_path_vec,
                                save_path::String;
                                voxel_edge_length_vec=nothing,
                                label_vec=nothing,
                                spectral_density_xlims = nothing,
                                autocovariance_fct_direction_x_y_lims = nothing,
                                autocovariance_fct_direction_clims = nothing,
                                spectral_density_heatmaps_clims = (0, 0.1),
                                spectral_density_heatmaps_x_y_lims = nothing,
                                plot_autocovariance_fct_bool = true,
                                plot_spectral_density_bool = true,
                                plot_local_volume_fraction_variance_bool = true,
                                plot_autocovariance_fct_direction_bool = true,
                                plot_spectral_density_along_directions_bool = true,
                                plot_spectral_density_heatmaps_bool = true
                                )

    #if desired plot autocovariance function
    if plot_autocovariance_fct_bool

        #initialize vector for plot dicts
        autocovariance_fct_plot_dict_vec = []

        #loop through data that will be plotted
        for i in eachindex(data_path_vec)

            #load plot dictionary
            autocovariance_fct_plot_dict = GU.load_h5_dict(data_path_vec[i]*"_autocovariance_fct.h5")

            #if desired adjust label and voxel edge length
            if  label_vec !== nothing
                autocovariance_fct_plot_dict["label"] = label_vec[i]
            end

            #if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                autocovariance_fct_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            #add current plot dict to vector
            push!(autocovariance_fct_plot_dict_vec, autocovariance_fct_plot_dict)

        end

        #plot the data
        plot_autocovariance_fct(autocovariance_fct_plot_dict_vec,
                        save_path;
                        title="Autocovariance function",
                        save_plot = true)

        #plot zoom into the data
        plot_autocovariance_fct(autocovariance_fct_plot_dict_vec,
                        save_path;
                        title="Autocovariance function",
                        save_plot = true,
                        ylims = [-0.07, 0.07])

    end

    #if desired plot spectral density
    if plot_spectral_density_bool

        #initialize vector for plot dicts
        spectral_density_plot_dict_vec = []

        #loop through data that will be plotted
        for i in eachindex(data_path_vec)

            #load plot dictionary
            spectral_density_plot_dict = GU.load_h5_dict(data_path_vec[i]*"_spectral_density.h5")

            #if desired adjust label and voxel edge length
            if label_vec !== nothing
                spectral_density_plot_dict["label"] = label_vec[i]
            end

            #if desired adjust label and voxel edge length
            if  voxel_edge_length_vec !== nothing
                spectral_density_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            #add current plot dict to vector
            push!(spectral_density_plot_dict_vec, spectral_density_plot_dict)

        end

        #plot the data
        plot_spectral_density(spectral_density_plot_dict_vec,
                        save_path;
                        save_plot = true)

        #if specified, plot a zoom into the spectral sensity
        if spectral_density_xlims !== nothing

            plot_spectral_density(spectral_density_plot_dict_vec,
                            save_path;
                            save_plot = true,
                            xlims = spectral_density_xlims)
        end

    end


    #if desired plot local volume fraction variance
    if plot_local_volume_fraction_variance_bool

        #initialize vector for plot dicts
        local_volume_fraction_variance_plot_dict_vec = []

        #loop through data that will be plotted
        for i in eachindex(data_path_vec)

            #load plot dictionary
            local_volume_fraction_variance_plot_dict = GU.load_h5_dict(
                                                            data_path_vec[i]*"_volume_fraction_variance.h5")

            #if desired adjust label and voxel edge length
            if  label_vec !== nothing
                local_volume_fraction_variance_plot_dict["label"] = label_vec[i]
            end

            #if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                local_volume_fraction_variance_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            #add current plot dict to vector
            push!(local_volume_fraction_variance_plot_dict_vec, local_volume_fraction_variance_plot_dict)

        end

        #plot the data
        plot_volume_fraction_variance(local_volume_fraction_variance_plot_dict_vec,
                        save_path;
                        save_plot = true)

    end


    #if desired plot autocovariance function heatmaps
    if plot_autocovariance_fct_direction_bool

        #loop through data that will be plotted
        for i in eachindex(data_path_vec)

            #load plot dictionary
            autocovariance_fct_direction_plot_dict = GU.load_h5_dict(
                                                            data_path_vec[i]*"_autocovariance_fct_direction_complete.h5")

            #if desired adjust label and voxel edge length
            if label_vec  !== nothing
                autocovariance_fct_direction_plot_dict["label"] = label_vec[i]
            end

            #if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                autocovariance_fct_direction_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            #plot autocovariance function heatmaps along x-y, x-z and y-z planes
            for j in 1:3

                plot_autocovariance_fct_heatmap(autocovariance_fct_direction_plot_dict,
                        save_path*"_"*autocovariance_fct_direction_plot_dict["label"]*"_"*string(j)*"_fixed";
                        save_plot = true,
                        clims = autocovariance_fct_direction_clims,
                        x_y_lims = autocovariance_fct_direction_x_y_lims,
                        sampling_vector_component_to_fix = j,
                        sampling_vector_value_fixed = 0)
            end

        end
        
    end


    #if desired plot spectral density along three different directions
    if plot_spectral_density_along_directions_bool

        #set vector with names to load the files
        naming_vec = string.( [[1,0,0], 
                                [0,1,0],
                                [0,0,1],
                                [1,-1,0],
                                [1,1,1],
                                [1,1,-2]] )

        #set vector with the two collections of directions that will be plotted in two different plots
        range_vec = [1:3, 4:6]

        #set the names of the two plots
        plot_name_vec = ["_along_axes", "_rotated_axes"]

        #loop through data that will be plotted
        for i in eachindex(data_path_vec)

            #create two plots per sample, one along the axes and one along rotated coordinate axes
            for plot_per_sample_nr in 1:2

                #initialize vector for plot dicts
                spectral_density_direction_plot_dict_vec = []

                spectral_density_direction_plot_dict = Dict()

                #loop through the three directions either along the axes or along rotated coordinate axes
                for j in range_vec[plot_per_sample_nr]

                    #load plot dictionary
                    spectral_density_direction_plot_dict = GU.load_h5_dict(
                                                data_path_vec[i]*"_"*naming_vec[j]*"_spectral_density_direction.h5")

                    #if desired adjust label
                    if label_vec  !== nothing
                        spectral_density_direction_plot_dict["label"] = label_vec[i]*" "*naming_vec[j]
                    end

                    #if desired adjust voxel edge length
                    if voxel_edge_length_vec !== nothing
                        spectral_density_direction_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
                    end

                    #add current plot dict to vector
                    push!(spectral_density_direction_plot_dict_vec, spectral_density_direction_plot_dict)

                end

                #set saving path by extracting the current sample name from data path
                save_path_specific = ( first( save_path, findlast('\\', save_path) )
                                *SubString(data_path_vec[i], (findlast('\\', data_path_vec[i]) + 1), length(data_path_vec[i]))
                                *"_direction"
                                *plot_name_vec[plot_per_sample_nr]
                                )

                #plot the data
                plot_spectral_density(spectral_density_direction_plot_dict_vec,
                                        save_path_specific;
                                        save_plot = true)

                #if specified, plot a zoom into the spectral sensity
                if spectral_density_xlims !== nothing
                
                    plot_spectral_density(spectral_density_direction_plot_dict_vec,
                    save_path_specific;
                                save_plot = true,
                                xlims = spectral_density_xlims)
                end

            end

        end

    end


    #if desired plot spectral density heatmaps
    if plot_spectral_density_heatmaps_bool

        #loop through data that will be plotted
        for i in eachindex(data_path_vec)

            #load plot dictionary
            spectral_density_array_dict = GU.load_h5_dict(data_path_vec[i]*"_spectral_density_array.h5")

            #if desired adjust label and voxel edge length
            if label_vec  !== nothing
                spectral_density_array_dict["label"] = label_vec[i]
            end

            #if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                spectral_density_array_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            #plot spectral density heatmaps along x-y, x-z and y-z planes
            for j in 1:3

                plot_spectral_density_heatmap(spectral_density_array_dict,
                        save_path*"_"*spectral_density_array_dict["label"]*"_"*string(j)*"_fixed";
                        save_plot = true,
                        clims = spectral_density_heatmaps_clims,
                        x_y_lims = spectral_density_heatmaps_x_y_lims,
                        wavevector_component_to_fix = j,
                        wavevector_value_fixed = 0)
            end

        end
        
    end

    return

end