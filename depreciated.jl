# there are some depreciated functions

"""
determine uncertainty on the local volume fraction for given window size
"""
function get_uncertainty_local_volume_fract_variance_random_windows_depreciated(nr_dimensions_data::Int64,
                                                    edge_length_window::Int64,
                                                    volume_fract_tot::Float64,
                                                    sampled_voxels_array::Array{Int64},
                                                    local_volume_fract_variance::Float64)
    
   

    # sometimes the uncertainty becomes negative which can only them from the nr of independent samples
    # therefore print this number for each window size
    println("window length: "*string(edge_length_window)*", nr_independent_samples: "*string(nr_independent_samples))
    println("local_volume_fract_variance: "*string(local_volume_fract_variance))

    # calculate  uncertainty on the local volume fraction
    uncertainty_local_volume_fract_variance = sqrt( 1/( window_volume^3 * nr_independent_samples
                                                        * volume_fract_tot * (1 - volume_fract_tot )^2 )
                                                    + 2 * local_volume_fract_variance 
                                                        / (nr_independent_samples - 1) )

    return uncertainty_local_volume_fract_variance

end



"""
plot window volume times local volume fraction uncertainty as a function of window edge length
to determine whether structure is hyperuniform.
If the plot tends to zero for edge length -> infintity, then the structure is hyperuniform
"""
function plot_hyperuniform_criterion_depreciated(nr_dimensions_data::Int64,
                                    window_edge_length_vec::Vector{Int64}, 
                                    local_volume_fract_variance_vec::Vector{Measurements.Measurement};
                                    title="Hyperuniformity test",
                                    save_plot = false,
                                    save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\\",
                                    save_filename="hyperuniform_test.pdf")

    # get vector of window volumes
    window_volume_vec = window_edge_length_vec .^ nr_dimensions_data

    # plot window volume times local volume fraction variance against window edge length
    hyperuniform_plot = Plots.plot(window_edge_length_vec, 
                Measurements.value.(local_volume_fract_variance_vec) .* window_volume_vec, 
                yerr = Measurements.uncertainty.(local_volume_fract_variance_vec) .* window_volume_vec,
                xlabel="window edge length "*Latex.L"L"*" / voxel edge length",
                ylabel=Fmt.format(Latex.L"\sigma_\tau^2(L) \cdot L^{:d} / ", 
                                    nr_dimensions_data)
                                *"voxel edge length"
                                *Fmt.format(Latex.L"^{:d}", nr_dimensions_data) ,
                legend = false, dpi=250, title=title)

    # if specified by the argument, save the plot
    if save_plot
        Plots.savefig(save_path*save_filename)

    else
        Plots.display(hyperuniform_plot)
    end

    return

end


"""
perform the tensor product of two 3d arrays
"""
function tensor_product_3d(array1::Array, array2::Array)

    # get the sizes of both arrays
    size1 = size(array1)
    size2 = size(array2)

    # initialize output array
    array_out = zeros(size1[1]*size2[1], size1[2]*size2[2], size1[3]*size2[3])

    # loop through all dimensions of both arrays and calculate tensor product element-wise
    for i in eachindex(array1[:,1,1])
        for j in eachindex(array1[1,:,1])
            for k in eachindex(array1[1,1,:])
                for l in eachindex(array2[:,1,1])
                    for m in eachindex(array2[1,:,1])
                        for n in eachindex(array2[1,1,:])

                            array_out[size1[1]*i-1+l, size1[2]*j-1+m, size1[3]*k-1+n] = (array1[i,j,k] 
                                                                                       * array1[l,m,n])
                        end
                    end
                end
            end
        end

    end

    return array_out

end



"""
reduce number of array elements along every dimension by a given factor
"""
function reduce_array(array_to_reduce, reduction_factor)

    # consider only those elements on a grid with the reduction factor as a step size
    reduced_array = array_to_reduce[1:reduction_factor:end, 
                                    1:reduction_factor:end, 
                                    1:reduction_factor:end]

    return reduced_array

end




function plot_hyperuniform_criterion(nr_dimensions_data::Int64,
    window_edge_length_vec::Vector{Int64}, 
    local_volume_fract_variance_vec::Vector{Measurements.Measurement};
    title="Hyperuniformity test",
    save_plot = false,
    save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\\",
    save_filename="hyperuniform_test.pdf")

# get vector of window volumes
window_volume_vec = window_edge_length_vec .^ nr_dimensions_data

# plot window volume times local volume fraction variance against window edge length
hyperuniform_plot = Plots.plot(window_edge_length_vec, 
Measurements.value.(local_volume_fract_variance_vec) .* window_volume_vec, 
ribbon = Measurements.uncertainty.(local_volume_fract_variance_vec) .* window_volume_vec,
xlabel="window edge length "*Latex.L"L"*" / voxel edge length",
ylabel=Fmt.format(Latex.L"\sigma_\tau^2(L) \cdot L^{:d} / ", 
    nr_dimensions_data)
*"voxel edge length"
*Fmt.format(Latex.L"^{:d}", nr_dimensions_data) ,
legend = false, dpi=250, title=title)

# if specified by the argument, save the plot
if save_plot
Plots.savefig(save_path*save_filename)

else
Plots.display(hyperuniform_plot)
end

return

end




"""
get a vector with odd measuring window edge lengths
"""
function get_window_edge_length_vec(mean_edge_length_data::Int64; nr_window_sizes::Int64 = 100 )

    # set maximal window edge length
    max_window_edge_length = Int( round( mean_edge_length_data/4 ) )

    # get vector with odd integer window lengths up to mean_edge_length_data/4
    window_edge_length_vec = ( 2 .* Int.( round.( collect(
                                 LinRange(1, max_window_edge_length/2, nr_window_sizes) 
                                 ) ) )  .+ 1 )

    return window_edge_length_vec

end



"""
determine spectral density as a function of the wavenumber k for isotrope 3d binary
media
"""
function get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data::Int64,
                size_data::Tuple, 
                volume_fract_tot::Float64,
                data_binary::Array{Float64};
                nr_sampling_distances::Int64 = get_nr_sampling_distances(mean_edge_length_data),
                nr_measurements_per_distance::Int64 = 10000,
                nr_wavenumbers::Int64 = 200,
                sampling_distance_vec = get_sampling_distance_vec(mean_edge_length_data; 
                                            nr_sampling_distances = nr_sampling_distances),
                autocovariance_fct_vec = get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                            size_data, 
                                            volume_fract_tot,
                                            data_binary;
                                            nr_sampling_distances = nr_sampling_distances,
                                            nr_measurements_per_distance = nr_measurements_per_distance,
                                            sampling_distance_vec = sampling_distance_vec)[2])
                                     

    # get vector of window edge lengths that will be measured
    wavenumber_vec = get_wavenumber_vec(sampling_distance_vec; nr_wavenumbers = nr_wavenumbers)

    # create vector where for each wavenumber the spectral density and its
    # uncertainty will be stored
    spectral_density_vec = Vector{Complex{Measurements.Measurement}}(undef, nr_wavenumbers)

    # for each wavenumber get spectral density and its uncertainty
    for i in eachindex(wavenumber_vec)

        # get vector of local volume fractions and the number of independent samples
        spectral_density_vec[i] = get_spectral_density_isotrope(autocovariance_fct_vec,
                                                                wavenumber_vec[i];
                                                                sampling_distance_vec=sampling_distance_vec)

        # print current calculation status
        println("wavenumber "*string(wavenumber_vec[i])*" done")

    end

    return[wavenumber_vec, spectral_density_vec]
    
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

    # set x label which is the same for all plots
    xlabel = "wavenumber " * Latex.L"k / ( 1/ \mathrm{nm} )"

    # create empty plots for real, imaginary parts and absolute value
    re_plot = Plots.plot(xlabel=xlabel,
                                    ylabel=Latex.L"\mathrm{Re}(\tilde{\chi} (k)) ",
                                    legend = true, dpi=250, title=title)

    im_plot = Plots.plot(xlabel=xlabel,
                                    ylabel=Latex.L"\mathrm{Im}(\tilde{\chi} (k)) ",
                                    legend = true, dpi=250, title=title)

    abs_plot = Plots.plot(xlabel=xlabel,
                                    ylabel=Latex.L"\mathrm{Abs}(\tilde{\chi} (k)) ",
                                    legend = true, dpi=250, title=title)

    # if desired, set xlims
    if xlims !== nothing 

        Plots.plot!(re_plot, xlims=xlims)
        Plots.plot!(im_plot, xlims=xlims)
        Plots.plot!(abs_plot, xlims=xlims)

    end

    # loop through plot dictionaries to add plots
    for plot_dict in plot_dict_vec

        # scale wavenumbers with inverse nanometers
        x_vec = plot_dict["wavenumber_vec"] ./ plot_dict["voxel_edge_length"]

        # plot real part of spectral density
        Plots.plot!(re_plot, x_vec, 
                    real.( Measurements.value.(plot_dict["spectral_density_vec"]) ), 
                    ribbon = real.( Measurements.uncertainty.(plot_dict["spectral_density_vec"]) ),
                    label = plot_dict["label"])

        # plot imaginary part of spectral density
        Plots.plot!(im_plot, x_vec, 
                    imag.( Measurements.value.(plot_dict["spectral_density_vec"]) ), 
                    ribbon = imag.( Measurements.uncertainty.(plot_dict["spectral_density_vec"]) ),
                    label = plot_dict["label"])

        # plot absolute value of spectral density
        abs_spectral_density_vec = abs.( plot_dict["spectral_density_vec"]  )

        Plots.plot!(abs_plot, x_vec, 
                    Measurements.value.(abs_spectral_density_vec), 
                    ribbon = Measurements.uncertainty.(abs_spectral_density_vec),
                    label = plot_dict["label"])


    end

    # if specified by the argument, save the plots
    if save_plot && xlims!== nothing
        Plots.savefig(re_plot, save_path*"_spectral_density_re_zoom.png")
        Plots.savefig(im_plot, save_path*"_spectral_density_im_zoom.png")
        Plots.savefig(abs_plot, save_path*"_spectral_density_abs_zoom.png")

    elseif save_plot
        Plots.savefig(re_plot, save_path*"_spectral_density_re.png")
        Plots.savefig(im_plot, save_path*"_spectral_density_im.png")
        Plots.savefig(abs_plot, save_path*"_spectral_density_abs.png")

    # otherwise display the plots
    else
        Plots.display(re_plot)
        Plots.display(im_plot)
        Plots.display(abs_plot)

    end

    return

end


"""
Get array with vectors for which the autocovariance function will be calculated
"""
function get_sampling_vec_array(size_data::Tuple)

    # determine the maximal sampling distances along the three axes
    max_sampling_distances = Int.( floor.( (size_data .- 1) ./ 2 ))

    # determine size of array where sampling vectors will be stored in
    # along one axis (z direction is chosen here) only positive directions are considered,
    # because negative ones would yield redundant information
    # the fourth dimension contains the 3 vector entries
    sampling_vec_array_size = (2*max_sampling_distances[1] + 1 , 
                                2*max_sampling_distances[2] + 1, 
                                max_sampling_distances[3] + 1,
                                3 )

    # initialize array where sampling vectors will be stored in
    sampling_vec_array = Array{Int64}(undef, sampling_vec_array_size...)

    # fill array with sampling vectors
    for i in -max_sampling_distances[1]:max_sampling_distances[1]
        for j in -max_sampling_distances[2]:max_sampling_distances[2]
            for k in 0:max_sampling_distances[3]

                sampling_vec_array[i + max_sampling_distances[1] + 1,
                                    j + max_sampling_distances[2] + 1,
                                    k + 1,
                                    :] = [i,j,k]

            end
        end
    end

    return sampling_vec_array

end




"""
Get vector of vectors of sampled wavenumbers
"""
function get_sampled_wavenumbers_vec_vec(autocovariance_fct_array_values::Array)

    # get vectors of wavenumbers along all three dimensions
    sampled_wavenumbers_vec_vec = collect( FFTW.rfftfreq.( 
                                            size(autocovariance_fct_array_values) ) )

    # convert to Float64
    sampled_wavenumbers_vec_vec_float = []

    for sampled_wavenumbers_vec in sampled_wavenumbers_vec_vec
        push!(sampled_wavenumbers_vec_vec_float, Float64.( sampled_wavenumbers_vec ))

    end

    return  sampled_wavenumbers_vec_vec_float
    
end


"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density(size_data::Tuple, 
                volume_fract_tot::Float64,
                data_binary::Array{Float64};
                nr_measurements_per_direction::Int64 = 1000,
                sampling_vec_array = get_sampling_vec_array(size_data),
                autocovariance_fct_array = get_autocovariance_fct_by_sampling_vec_array(size_data, 
                            volume_fract_tot,
                            data_binary;
                            sampling_vec_array = sampling_vec_array,
                            nr_measurements_per_direction = nr_measurements_per_direction)[3])
    
    # get values of autocovariance function
    autocovariance_fct_array_values = Measurements.value.( autocovariance_fct_array )

    # determine fourier transform of autocovariance function values
    spectral_density_array = FFTW.rfft(autocovariance_fct_array_values)

    # get tuple of vectors of sampled wavenumbers
    sampled_wavenumbers_vec_vec = get_sampled_wavenumbers_vec_vec(
                            autocovariance_fct_array_values)

    # get array of sampled sampled wavevectors 
    sampled_wavevectors_array = get_vector_array(sampled_wavenumbers_vec_vec)

    return[sampled_wavenumbers_vec_vec, sampled_wavevectors_array, spectral_density_array]
    
end



"""
get vector of sampled wavenumbers
"""
function get_sampled_wavenumbers_vec(direction_vec::Vector{Int64}, autocovariance_fct_vec)

    # determine geometrical length of direction vector
    sampling_distance = sqrt(sum( direction_vec .^2 ))

    # determine sampled wavenumbers vec, where the sampling rate is the inverse of the sampling distance
    sampled_wavenumbers_vec = FFTW.rfftfreq( length( autocovariance_fct_vec ), 1/sampling_distance )

    # convert to float before returning
    return Float64.( sampled_wavenumbers_vec )
end


"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density_along_direction(size_data::Tuple, 
                volume_fract_tot::Float64,
                data_binary::Array{Float64};
                direction_vec::Vector{Int64} = [0,0,1],
                nr_measurements_per_direction::Int64 = 1000,
                sampling_vec_array = get_sampling_vec_array(size_data),
                autocovariance_fct_array = get_autocovariance_fct_by_sampling_vec_array(size_data, 
                            volume_fract_tot,
                            data_binary;
                            sampling_vec_array = sampling_vec_array,
                            nr_measurements_per_direction = nr_measurements_per_direction)[3])
    
    # extract autocovariance function vector along given direction
    autocovariance_fct_vec = get_autocovariance_fct_along_direction_vec(direction_vec, autocovariance_fct_array)

    # determine fourier transform of autocovariance function values along direction
    spectral_density_value_vec = FFTW.rfft( Measurements.value.(autocovariance_fct_vec) )

    # determine uncertainty of spectral_density. Since Fourier transform is linear,
    # the output's uncertainty is the Fourier transform of the input's uncertainty
    spectral_density_uncertainty_vec  = FFTW.rfft( Measurements.uncertainty.(autocovariance_fct_vec) )

    # combine values and uncertainties into measurement type
    spectral_density_vec = complex.(Measurements.measurement.( 
                                        real.(spectral_density_value_vec), 
                                        real.(spectral_density_uncertainty_vec)), 
                                    Measurements.measurement.(
                                        imag.(spectral_density_value_vec), 
                                        imag.(spectral_density_uncertainty_vec)))
 
    # get vector of sampled wavenumbers
    sampled_wavenumbers_vec = get_sampled_wavenumbers_vec(direction_vec,
                                            autocovariance_fct_vec)

    return sampled_wavenumbers_vec, spectral_density_vec
    
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
    clims = nothing,
    wavevector_component_to_fix::Int64 = 3,
    wavevector_value_fixed = 0)

    # discriminate between different wavevector components that are fixed
    if wavevector_component_to_fix == 1

        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["sampled_wavenumbers_vec_vec"][2]
        wavenumber_vec_y = plot_dict["sampled_wavenumbers_vec_vec"][3]

        # find index of fixed wavenumber value
        wavevector_fixed_index = argmin( abs.( 
                                        plot_dict["sampled_wavenumbers_vec_vec"][1]
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][wavevector_fixed_index,:,:] 
        
        # set labels and title for the plot
        xlabel = Latex.L"k_y / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_z / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_x = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["sampled_wavenumbers_vec_vec"][1][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")


    elseif wavevector_component_to_fix == 2
        
        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["sampled_wavenumbers_vec_vec"][1]
        wavenumber_vec_y = plot_dict["sampled_wavenumbers_vec_vec"][3]

        # find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
                                    plot_dict["sampled_wavenumbers_vec_vec"][2]
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][:,wavevector_fixed_index,:] 

        
        # set labels and title for the plot
        xlabel = Latex.L"k_x / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_z / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_y = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["sampled_wavenumbers_vec_vec"][2][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")

    elseif wavevector_component_to_fix == 3
        
        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["sampled_wavenumbers_vec_vec"][1]
        wavenumber_vec_y = plot_dict["sampled_wavenumbers_vec_vec"][2]

        # find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
                                    plot_dict["sampled_wavenumbers_vec_vec"][3] 
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][:,:,wavevector_fixed_index] 

        # set labels and title for the plot
        xlabel = Latex.L"k_x / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_y / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_z = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["sampled_wavenumbers_vec_vec"][3][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")

    else
        @error ("Wavevector component to fix must 
                be 1, 2 or 3, but is "*string(wavevector_component_to_fix))
    end

    # scale x and y axes in units of 1/nm
    x_axis = FFTW.fftshift( wavenumber_vec_x ) / plot_dict["voxel_edge_length"]
    y_axis = FFTW.fftshift(wavenumber_vec_y ) / plot_dict["voxel_edge_length"]

    # permute dimensions of spectral density array, such that they match the axes
    spectral_density_2d_permuted_array = FFTW.fftshift(permutedims(spectral_density_2d_array) )

    # create plots
    abs_plot = Plots.heatmap(x_axis,
                                y_axis,
                                abs.(spectral_density_2d_permuted_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = Latex.L"\mathrm{Abs}( \tilde{\chi} (\vec{k}) ) " ,
                                legend = true, dpi=250, title=title,
                                c = :bluesreds,
                                aspect_ratio = :equal)

    re_plot = Plots.heatmap(x_axis,
                                y_axis,
                                real.(spectral_density_2d_permuted_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = Latex.L"\mathrm{Re}( \tilde{\chi} (\vec{k}) ) " ,
                                legend = true, dpi=250, title=title,
                                c = :bluesreds,
                                aspect_ratio = :equal)

    im_plot = Plots.heatmap(x_axis,
                                y_axis,
                                imag.(spectral_density_2d_permuted_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = Latex.L"\mathrm{Im}( \tilde{\chi} (\vec{k}) ) " ,
                                legend = true, dpi=250, title=title,
                                c = :bluesreds,
                                aspect_ratio = :equal)

    # set clims if desired
    if clims !== nothing
        Plots.heatmap!(abs_plot, clims = clims)
        Plots.heatmap!(re_plot, clims = clims)
        Plots.heatmap!(im_plot, clims = clims)
    end

    # if specified by the argument, save the plot
    if  save_plot
        Plots.savefig(abs_plot, save_path*"_spectral_density_abs.png")
        Plots.savefig(re_plot, save_path*"_spectral_density_re.png")
        Plots.savefig(im_plot, save_path*"_spectral_density_im.png")

    # otherwise display the plot
    else
        Plots.display(abs_plot)

    end

    return
end



"""
get wavenumber vector such that the discrete Fourier transform can be determined
properly based on the sampled distances
"""
function get_wavenumber_vec(sampling_distance_vec::Vector;
                                nr_wavenumbers::Int64 = 200)

    # determine maximal wavenumber based on the Nyquist–Shannon sampling theorem
    # (Nyquist frequency, https://en.wikipedia.org/wiki/Nyquist%E2%80%93Shannon_sampling_theorem)
    wavenumber_max = 1 / ( 2 * Statistics.mean( sampling_distance_vec[2:end] 
                                                    .- sampling_distance_vec[1:end-1] ) )

    # get minimal wavenumber from maximal sampled distance
    wavenumber_min = 2 * pi / maximum(sampling_distance_vec)

    # get vector of float wavenumbers
    wavenumber_vec = collect(LinRange(wavenumber_min, wavenumber_max, nr_wavenumbers))

    return wavenumber_vec
    
end


"""
get vector of sampled wavenumbers
"""
function get_sampled_wavenumbers_vec(direction_vec::Vector, 
                                    autocovariance_fct_along_direction_vec::Vector)

    # determine geometrical length of direction vector
    sampling_distance = sqrt(sum( direction_vec .^2 ))

    # determine nr of sampling distances
    nr_sampling_distances = length(autocovariance_fct_along_direction_vec)

    # determine fundamental wavenumber
    fundamental_wavenumber = 2*pi/(sampling_distance*nr_sampling_distances)

    # determine sampled wavenumbers vec
    sampled_wavenumbers_vec = collect(1:floor(nr_sampling_distances/2)) .* fundamental_wavenumber

    return sampled_wavenumbers_vec
end


"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density_array_by_fft(size_data::Tuple, 
                volume_fract_tot::Float64,
                data_binary::Array{Float64};
                nr_measurements_per_direction::Int64 = 1000,
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(size_data),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                autocovariance_fct_array = get_autocovariance_fct_by_sampling_vec_array(size_data, 
                            volume_fract_tot,
                            data_binary;
                            sampling_vec_array = sampling_vec_array,
                            nr_measurements_per_direction = nr_measurements_per_direction)[3],
                save_result = false,
                save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\sample_name",
                voxel_edge_length = 10,
                label = "some structure")

    # save autocovariance_fct array and sampling distances to dictionary
    autocovariance_fct_direction_dict = Dict("sampling_vec_array" => sampling_vec_array,
                                            "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                                            "autocovariance_fct_array" => autocovariance_fct_array,
                                            "voxel_edge_length" => voxel_edge_length,
                                            "label" => label)

    # autocovariance fct was only calculated for half space with z>0 which, due to its mirror
    # symmetry, is sufficient. For the FFT it is however necessary (I think?) to get the values
    # for the entire space, which is calculated here
    (complete_sampling_distance_vec_vec, 
        complete_sampling_vec_array, 
        complete_autocovariance_fct_array) = mirror_autoautocovariance_fct_by_sampling_vec_array(
                                                autocovariance_fct_direction_dict)

    # get values of autocovariance function
    complete_autocovariance_fct_array_values = Measurements.value.( complete_autocovariance_fct_array )

    # determine fourier transform of autocovariance function values
    spectral_density_array = FFTW.rfft(complete_autocovariance_fct_array_values)

    # get tuple of vectors of sampled wavenumbers
    wavenumber_vec_vec = get_wavenumber_vec_vec(complete_autocovariance_fct_array_values)

    # get array of sampled sampled wavevectors 
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    # save results if desired
    if save_result
        
        # create dict to save
        plot_dict = Dict("wavevector_array" => wavevector_array,
                            "wavenumber_vec_vec" => wavenumber_vec_vec,
                            "spectral_density_array" => spectral_density_array,
                            "voxel_edge_length" => voxel_edge_length,
                            "label" => label)

        save_dict_to_h5(copy(plot_dict);
                        save_path=save_path*"_spectral_density_array.h5")

    end

    return[wavenumber_vec_vec, wavevector_array, spectral_density_array]
    
end



"""
load voxel size corrected data from h5 file
"""
function load_binary_data(data_path_h5::String)

    # load data dictionary (this is how h5 files are structured)
    data_binary_dict = FileIO.load(data_path_h5)

    # get data from dictionary. It is expected under the "data" key
    data_binary = data_binary_dict["data_binary"]

    return data_binary

end





"""
Update bond stretching energy for a given bond
"""
function update_local_bond_stretching_energy_keating(graph_dict::Dict, bond::Tuple{Int64, Int64})

    # get bond stretching energy
    bond_stretching_energy = (3/16 * ( 
            graph_dict["spatial_network"][bond...]["distance_squared"] - 1 
                                            )^2 ) 

    # save to dict
    graph_dict["spatial_network"][bond...]["bond_stretching_energy"] = bond_stretching_energy
    
    return graph_dict

end


"""
Update bond bending energy for a given vertex
"""
function update_local_bond_bending_energy_keating!(graph_dict::Dict, vertex_label::Int64)

    # get vector of neighbor labels
    neighbor_label_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex_label))

    # initialize bond bending sum
    bond_bending_sum = 0

    # loop through all bond combinations
    for j in 1:graph_dict["coordination_nr"]

        for k in j+1:graph_dict["coordination_nr"]

            bond_bending_sum += ( 3/8 * graph_dict["bond_bending_const"] 
                * ( LinearAlgebra.dot( sign(neighbor_label_vec[j] - vertex_label) .* 
                            graph_dict["spatial_network"][vertex_label, neighbor_label_vec[j]]["vector"], 
                            sign(neighbor_label_vec[k] - vertex_label) .* 
                            graph_dict["spatial_network"][vertex_label, neighbor_label_vec[k]]["vector"]
                                     ) + 1/3 )^2 )
            
        end

    end

    # calculate bond bending energy
    bond_bending_energy = 3/8 * graph_dict["bond_bending_const"] * bond_bending_sum
    
    # save to dict
    graph_dict["spatial_network"][vertex_label]["bond_bending_energy"] = bond_bending_energy

    return graph_dict

end


"""
Calculate the total energy of a spatial network
"""
function total_energy_keating(graph_dict::Dict)

    total_energy = 0

    # loop through all vertices and sum bond bending energies
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        total_energy += graph_dict["spatial_network"][vertex]["bond_bending_energy"]
    end

    # loop through all bonds and sum bond stretching energies
    for edge in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        total_energy += graph_dict["spatial_network"][edge...]["bond_stretching_energy"]
    end

    return total_energy
end



"""
Calculate the local Keating energy for a given vertex from 
the bond bending and stretching energies stored in the dictionary by fully
considering its bonds and not sharing their energy between the two vertices
"""
function local_energy_keating(vertex_label::Int64, graph_dict::Dict)

    # initialize local energy
    local_energy = graph_dict["spatial_network"][vertex_label]["bond_bending_energy"]

    # sum bond stretching energy contributions by considering that each bond
    # is shared by two vertices
    for neighbor in MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex_label)

        local_energy += graph_dict["spatial_network"][vertex_label, neighbor]["bond_stretching_energy"]

    end
    
    return local_energy

end



"""
Calculate the local Keating energy for a given vertex from 
the bond bending and stretching energies stored in the dictionary by
considering the bonds as shared and using half their bond energies
"""
function local_energy_keating_shared_bonds(vertex_label::Int64, graph_dict::Dict)

    # initialize local energy
    local_energy = graph_dict["spatial_network"][vertex]["bond_bending_energy"]

    # sum bond stretching energy contributions by considering that each bond
    # is shared by two vertices
    for neighbor in MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex_label)

        local_energy += 1/2 * graph_dict["spatial_network"][vertex_label, neighbor]["bond_stretching_energy"]

    end
    
    return local_energy
end

    
"""
add bond bending and stretching energies to all vertices and bonds
"""
function add_energies_to_spatial_network!(graph_dict;
    update_local_vertex_energy_fct! = update_local_bond_bending_energy_keating!,
    update_local_bond_energy_fct! = update_local_bond_stretching_energy_keating!,
    total_energy_fct = total_energy_keating)

    # add bond bending energies to vertices
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        graph_dict = update_local_vertex_energy_fct!(graph_dict, vertex)

    end

    # add bond stretching energies to bonds
    for edge in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        graph_dict = update_local_bond_energy_fct!(graph_dict, edge)

    end

    # add total energy to graph dict
    graph_dict["total_energy"] = total_energy_fct(graph_dict)

    return graph_dict

end




"""
calculate the total bond bending and stretching Keating
energy of a given network graph
"""
function total_bending_and_stretching_energy_keating(graph_dict::Dict)

    # first calculate bond stretching energy by looping through edges
    bond_stretching_sum = 0

    for edge in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        # get current edge length
        edge_length = graph_dict["spatial_network"][edge...]["distance"]

        # add bond stretching energy of current bond to sum
        bond_stretching_sum += (edge_length^2 - 1 )^2

    end

    # multiply sum with prefactors to determine bond stretching energy
    bond_stretching_energy = 3/16 * bond_stretching_sum

    # calculate bond-bending energy

    
    # get iterator of bond combinations
    bond_combinations_iter = Combinatorics.combinations(
                        collect(1:graph_dict["coordination_nr"]), 2)

    bond_bending_sum = 0

    for current_vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        # get list of neighbors
        neighbors_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], current_vertex))

        # get vectors pointing to neighbors_vec
        bond_vectors_mat = Matrix{Float64}(undef, graph_dict["nr_dimensions"], graph_dict["coordination_nr"])

        for i in eachindex(neighbors_vec)

            # get vector pointing from the current vertex to neighbor
            # where the direction needs to be flipped if the vertex code of
            # the current vertex is inferior to the one of its neighbor
            bond_vectors_mat[:,i] = (sign(neighbors_vec[i] - current_vertex)
                        .* graph_dict["spatial_network"][current_vertex, neighbors_vec[i]]["vector"])

        end

        # loop through all bond combinations to determine corresponding bond-bending energy
        for bond_combination in bond_combinations_iter

            bond_bending_sum += ( LinearAlgebra.dot(bond_vectors_mat[:,bond_combination[1]], 
                                                bond_vectors_mat[:,bond_combination[2]]) 
                                + 1/3 )^2

        end

    end

    # multiply sum with prefactors to determine bond bending energy
    bond_bending_energy = 3/8 * graph_dict["bond_bending_const"]  * bond_bending_sum

    # get total energy
    total_energy = bond_stretching_energy + bond_bending_energy

    return [total_energy, bond_stretching_energy, bond_bending_energy]
end


"""
This function performs a bond switch on a graph.
The argument switched_bond is a tuple of two integers
which is the edge type of the MetaGraphsNext package
"""
function switch_bond_wrong!(graph_dict::Dict, switched_bond::Tuple{Int64, Int64} )

    # break the original bond
    MetaGraphsNext.rem_edge!(graph_dict["spatial_network"], switched_bond...)

    # find the other vertex's neighbors that are the closest to the current vertex
    closest_other_vertices_neighbor_vec = Vector{Int64}(undef, 2)
    vector_to_closest_other_vertices_neighbor_vec = Vector{Vector{Float64}}(undef, 2)
    distance_to_closest_other_vertices_neighbor_vec = Vector{Float64}(undef, 2)

    for i in 1:2

        # get the other vertex's neighbors
        other_vertex_neighbors_vec = MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], switched_bond[(3-i)])

        # get the vertexic position of bond vertex and store it three times
        vertexic_position_vec = graph_dict["spatial_network"][switched_bond[i]]["position"]

        # determine the closest one of the other vertex's neighbors
        closest_distance = Inf

        for neighbor in other_vertex_neighbors_vec

            # determine vector to current neighbor
            if switched_bond[i] < neighbor
                vector_to_neighbor =  get_distance_vector_pbc(
                        vertexic_position_vec,
                        graph_dict["spatial_network"][neighbor]["position"],
                        graph_dict["supercell_edge_length"] )
            else
                vector_to_neighbor =  get_distance_vector_pbc(
                        graph_dict["spatial_network"][neighbor]["position"],
                        vertexic_position_vec,
                        graph_dict["supercell_edge_length"] )
            end

            # determine length of vector to neighbor
            distance_to_neighbor = LinearAlgebra.norm(vector_to_neighbor)

            # store current neighbor if its closer than the previous neighbors
            if distance_to_neighbor < closest_distance
                closest_other_vertices_neighbor_vec[i] = neighbor
                vector_to_closest_other_vertices_neighbor_vec[i] = vector_to_neighbor
                distance_to_closest_other_vertices_neighbor_vec[i] = distance_to_neighbor
            end

        end

    end

    # create the two new bonds
    for i in 1:2

        graph_dict["spatial_network"][switched_bond[i], closest_other_vertices_neighbor_vec[i]] = Dict(
            "vector" => vector_to_closest_other_vertices_neighbor_vec[i], 
            "distance_squared" => distance_to_closest_other_vertices_neighbor_vec[i]^2 )
    end

    return graph_dict

end



"""
Get four vertices in one line around a central bond,
that is one vertex on each side of the bond
"""
function get_four_vertices_around_bond(graph_dict::Dict, 
                                switched_bond::Tuple{Int64, Int64})

    # store four vertices that sit in one line with the switched bond in the center
    vertices_in_line = zeros(Int64, 4)
    vertices_in_line[2:3] = collect(switched_bond)

    # loop through bond vertices
    for i in 1:2

        # pick a random neighbor to which the bond will be cut        
        vertices_in_line[-2+3*i] = setdiff(collect(MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], switched_bond[i]
                            )) , vertices_in_line)[rand(1:graph_dict["coordination_nr"])]

    end

    return vertices_in_line
end



"""
Approximately relax a single vertex by efficiently determining
the approximated coordinate shift. The corresponding methdo is explained in
10.1142/S0217984987000065
"""
function relax_single_vertex_keating_efficiently!(graph_dict::Dict,
    vertex_to_relax::Int64;
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    update_total_energy::Bool = false)

    # get energy gradient at current vertexic position
    gradient = gradient_keating_efficient(graph_dict, vertex_to_relax)

    # get energy hessian at current vertexic position
    hessian = hessian_keating_efficient(graph_dict, vertex_to_relax)

    # calculate translation vector to approximate energy minimum
    translation_vector = .- LinearAlgebra.inv(hessian)*gradient

    # move vertex 
    graph_dict = move_vertex!(graph_dict, 
                            vertex_to_relax, 
                            translation_vector;
                            update_total_energy = update_total_energy)


    return graph_dict

end


"""
Get mesh from network
"""
function get_mesh_from_network(graph_dict::Dict; bond_radius::Real = 0.05)

    cylinders_merged = nothing

    # loop through bonds
    for bond in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        # get bond's start and target positions and its direction vector
        start_pos = graph_dict["spatial_network"][bond[1]]["position"]
        target_pos = graph_dict["spatial_network"][bond[2]]["position"]
        # direction_vec = graph_dict["spatial_network"][bond...]["vector"]

        # create cylinder surface object
        cylinder_surface = Meshes.CylinderSurface(start_pos, target_pos, bond_radius)
        

        if cylinders_merged == 0
            cylinders_merged = cylinder_surface
        else
            # merge current cylinder surface with previous ones
            cylinders_merged = Meshes.merge(cylinders_merged, cylinder_surface)
        end

    end

    # discretize object of merged cylinders
    cylinder_mesh = Meshes.discretize(cylinders_merged)

    return cylinder_mesh
end


"""
Fully relax a cluster of vertices
"""
function relax_cluster_keating!(graph_dict::Dict,
    cluster_dict::Dict; 
    nr_max_relaxation_cycles::Int64 = 25,
    break_at_relative_cluster_energy_change::Float64 = 0.001,
    reject_during_relaxation_cycle_threshold::Int64 = 5,
    initial_cluster_energy  = cluster_dict["cluster_energy"],
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    update_total_energy::Bool = false,
    track_cluster_energy::Bool = false)

    # track cluster energy during relaxatio if desired
    if track_cluster_energy
        cluster_energy_vec = [cluster_dict["cluster_energy"]]
    end

    # perform the given number of relaxation cycles
    for cycle_nr in 1:nr_max_relaxation_cycles

        # store previous cluster energy
        previous_cluster_energy = cluster_dict["cluster_energy"]

        graph_dict, cluster_dict = relax_cluster_one_cycle_keating!(graph_dict, 
        cluster_dict;
        relax_efficiently = relax_efficiently,
        relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
        relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
        update_cluster_energy = true )

        # save cluster energy if desired
        if track_cluster_energy
            push!(cluster_energy_vec, cluster_dict["cluster_energy"]) 
        end

        # break if cluster energy changes less than the given threshold
        relative_cluster_energy_change = (
            abs((previous_cluster_energy - cluster_dict["cluster_energy"])
                    /cluster_dict["cluster_energy"]))

        if (relative_cluster_energy_change < break_at_relative_cluster_energy_change 
                && cycle_nr > reject_during_relaxation_cycle_threshold)
            println("Breaking at cycle nr "*string(cycle_nr))
            break
        end

        # if cycle nr is above the given threshold, check if the relaxation can 
        # be rejected before full relaxation by estimating the final energy
        if cycle_nr > reject_during_relaxation_cycle_threshold

            # to be implemented

        end

    end

    # update total energy if desired
    if update_total_energy
        graph_dict["total_energy"] = (graph_dict["total_energy"] 
                                    + cluster_dict["cluster_energy"]
                                    - initial_cluster_energy)

        graph_dict["total_energy_up_to_date"] = true
    else
        graph_dict["total_energy_up_to_date"] = false
    end
    
    if track_cluster_energy
        return [graph_dict, cluster_energy_vec]
    else
        return graph_dict
    end
end


"""
Pick a random bond that has not been declined since the last accepted move
"""
function get_random_bond(graph_dict::Dict; declined_bonds = [], seed = Nothing)

    # set seed if desired
    if seed !== Nothing
        Random.seed!(seed)
    end

    # determine nr of bonds
    nr_bonds = graph_dict["nr_vertices"] * graph_dict["coordination_nr"] / 2 

    # check if all bonds have been attempted already
    if length(declined_bonds) == nr_bonds
            @warn "All bonds have been attempted without success"
        random_bond = []

    # if the list of declined bonds is already very long
    # pick one of the remaining ones
    elseif length(declined_bonds) > nr_bonds/2
        all_bonds_vec = collect(
                MetaGraphsNext.edge_labels(graph_dict["spatial_network"]))

        random_bond = rand(all_bonds_vec)

    # otherwise get random bond without listing all bonds
    else

        # pick a random vertex
        vertex_1 = rand(1:graph_dict["nr_vertices"])

        # pick a random neighbor
        vertex_2 = collect(MetaGraphsNext.neighbor_labels(
                            graph_dict["spatial_network"], vertex_1)
                            )[rand(1:graph_dict["coordination_nr"])]

        # create bond
        random_bond = Tuple(sort([vertex_1, vertex_2]))

        # find new bond if current one was already declined
        if random_bond in declined_bonds
            random_bond = get_random_bond(graph_dict; declined_bonds = declined_bonds)
        end
    end

    return random_bond
end


"""
This function performs a bond switch on a graph.
The argument switched_bond is a tuple of two integers
which is the edge type of the MetaGraphsNext package
"""
function switch_bond!(graph_dict::Dict,
    switched_bond::Tuple{Int64, Int64} )

    # find the other vertex's neighbors that are the closest to the current vertex
    new_bond_vertex_vec = Vector{Int64}(undef, 2)
    vector_to_new_bond_vertex_vec = Vector{Vector{Float64}}(undef, 2)
    distance_to_new_bond_vertex_vec = Vector{Float64}(undef, 2)

    # get vectors of original neighbors
    original_neighbors_vec_vec = [collect(MetaGraphsNext.neighbor_labels(
        graph_dict["spatial_network"], switched_bond[1]) ),
        collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[2]) )]

    for i in 1:2

        # get the vertex position of bond vertex
        vertex_position_vec = graph_dict["spatial_network"][switched_bond[i]]["position"]

        # get the other bond vertex's neighbors excluding 
        # the bond vertex and the bond vertex's neighbors
        considered_new_bond_vertices_vec = setdiff(original_neighbors_vec_vec[3-i], 
                                            switched_bond[i], 
                                            original_neighbors_vec_vec[i])

        # break if there are no possible new bond vertices
        if considered_new_bond_vertices_vec == []
            new_bond_vec = []
            return [graph_dict, new_bond_vec]
        
        # otherwise, pick a random new bond vertex
        else
            new_bond_vertex_vec[i] = rand(considered_new_bond_vertices_vec)
        end

        # determine vector to new bond vertex
        if switched_bond[i] < new_bond_vertex_vec[i]
            vector_to_new_bond_vertex_vec[i] = get_distance_vector_pbc(
                    vertex_position_vec,
                    graph_dict["spatial_network"][new_bond_vertex_vec[i]]["position"],
                    graph_dict["supercell_edge_length"] )
        else
            vector_to_new_bond_vertex_vec[i] = get_distance_vector_pbc(
                    graph_dict["spatial_network"][new_bond_vertex_vec[i]]["position"],
                    vertex_position_vec,
                    graph_dict["supercell_edge_length"] )
        end

        # determine length of vector to new bond vertex
        distance_to_new_bond_vertex_vec[i] = LinearAlgebra.norm(vector_to_new_bond_vertex_vec[i])

    end

    # create vector to save new bonds
    new_bond_vec = Vector{Tuple{Int64, Int64}}(undef, 2)

    # for each bond vertex, break bond to one neighbor and reconnect to
    # random neighbor of the other vertex
    for i in 1:2

        MetaGraphsNext.rem_edge!(graph_dict["spatial_network"],
            switched_bond[i], new_bond_vertex_vec[3-i])

        new_bond_vec[i] = (switched_bond[i], new_bond_vertex_vec[i])

        graph_dict["spatial_network"][new_bond_vec[i]...] = Dict(
            "vector" => vector_to_new_bond_vertex_vec[i], 
            "distance_squared" => distance_to_new_bond_vertex_vec[i]^2 )
    end

    # note, that total energy is not up to date any more
    graph_dict["total_energy_up_to_date"] = false

    return [graph_dict, new_bond_vec]

end



"""
Get effective hyperuniformity parameter which is the structure factor
at zero momentum normalized by the height of the first peak in the structure factor
as defined in equation 251 in 10.1016/j.physrep.2018.03.001
"""
function get_effective_hyperuniformity_parameter(structure_factor_dict::Dict)

    # locate first peak of structure factor
    pks, vals = Peaks.findmaxima(structure_factor_dict["structure_factor_vec"])

    # cut structure factor data at momentum just above first peak
    structure_factor_cut_vec = structure_factor_dict["structure_factor_vec"][1:pks[1]+1]
    wavenumber_cut_vec = structure_factor_dict["wavenumber_vec"][1:pks[1]+1]

    # set the order of the fitted polynomial
    polynomial_order = 3

    # fit polynomial of given order to cut data
    polynomial_fit = Polynomials.fit(wavenumber_cut_vec, 
                                    structure_factor_cut_vec,
                                    polynomial_order)

    # get extrapolated structure factor at zero momentum
    structure_factor_zero_momentum = polynomial_fit(0)

    # get the two critical momenta where the fitted structure factor is extremal
    critical_momenta = (
    ((-polynomial_fit[2]) 
        .+ [-1, +1 ] .* (sqrt(polynomial_fit[2]^2-3*polynomial_fit[3]*polynomial_fit[1])) )
    ./ (3*polynomial_fit[3]) )

    # get fitted structure factor at (first) peak
    structure_factor_first_peak = maximum( polynomial_fit.(critical_momenta) )

    # get hyperuniformity parameter
    hyperuniformity_parameter = structure_factor_zero_momentum/structure_factor_first_peak

    return [hyperuniformity_parameter, polynomial_fit]
end



"""
Get vector of mean values of q_l (rotationally invariant Steinhardt local 
bond order parameters) for the entire network and for all parameters l up 
to l_max where l is the index of the spherical harmonic Y_{lm}.
"""
function get_q_l_total_network_mean_dict(graph_dict::Dict,
    l_max::Int64)

    # initialize dictionary of q_l averaged over entire network with all values
    # set to 0
    q_l_total_network_mean_dict = Dict{Int64, Float64}()

    for l in 0:l_max
        q_l_total_network_mean_dict[l] = 0.0

    end

    # loop through vertices
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        # get vector of steinhardt order parameters for current vertex
        q_l_averaged_single_vertex_dict = (
            get_q_l_averaged_single_vertex_dict(
                graph_dict,
                vertex,
                l_max))

        # for each l, add current vertex' contribution to sum of all vertices
        for l in 0:l_max
            q_l_total_network_mean_dict[l] += (1/graph_dict["nr_vertices"] 
                                        * q_l_averaged_single_vertex_dict[l])
    
        end

    end

    return q_l_total_network_mean_dict
end



"""
Save spatial network to a DOT format file 
"""
function save_spatial_network_to_dot(spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    save_path::String 
        = raw"..\structures\random_networks\\")

    # open new file
    open(save_path*filename*".gv", "w") do opened_file

        # write header
        write(opened_file, "graph T {\n")

        # loop through vertices
        for vertex in MetaGraphsNext.labels(spatial_network)

            # write vertex
            write(opened_file, Format.format("    {1} [position = [{2}, {3}, {4}]];\n",
                vertex,
                spatial_network[vertex]["position"][1],
                spatial_network[vertex]["position"][2],
                spatial_network[vertex]["position"][3]))

        end

        # loop through edges
        for edge in MetaGraphsNext.edge_labels(spatial_network)

            # write edge
            write(opened_file, Format.format("    {1} -- {2} [vector = [{3}, {4}, {5}], distance_squared = {6}];\n", 
            edge[1], edge[2], 
            spatial_network[edge...]["vector"][1],
            spatial_network[edge...]["vector"][2],
            spatial_network[edge...]["vector"][3],
            spatial_network[edge...]["distance_squared"]))

        end

        # write footer
        write(opened_file, "}\n")

    end

    return
end 


"""
Save spatial network to a DOT format file and the rest of graph_dict and
evolution_dict to an h5 file
"""
function save_graph_to_h5_and_dot(graph_dict::Dict,
    filename::String;
    evolution_dict = nothing,
    save_path::String 
        = raw"..\structures\random_networks\\")

    # save evolution dict if passed
    if evolution_dict !== nothing
        GU.save_dict_to_h5(evolution_dict;
            save_path=save_path*filename*"_evolution.h5")
    end

    # create copy of graph_dict to not change the original file
    graph_dict_to_save = deepcopy(graph_dict)

    # save graph to dot format
    save_spatial_network_to_dot(graph_dict["spatial_network"], filename; save_path=save_path)

    # remove spatial_network from graph_dict
    delete!(graph_dict_to_save, "spatial_network")

    # save graph dict
    GU.save_dict_to_h5(graph_dict_to_save; save_path=save_path*filename*".h5")

    return
end


"""
Load spatial network from a DOT format file 
"""
function load_spatial_network_from_dot(dict_path::String)

    # create an empty network graph where vertexic positions and edge vectors will be stored
    spatial_network = MetaGraphsNext.MetaGraph(Graphs.Graph(); 
                                        label_type = Int64,
                                        vertex_data_type = Dict{String, Any},
                                        edge_data_type = Dict{String, Any} )

    # read file contents, one line at a time 
    open(dict_path) do opened_file

        # read until end of file
        while ! eof(opened_file) 
        
            # read a new / next line for every iteration		 
            line_string = readline(opened_file)	

            # save edge to graph
            if occursin(" -- ", line_string)

                # get start vertex
                start_vertex_string, rest_string = split(line_string, " -- ")
                start_vertex = parse(Int64, start_vertex_string)

                # get end vertex
                end_vertex_string, rest_string = split(rest_string, " [vector = [")
                end_vertex = parse(Int64, end_vertex_string)

                # get vector and distance squared
                vector_string, rest_string = split(rest_string, "], distance_squared =")
                vector = parse.(Float64, split(vector_string, ", "))

                distance_squared = parse(Float64, rest_string[1:end-3])

                # add edge to graph
                spatial_network[start_vertex, end_vertex] = Dict("vector" => vector, "distance_squared" => distance_squared)
            
            # save vertex to graph
            elseif occursin("position", line_string)

                # get vertex and position
                vertex_string, rest_string = split(line_string, "[position = [")
                vertex = parse(Int64, vertex_string)
                position = parse.(Float64, split( rest_string[1:end-3], ", "))

                # add vertex to graph
                spatial_network[vertex] = Dict("position" => position)
            
            end
        end
    end
    
    return spatial_network
end


"""
Load graph and its properties from a DOT file and a h5 dictionary
"""
function load_graph_from_h5_and_dot(dict_path_without_format::String)

    # load spatial network in MGformat
    spatial_network = load_spatial_network_from_dot(
            dict_path_without_format*".gv")

    # load rest of graph dict
    graph_dict = GU.load_h5_dict(dict_path_without_format*".h5")

    # add spatial network key to graph dict
    graph_dict["spatial_network"] = spatial_network

    return graph_dict
end


"""
Convert a graph in MGformat to a dot format
"""
function convert_MGformat_to_dot(
    filename::String;
    save_path::String 
        = raw"..\structures\random_networks\\")

    # load graph dict
    graph_dict = load_graph_from_h5_and_MGformat(save_path*filename)

    # save graph to dot format
    save_spatial_network_to_dot(graph_dict["spatial_network"], filename; save_path=save_path)

    return
end


"""
Convert all files in a directory in MGformat to dot format
"""
function convert_all_files_in_directory_MGformat_to_dot(directory_path::String)

    # get all files in directory
    filenames = readdir(directory_path)

    # loop through files
    for filename in filenames

        if endswith(filename, ".mg")

            # convert file to dot format
            convert_MGformat_to_dot(filename[1:end-3]; save_path=directory_path)
        end

    end

    return
    
end


"""
Get extremum of a quadratic function given by three points
"""
function get_quadratic_fct_extremum(x_vec::Vector, y_vec::Vector)

    # get coefficients of quadratic function
    a = (x_vec[3] * (-y_vec[1] + y_vec[2]) 
    + x_vec[2] * (y_vec[1] - y_vec[3]) 
    +  x_vec[1] * (-y_vec[2] + y_vec[3])
    )/((x_vec[1] - x_vec[2]) * (x_vec[1] - x_vec[3]) * (x_vec[2] - x_vec[3]))

    b = (x_vec[3]^2 * (y_vec[1] - y_vec[2]) 
    + x_vec[1]^2 * (y_vec[2] - y_vec[3]) 
    + x_vec[2]^2 * (-y_vec[1] + y_vec[3])
    )/((x_vec[1] - x_vec[2]) * (x_vec[1] - x_vec[3]) * (x_vec[2] - x_vec[3]))

    c = (x_vec[1] * x_vec[3] * (-x_vec[1] + x_vec[3]) * y_vec[2] 
    + x_vec[2]^2 * (x_vec[3] * y_vec[1] - x_vec[1] * y_vec[3]) 
    + x_vec[2] * (-x_vec[3]^2 * y_vec[1] + x_vec[1]^2 * y_vec[3])
    )/((x_vec[1] - x_vec[2]) * (x_vec[1] - x_vec[3]) * (x_vec[2] - x_vec[3]))

    # get extremum of quadratic function
    x_extremum = -b/(2*a)
    y_extremum = a*x_extremum^2 + b*x_extremum + c

    return [x_extremum, y_extremum]
end


