
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


"""
Get hyperuniformity metric which is the structure factor
at zero momentum normalized by the height of the first peak in the structure factor
as defined in equation 251 in 10.1016/j.physrep.2018.03.001
"""
function get_hyperuniformity_metric(structure_factor_dict::Dict)

    # locate first peak of structure factor
    pks, vals = Peaks.findmaxima(structure_factor_dict["structure_factor_vec"])

    # cut structure factor data at momentum just above first peak
    structure_factor_cut_vec = structure_factor_dict["structure_factor_vec"][1:pks[2]+1]
    wavenumber_cut_vec = structure_factor_dict["wavenumber_vec"][1:pks[2]+1]

    # set the order of the fitted polynomial
    polynomial_order = 5

    # fit polynomial of given order to cut data
    polynomial_fit = Polynomials.fit(wavenumber_cut_vec, 
                                    structure_factor_cut_vec,
                                    polynomial_order)

    # get extrapolated structure factor at zero momentum
    structure_factor_zero_momentum = polynomial_fit(0)

    # get critical momenta which is roots of first derivative of polynomial
    polynomial_derivative = Polynomials.derivative(polynomial_fit)
    critical_momenta = Polynomials.roots(polynomial_derivative)

    # get real critical momenta
    critical_momenta_real = real.(critical_momenta[imag.(critical_momenta) .== 0])

    # get fitted structure factor at highest peak
    structure_factor_first_peak = maximum( polynomial_fit.(critical_momenta_real) )

    # get hyperuniformity metric
    hyperuniformity_metric = structure_factor_zero_momentum/structure_factor_first_peak

    return [hyperuniformity_metric, polynomial_fit]
end



"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density_array_by_fft(complete_autocovariance_fct_direction_dict::Dict;
            save_result = false,
            save_path = raw"..\analysis_data\sample_name",
            voxel_edge_length = nothing,
            label = nothing)

    # get values of autocovariance function
    complete_autocovariance_fct_array_values = Measurements.value.( 
                complete_autocovariance_fct_direction_dict["autocovariance_fct_array"] )

    # determine fourier transform of autocovariance function values
    spectral_density_array_fft_output = FFTW.rfft(complete_autocovariance_fct_array_values)

    # bring spectral densities into an order from negative to positive wavenumbers 
    spectral_density_array = FFTW.fftshift(spectral_density_array_fft_output, [2, 3])

    # get tuple of vectors of sampled wavenumbers
    wavenumber_vec_vec = get_wavenumber_vec_vec(complete_autocovariance_fct_array_values)

    # get array of sampled sampled wavevectors 
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    # create dict to save
    spectral_density_dict = Dict("wavevector_array" => wavevector_array,
        "wavenumber_vec_vec" => wavenumber_vec_vec,
        "spectral_density_array" => spectral_density_array,
        "voxel_edge_length" => complete_autocovariance_fct_direction_dict["voxel_edge_length"],
        "label" => complete_autocovariance_fct_direction_dict["label"])

    # if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict), save_path*"_spectral_density_array.h5")

    end

    return spectral_density_dict
    
end



"""
determine spectral density as a function of wavevector.
Unfortunately in the moment this function takes forever (weeks?).
Instead, I think I need to use the FFT
"""
function get_spectral_density_by_wavevector_array(
                structure_dict::Dict;
                nr_measurements_per_direction::Int64 = 1000,
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(structure_dict;
                                            sampling_vec_array = sampling_vec_array,
                                            nr_measurements_per_direction = nr_measurements_per_direction),
                save_result = false,
                save_path = raw"..\analysis_data\sample_name",
                voxel_edge_length = nothing,
                label = nothing)

    # correct the sampling distance vec along the third dimension, where due to the mirror
    # symmetry of the autocovariance fct only positive z values where considered
    corrected_sampling_distance_vec_vec = [
                            autocovariance_fct_direction_dict["sampling_distance_vec_vec"][1], 
                            autocovariance_fct_direction_dict["sampling_distance_vec_vec"][2],
                            vcat(.- reverse(autocovariance_fct_direction_dict["sampling_distance_vec_vec"][3][2:end]),
                                autocovariance_fct_direction_dict["sampling_distance_vec_vec"][3]) ]
    
    # get vector of sampled wavenumbers along the three coordinate directions
    wavenumber_vec_vec = []

    for i in 1:3

        push!(wavenumber_vec_vec, get_wavenumber_vec(corrected_sampling_distance_vec_vec[i]) )

    end

    # create array of wavevectors
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    # determine size of spectral density array
    spectral_density_array_size = size(wavevector_array)[1:3]

    # initialize spectral density array 
    spectral_density_array = Array{Complex{Measurements.Measurement}}(undef, spectral_density_array_size...)

    # for each wavevector, determine the spectral density
    for i in eachindex(wavenumber_vec_vec[1])
        for j in eachindex(wavenumber_vec_vec[2])
            for k in eachindex(wavenumber_vec_vec[3])

                spectral_density_array[i,j,k] = get_spectral_density(
                                                    wavevector_array[i,j,k,:], 
                                                    autocovariance_fct_direction_dict["sampling_distance_vec_vec"], 
                                                    autocovariance_fct_direction_dict["autocovariance_fct_array"])

            end
        end

        println("wavenumber "*string(wavenumber_vec_vec[1][i]*" along x done") )
    end

    # create dict to save
    spectral_density_dict = Dict("wavevector_array" => wavevector_array,
                                "wavenumber_vec_vec" => wavenumber_vec_vec,
                                "spectral_density_array" => spectral_density_array,
                                "nr_measurements_per_direction" 
                                        => autocovariance_fct_direction_dict["nr_measurements_per_direction"],
                                "voxel_edge_length" => autocovariance_fct_direction_dict["voxel_edge_length"],
                                "label" => autocovariance_fct_direction_dict["label"])

    # if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict),
            save_path*"_spectral_density_array.h5")

    end

    return spectral_density_dict

end



"""
Calculate spectral density from autocovariance fct along a given direction.
Somehow this function does not work yet, maybe because the calculation of exponentials for
Measurement type data simply takes too long.
"""
function get_spectral_density_along_direction_by_wavenumber_vec(structure_dict::Dict,
                direction_vec::Vector;
                nr_measurements_per_direction::Int64 = 1000,
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(structure_dict;
                                            sampling_vec_array = sampling_vec_array,
                                            nr_measurements_per_direction = nr_measurements_per_direction),
                save_result = false,
                save_path = raw"..\analysis_data\sample_name_direction",
                voxel_edge_length = nothing,
                label = nothing)
    
    # extract autocovariance function vector along given direction
    autocovariance_fct_along_direction_vec = get_autocovariance_fct_along_direction_vec(direction_vec, 
                                                                autocovariance_fct_direction_dict["autocovariance_fct_array"])

    # get sampling distances along given direction
    sampling_distance_along_direction_vec = get_sampling_distance_along_direction_vec(direction_vec, 
                                                                autocovariance_fct_along_direction_vec)

    # get vector of wavenumbers where spectral density will be calculated
    wavenumber_vec = get_wavenumber_vec(sampling_distance_along_direction_vec)

    # initialize spectral density vector
    spectral_density_vec = Vector{Measurements.Measurement}(undef, length(wavenumber_vec))

    # for each wavenumber, determine spectral density
    for i in eachindex(wavenumber_vec)
        spectral_density_vec[i] = get_spectral_density_along_direction(wavenumber_vec[i], 
                                            sampling_distance_along_direction_vec, 
                                            autocovariance_fct_along_direction_vec)

    end

    # create dict to save
    spectral_density_dict = Dict("wavenumber_vec" => wavenumber_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "nr_measurements_per_direction" => nr_measurements_per_direction,
                        "direction_vec" => direction_vec,
                        "voxel_edge_length" => autocovariance_fct_direction_dict["voxel_edge_length"],
                        "label" => autocovariance_fct_direction_dict["label"])

    # if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict),
            save_path*"_spectral_density_direction.h5")

    end

    return spectral_density_dict
    
end


"""
determine spectral density as a function of wavevector.
Unfortunately this function takes very long (couple of hours for a single wavevector)
"""
function get_spectral_density(wavevector::Vector{Float64}, 
                                sampling_distance_vec_vec::Vector, 
                                autocovariance_fct_array::Array)

    # initialize complex spectral density
    spectral_density = 0.0 + 0.0 * im

    # create an autocovariance function array that is mirrored in the first two dimensions
    autocovariance_fct_array_mirrored = reverse( reverse(autocovariance_fct_array, dims=1), dims=2) 

    for i in eachindex(sampling_distance_vec_vec[1])
        for j in eachindex(sampling_distance_vec_vec[2])

            # initialize complex spectral density sum along z direction
            spectral_density_z_component_sum = 0.0 + 0.0 * im

            for k in eachindex(sampling_distance_vec_vec[3])[2:end]

                spectral_density_z_component_sum += (exp(im*wavevector[3]*sampling_distance_vec_vec[3][k])
                                                            *autocovariance_fct_array_mirrored[i,j,k]
                                                    + exp(-im*wavevector[3]*sampling_distance_vec_vec[3][k])
                                                            *autocovariance_fct_array[i,j,k]   )
                                    
            end

            # add sum along z direction and other terms to sum along x and y directions
            spectral_density += (exp(-im*wavevector[1]*sampling_distance_vec_vec[1][i]) 
                                    * exp(-im*wavevector[2]*sampling_distance_vec_vec[2][j]) 
                                    * ( autocovariance_fct_array[i,j,1] + spectral_density_z_component_sum )
                                    )

            println("sampling distance "*string(sampling_distance_vec_vec[2][j])*" along y done")
        end

        println("sampling distance "*string(sampling_distance_vec_vec[1][i])*" along x done")
    end

    # multiply by the inverse of the sampling volume
    spectral_density *= 1/( (sampling_distance_vec_vec[1][2] - sampling_distance_vec_vec[1][1])
                            * (sampling_distance_vec_vec[2][2] - sampling_distance_vec_vec[2][1])
                            * (sampling_distance_vec_vec[3][2] - sampling_distance_vec_vec[3][1]) )


    return spectral_density
end



"""
get spectral density 
"""
function get_spectral_density_along_direction(wavenumber::Float64, 
                            sampling_distance_vec::Vector, 
                            autocovariance_fct_along_direction_vec::Vector)

    # determine sampling distance
    sampling_distance = sampling_distance_vec[2] - sampling_distance_vec[1]

    # calculate fourier transform for given wavenumber
    spectral_density = (1/sampling_distance 
                            * ( autocovariance_fct_along_direction_vec[1] 
                                + 2 * sum( autocovariance_fct_along_direction_vec[2:end] 
                                                .* cos.( wavenumber .* sampling_distance_vec[2:end] )  ) 
                                )
                        )

    return spectral_density

end


"""
get vector of sampling distances along the given direction
"""
function get_sampling_distance_along_direction_vec(direction_vec, 
                autocovariance_fct_along_direction_vec::Vector)
    
    # determine geometrical length of direction vector
    sampling_distance = sqrt(sum( direction_vec .^2 ))

    # get vector of sampling distances
    sampling_distance_vec = (collect(0:length(autocovariance_fct_along_direction_vec)-1) 
                                .* sampling_distance)

    return sampling_distance_vec

end



"""
function to calculate and save the following measures for a given 3d data set
- local volume fraction variance
- autocovariance function as a function of sampling distance (assuming an isotropic medium)
- spectral density as a function of sampling distance (assuming an isotropic medium)
- autocovariance function as a function of sampling vector (not assuming an isotropic medium)
- spectral density along 6 different directions (not assuming an isotropic medium)
"""
function save_statistical_measures(data_path::String,
                                    save_path::String;
                                    voxel_edge_length = nothing, 
                                    label = nothing,
                                    nr_sampling_distances = nothing,
                                    nr_measurements_per_distance = 10000,
                                    nr_window_sizes = 100,
                                    nr_measurements_per_direction = 1000,
                                    save_autocovariance_fct = true,
                                    save_spectral_density = true,
                                    save_local_volume_fraction_variance = true,
                                    save_autocovariance_fct_direction = true,
                                    save_spectral_density_along_directions = true,
                                    save_complete_autocovariance_fct_direction = true,
                                    save_spectral_density_array = true)

    # load structure dictionary which contains all essential information about the structure
    structure_dict = GU.load_h5_dict(data_path)

    # set voxel edge length and label from structure dict if not specified in the arguments
    if voxel_edge_length === nothing
        voxel_edge_length = structure_dict["voxel_edge_length"]
    end
    if label === nothing
        label = structure_dict["label"]
    end

    # if not specified, determine nr of sampling distances for autocovariance function
    # and spectral density
    if nr_sampling_distances === nothing
        nr_sampling_distances = get_nr_sampling_distances(structure_dict["mean_edge_length_data"] )
    end

    # save autocovariance function if desired
    if save_autocovariance_fct

        # get autocovariance function as a function of sampling distance
        autocovariance_fct_isotrope_dict = get_autocovariance_fct_isotrope_by_sampling_distance_vec(
                                                structure_dict;
                                                nr_sampling_distances = nr_sampling_distances,
                                                nr_measurements_per_distance = nr_measurements_per_distance,
                                                save_result = true,
                                                save_path = save_path,
                                                voxel_edge_length = voxel_edge_length,
                                                label = label)

    end


    # save spectral density if desired
    if save_spectral_density

        # check if autocovariance fct was calculated within this function
        if !save_autocovariance_fct

            # if autocovariance function was not calculated within this function, check if it was
            # calculated before and if so, load the corresponding dictionary
            if isfile(save_path*"_autocovariance_fct.h5")

                # load autocovariance function dict
                autocovariance_fct_isotrope_dict = GU.load_h5_dict(save_path*"_autocovariance_fct.h5")

            else
                autocovariance_fct_isotrope_dict = get_autocovariance_fct_isotrope_by_sampling_distance_vec(
                                                structure_dict;
                                                nr_sampling_distances = nr_sampling_distances,
                                                nr_measurements_per_distance = nr_measurements_per_distance,
                                                save_result = false,
                                                save_path = save_path,
                                                voxel_edge_length = voxel_edge_length,
                                                label = label)
            end
        end
    
        # calculate spectral_density
        spectral_density_isotrope_dict = get_spectral_density_isotrope_by_wavenumber_vec(
                structure_dict;
                nr_sampling_distances = length(autocovariance_fct_isotrope_dict["sampling_distance_vec"]),
                nr_measurements_per_distance = autocovariance_fct_isotrope_dict["nr_measurements_per_distance"],
                sampling_distance_vec = autocovariance_fct_isotrope_dict["sampling_distance_vec"],
                autocovariance_fct_dict = autocovariance_fct_isotrope_dict,
                save_result = true,
                save_path = save_path,
                voxel_edge_length = voxel_edge_length,
                label = label)
                    
    end


    # save local volume fraction if desired
    if save_local_volume_fraction_variance

        # determine local volume fraction variance vector by using measuring windows
        local_volume_fract_variance_dict = get_local_volume_fract_variance_by_window_vec(
                                                            structure_dict;
                                                            nr_window_sizes = nr_window_sizes,
                                                            window_positioning="random",
                                                            window_shape="spherical",
                                                            save_result = true,
                                                            save_path = save_path,
                                                            voxel_edge_length = voxel_edge_length,
                                                            label = label)
    
    end


    # save autocovariance function as a function of sampling direction
    if save_autocovariance_fct_direction

        # get autocovariance function array
        autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                            structure_dict;
                            nr_measurements_per_direction = nr_measurements_per_direction,
                            save_result = true,
                            save_path = save_path,
                            voxel_edge_length = voxel_edge_length,
                            label = label)
    
    end


    # save spectral density along certain directions if desired
    if save_spectral_density_along_directions

        # set vectors along which spectral density will be measured
        direction_vec_vec = [[1,0,0], 
                            [0,1,0],
                            [0,0,1],
                            (1/sqrt(2)) .* [1,-1,0],
                            (1/sqrt(3)) .* [1,1,1],
                            (1/sqrt(6)) .* [1,1,-2]]

        # set vector with names to save the files
        naming_vec = string.( [[1,0,0], 
                                [0,1,0],
                                [0,0,1],
                                [1,-1,0],
                                [1,1,1],
                                [1,1,-2]] )

        # check if autocovariance function was previously calculated in this function
        if !save_autocovariance_fct_direction

            # check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction.h5" )
            
                # load autocovariance function per direction dict
                autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
            
            else
                # get autocovariance function array
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                    structure_dict;
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    save_result = false,
                    save_path = save_path,
                    voxel_edge_length = voxel_edge_length,
                    label = label)

            end
        end

        for i in eachindex(direction_vec_vec)

            # determine spectral density
            spectral_density_along_direction_dict = get_spectral_density_along_direction_by_wavenumber_vec(
                    structure_dict,
                    direction_vec_vec[i];
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    sampling_distance_vec_vec = autocovariance_fct_direction_dict["sampling_distance_vec_vec"],
                    autocovariance_fct_dict = autocovariance_fct_direction_dict,
                    save_result = true,
                    save_path = save_path*"_"*naming_vec[i],
                    voxel_edge_length = voxel_edge_length,
                    label = label*" "*naming_vec[i])


        end
    end


    # save complete autocovariance function as a function of sampling direction for all spacial directions,
    # not only the half space considered previously
    if save_complete_autocovariance_fct_direction

        # check if autocovariance function was previously calculated in this function
        if !save_autocovariance_fct_direction

            # check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction.h5" )
            
                # load autocovariance function per direction dict
                autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
            
            else
                # get autocovariance function array
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                    structure_dict;
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    save_result = false,
                    save_path = save_path,
                    voxel_edge_length = voxel_edge_length,
                    label = label)

            end
        end
    
        # calculate complete autocovariance function
        complete_autocovariance_fct_direction_dict = get_complete_autocovariance_fct_by_sampling_vec_array(
            autocovariance_fct_direction_dict;
            save_result = true,
            save_path = save_path)
                    
    end


    # save complete autocovariance function as a function of sampling direction for all spacial directions,
    # not only the half space considered previously
    if save_spectral_density_array

        # check if complete autocovariance function was previously calculated in this function
        if save_complete_autocovariance_fct_direction

            # check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction_complete.h5" )
            
                # load complete autocovariance function per direction dict
                complete_autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_complete.h5")
            
            else

                # check if autocovariance function was previously calculated in this function
                if !save_autocovariance_fct_direction
                
                    # check if data can be loaded from file
                    if isfile( save_path*"_autocovariance_fct_direction.h5" )
                    
                        # load autocovariance function per direction dict
                        autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
                    
                    else
                        # get autocovariance function array
                        autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                            structure_dict;
                            nr_measurements_per_direction = nr_measurements_per_direction,
                            save_result = false,
                            save_path = save_path,
                            voxel_edge_length = voxel_edge_length,
                            label = label)
                    
                    end
                end

                # get complete autocovariance function array
                complete_autocovariance_fct_direction_dict = get_complete_autocovariance_fct_by_sampling_vec_array(
                    autocovariance_fct_direction_dict;
                    save_result = false,
                    save_path = save_path)

            end
        end
    
        # calculate spectral_density
        spectral_density_array_dict = get_spectral_density_array_by_fft(complete_autocovariance_fct_direction_dict;
                                                                    save_result = true,
                                                                    save_path = save_path,
                                                                    voxel_edge_length = voxel_edge_length,
                                                                    label = label)

        # calculate spectral_density
        spectral_density_array_dict = get_spectral_density_by_wavevector_array_fft(structure_dict;
            nr_measurements_per_direction = nr_measurements_per_direction,
            sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
            sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
            autocovariance_fct_direction_dict = autocovariance_fct_direction_dict,
            save_complete_autocovariance_fct_direction = save_complete_autocovariance_fct_direction,
            complete_autocovariance_fct_direction_dict = complete_autocovariance_fct_direction_dict,
            save_result = true,
            save_path = save_path)
                    
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

    # if desired plot autocovariance function
    if plot_autocovariance_fct_bool

        # initialize vector for plot dicts
        autocovariance_fct_plot_dict_vec = []

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            autocovariance_fct_plot_dict = GU.load_h5_dict(data_path_vec[i]*"_autocovariance_fct.h5")

            # if desired adjust label and voxel edge length
            if  label_vec !== nothing
                autocovariance_fct_plot_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                autocovariance_fct_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # add current plot dict to vector
            push!(autocovariance_fct_plot_dict_vec, autocovariance_fct_plot_dict)

        end

        # plot the data
        plot_autocovariance_fct(autocovariance_fct_plot_dict_vec,
                        save_path;
                        title="Autocovariance function",
                        save_plot = true)

        # plot zoom into the data
        plot_autocovariance_fct(autocovariance_fct_plot_dict_vec,
                        save_path;
                        title="Autocovariance function",
                        save_plot = true,
                        ylims = [-0.07, 0.07])

    end

    # if desired plot spectral density
    if plot_spectral_density_bool

        # initialize vector for plot dicts
        spectral_density_plot_dict_vec = []

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            spectral_density_plot_dict = GU.load_h5_dict(data_path_vec[i]*"_spectral_density.h5")

            # if desired adjust label and voxel edge length
            if label_vec !== nothing
                spectral_density_plot_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if  voxel_edge_length_vec !== nothing
                spectral_density_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # add current plot dict to vector
            push!(spectral_density_plot_dict_vec, spectral_density_plot_dict)

        end

        # plot the data
        plot_spectral_density(spectral_density_plot_dict_vec,
                        save_path;
                        save_plot = true)

        # if specified, plot a zoom into the spectral sensity
        if spectral_density_xlims !== nothing

            plot_spectral_density(spectral_density_plot_dict_vec,
                            save_path;
                            save_plot = true,
                            xlims = spectral_density_xlims)
        end

    end


    # if desired plot local volume fraction variance
    if plot_local_volume_fraction_variance_bool

        # initialize vector for plot dicts
        local_volume_fraction_variance_plot_dict_vec = []

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            local_volume_fraction_variance_plot_dict = GU.load_h5_dict(
                                                            data_path_vec[i]*"_volume_fraction_variance.h5")

            # if desired adjust label and voxel edge length
            if  label_vec !== nothing
                local_volume_fraction_variance_plot_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                local_volume_fraction_variance_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # add current plot dict to vector
            push!(local_volume_fraction_variance_plot_dict_vec, local_volume_fraction_variance_plot_dict)

        end

        # plot the data
        plot_volume_fraction_variance(local_volume_fraction_variance_plot_dict_vec,
                        save_path;
                        save_plot = true)

    end


    # if desired plot autocovariance function heatmaps
    if plot_autocovariance_fct_direction_bool

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            autocovariance_fct_direction_plot_dict = GU.load_h5_dict(
                                                            data_path_vec[i]*"_autocovariance_fct_direction_complete.h5")

            # if desired adjust label and voxel edge length
            if label_vec  !== nothing
                autocovariance_fct_direction_plot_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                autocovariance_fct_direction_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # plot autocovariance function heatmaps along x-y, x-z and y-z planes
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


    # if desired plot spectral density along three different directions
    if plot_spectral_density_along_directions_bool

        # set vector with names to load the files
        naming_vec = string.( [[1,0,0], 
                                [0,1,0],
                                [0,0,1],
                                [1,-1,0],
                                [1,1,1],
                                [1,1,-2]] )

        # set vector with the two collections of directions that will be plotted in two different plots
        range_vec = [1:3, 4:6]

        # set the names of the two plots
        plot_name_vec = ["_along_axes", "_rotated_axes"]

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # create two plots per sample, one along the axes and one along rotated coordinate axes
            for plot_per_sample_nr in 1:2

                # initialize vector for plot dicts
                spectral_density_direction_plot_dict_vec = []

                spectral_density_direction_plot_dict = Dict()

                # loop through the three directions either along the axes or along rotated coordinate axes
                for j in range_vec[plot_per_sample_nr]

                    # load plot dictionary
                    spectral_density_direction_plot_dict = GU.load_h5_dict(
                                                data_path_vec[i]*"_"*naming_vec[j]*"_spectral_density_direction.h5")

                    # if desired adjust label
                    if label_vec  !== nothing
                        spectral_density_direction_plot_dict["label"] = label_vec[i]*" "*naming_vec[j]
                    end

                    # if desired adjust voxel edge length
                    if voxel_edge_length_vec !== nothing
                        spectral_density_direction_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
                    end

                    # add current plot dict to vector
                    push!(spectral_density_direction_plot_dict_vec, spectral_density_direction_plot_dict)

                end

                # set saving path by extracting the current sample name from data path
                save_path_specific = ( first( save_path, findlast('\\', save_path) )
                                *SubString(data_path_vec[i], (findlast('\\', data_path_vec[i]) + 1), length(data_path_vec[i]))
                                *"_direction"
                                *plot_name_vec[plot_per_sample_nr]
                                )

                # plot the data
                plot_spectral_density(spectral_density_direction_plot_dict_vec,
                                        save_path_specific;
                                        save_plot = true)

                # if specified, plot a zoom into the spectral sensity
                if spectral_density_xlims !== nothing
                
                    plot_spectral_density(spectral_density_direction_plot_dict_vec,
                    save_path_specific;
                                save_plot = true,
                                xlims = spectral_density_xlims)
                end

            end

        end

    end


    # if desired plot spectral density heatmaps
    if plot_spectral_density_heatmaps_bool

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            spectral_density_array_dict = GU.load_h5_dict(data_path_vec[i]*"_spectral_density_array.h5")

            # if desired adjust label and voxel edge length
            if label_vec  !== nothing
                spectral_density_array_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                spectral_density_array_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # plot spectral density heatmaps along x-y, x-z and y-z planes
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



"""
Get the autocovariance function for 3d media with periodic
boundary conditions
"""
function get_autocovariance_fct(sampling_vec::Vector{Int64},
    structure_dict::Dict)

    # initialize vector from which the two point prob. fct. will be calculated later
    two_point_prob_fct_summand_vec = Vector{Float64}(undef, prod(structure_dict["size_data"]) )
    current_index = 1

    # loop through all voxels
    for i in 1:structure_dict["size_data"][1]
        for j in 1:structure_dict["size_data"][2]
            for k in 1:structure_dict["size_data"][3]

                # get indices of current voxel and the one at given sampling vector to it
                x1 = (i,j,k)
                x2 = (mod.(x1 .+ sampling_vec .- 1, structure_dict["size_data"]) .+ 1)

                # calculate the contribution to the two point prob. fct. from these coodinates
                two_point_prob_fct_summand_vec[current_index] = structure_dict["data_binary"][x1...] * structure_dict["data_binary"][x2...]

                current_index += 1

            end
        end
    end

    # calculate 2 point prob. function
    two_point_prob_fct = Statistics.mean( two_point_prob_fct_summand_vec )

    # determine autocovariance function
    autocovariance_fct = two_point_prob_fct - structure_dict["volume_fract_tot"]^2

    return autocovariance_fct
end



"""
Fully relax a cluster of vertices. The cluster energy will always be updated
"""
function relax_cluster_keating!(graph_dict::Dict,
    cluster_dict::Dict, 
    evolution_dict::Dict;
    threshold_cluster_energy = Inf,
    update_total_energy::Bool = false,
    print_progress::Bool = false)

    # make sure that cluster energy is up to date
    if !cluster_dict["cluster_energy_up_to_date"]
        cluster_dict["cluster_energy"] = get_cluster_energy(graph_dict, cluster_dict)
    end

    # store initial cluster energy
    unrelaxed_cluster_energy  = cluster_dict["cluster_energy"]

    # make sure that total energy is up to date if it will be updated later
    if !graph_dict["total_energy_up_to_date"] && update_total_energy
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
    end

    # if network is supposed to be relaxed globally, store initial shell nr
    # and set threshold cycle for global relaxation
    if haskey(evolution_dict, "relax_globally_after_threshold_cycle")
        initial_shell_nr = evolution_dict["shell_nr"]
        threshold_cycle_global_relaxation = (evolution_dict["reject_during_relaxation_cycle_threshold"]*2+1)
    else
        threshold_cycle_global_relaxation = evolution_dict["nr_max_relaxation_cycles"] + 1
    end

    # perform the given number of relaxation cycles
    for cycle_nr in 1:evolution_dict["nr_max_relaxation_cycles"]

        # from the threshold cylcle on, relax globally, if desired
        if cylce_nr == threshold_cycle_global_relaxation
            # the following gives shell_nr = 8 for 216 vertices
            evolution_dict["shell_nr"] = Int(ceil(log(graph_dict["nr_vertices"])))+2
        end

        # only update cluster energy, if this is needed to get cluster energy change
        if cycle_nr <= evolution_dict["reject_during_relaxation_cycle_threshold"]-1
            update_cluster_energy = false
        else
            update_cluster_energy = true

            # store previous cluster force and energy
            previous_cluster_force = cluster_dict["cluster_force"]
            previous_cluster_energy = cluster_dict["cluster_energy"]
        end

        # relax cluster for one cycle
        graph_dict, cluster_dict = relax_cluster_one_cycle_keating!(graph_dict, 
        cluster_dict,
        evolution_dict;
        update_total_energy = false,
        update_cluster_energy = update_cluster_energy )

        # if cycle nr is above the given threshold, check if the relaxation can 
        # be breaked when it becomes clear that the total energy will exceed the threshold
        # or because of small relative energy change
        if cycle_nr > evolution_dict["reject_during_relaxation_cycle_threshold"]

            # get vector of last two cluster forces
            cluster_force_vec = [previous_cluster_force, cluster_dict["cluster_force"]]

            # get vector of last two cluster energies
            cluster_energy_vec =[previous_cluster_energy, cluster_dict["cluster_energy"]]

            # estimate relaxed cluster energy
            prefactor_force_squared, relaxed_cluster_energy = get_energy_relaxation_coefficients(
                cluster_force_vec, cluster_energy_vec)

            if print_progress
                println("Prefactor force squared: "*string(prefactor_force_squared))
                println("Relaxed cluster energy: "*string(relaxed_cluster_energy))
            end

            # break if estimated energy change exceeds the given threshold
            if (prefactor_force_squared < 1
                && relaxed_cluster_energy > 1.05*threshold_cluster_energy) 
                
                if print_progress
                    println("Relaxed energy exceeds threshold: breaking at cycle nr "*string(cycle_nr))
                end
                break
            end

            # break if cluster energy changes less than the given threshold
            relative_cluster_energy_change = (
                abs((previous_cluster_energy - cluster_dict["cluster_energy"])
                        /cluster_dict["cluster_energy"]))

            if relative_cluster_energy_change < evolution_dict["break_at_relative_cluster_energy_change"] 
    
                if print_progress
                    println("Negligeable energy change: breaking at cycle nr "*string(cycle_nr))
                end
                break
            end

        end

    end

    # restore initial shell nr in case it was altered during relaxation
    if haskey(evolution_dict, "relax_globally_after_threshold_cycle")
        evolution_dict["shell_nr"] = initial_shell_nr
    end

    # update total energy if desired
    if update_total_energy
        graph_dict["total_energy"] = (graph_dict["total_energy"] 
                                    + cluster_dict["cluster_energy"]
                                    - unrelaxed_cluster_energy)

        graph_dict["total_energy_up_to_date"] = true
    else
        graph_dict["total_energy_up_to_date"] = false
    end
    
    return [graph_dict, cluster_dict]
end


"""
Perform a Monte Carlo move without thermal fluctuations
by switching a bond, relaxing the network and then accepting the network
with Metropolis acceptance probability
"""
function monte_carlo_move!(graph_dict::Dict, 
    evolution_dict::Dict,
    temperature; 
    switched_chain::Tuple{Int64, Int64, Int64, Int64} = get_random_chain(graph_dict),
    print_progress::Bool = false)

    # save original graph dict 
    initial_graph_dict = deepcopy(graph_dict)

    # get initial cluster before bond switch
    initial_cluster_dict = get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_chain; 
                                    shell_nr = evolution_dict["shell_nr"])

    # make sure that total energy is up to date
    if !graph_dict["total_energy_up_to_date"]
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
        graph_dict["total_energy_up_to_date"] = true
    end

    # initialize cluster weights which will contribute to the acceptance probability
    # when thermal fluctuations are included
    cluster_relaxation_weight = 1
    cluster_excitation_weight = 1

    # set threshold cluster energy for Metropolis acceptance probability
    # if there are no thermal fluctuations considered
    if !evolution_dict["thermal_fluctuations"]
        threshold_cluster_energy = (initial_cluster_dict["cluster_energy"] 
                                    - temperature * log(rand()))
    else
        threshold_cluster_energy = Inf
    end

    # if there are thermal fluctuations, relax cluster first and calculate
    # weights of the corresponding shifts
    if evolution_dict["thermal_fluctuations"]

        # deep copy initial cluster, such that it does not get modified
        cluster_dict = deepcopy(initial_cluster_dict)

        # relax cluster
        graph_dict, cluster_dict = relax_cluster_keating!(graph_dict,
            cluster_dict,
            evolution_dict;
            update_total_energy = false,
            print_progress = print_progress)

        # get cluster weight corresponding to the relaxation translations
        cluster_relaxation_weight = get_cluster_fluctuation_weight(initial_graph_dict, 
                                                            graph_dict, 
                                                            cluster_dict,
                                                            temperature)
    end

    # switch bonds
    graph_dict = switch_chain!(graph_dict, switched_chain)

    # get cluster after bond switch
    cluster_dict = get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_chain; 
                                    shell_nr = evolution_dict["shell_nr"])

    # relax cluster around switched chain and only update energy when there won't be
    # thermal fluctuations included afterward
    graph_dict, cluster_dict = relax_cluster_keating!(graph_dict,
        cluster_dict,
        evolution_dict;
        threshold_cluster_energy = threshold_cluster_energy,
        update_total_energy = false)


    # if desired, include thermal fluctuations by randomly shifting all cluster vertices
    if evolution_dict["thermal_fluctuations"]

        # excite cluster with thermal fluctuations and get the corresponding excitation
        # weight
        graph_dict, cluster_dict, cluster_excitation_weight = excite_cluster!(graph_dict,
                                                            cluster_dict,
                                                            temperature;
                                                            update_total_energy = false,
                                                            update_cluster_energy = true)

        # set random threshold energy increase for Metropolis acceptance probability
        threshold_cluster_energy = (initial_cluster_dict["cluster_energy"]
        - temperature 
        * (cluster_excitation_weight/cluster_relaxation_weight) * log(rand()))

    end

    # accept move if energy increase is below threshold
    move_accepted = false

    if cluster_dict["cluster_energy"] <= threshold_cluster_energy
        # update total energy
        graph_dict["total_energy"] = (graph_dict["total_energy"] 
        + cluster_dict["cluster_energy"] - initial_cluster_dict["cluster_energy"])  
        graph_dict["total_energy_up_to_date"] = true   

        move_accepted = true
    else
        graph_dict = initial_graph_dict
    end

    return [graph_dict, move_accepted]
end



"""
Define a anisotropy metric as the normalized variance of the structure factor at the
first peak of the structure factor of the diamond lattice
"""
function get_anisotropy_metric_from_structure_factor(
    structure_factor_angle_averaged_dict::Dict,
    diamond_structure_factor_peak_wavenumber::Float64 = 4.676810,
    diamond_structure_factor_peak_std::Float64 = 8.5573)

    # find the wavenumber that is the closest to the diamond peak wavenumber
    diamond_structure_factor_peak_wavenumber_index = argmin(abs.(
        structure_factor_angle_averaged_dict["wavenumber_vec"] 
        .- diamond_structure_factor_peak_wavenumber))

    # get the structure factor standard deviation around the diamond peak
    peak_structure_factor_std = Measurements.uncertainty(
        structure_factor_angle_averaged_dict["structure_factor_vec"][diamond_structure_factor_peak_wavenumber_index])

    # normalize this standard deviation by the structure factor standard deviation of the diamond peak
    anisotropy_metric_from_structure_factor = peak_structure_factor_std / diamond_structure_factor_peak_std

    return anisotropy_metric_from_structure_factor
end



"""
Define a anisotropy metric as the normalized variance of the spectral density at the
first peak of the spectral density of the diamond lattice
"""
function get_anisotropy_metric_from_spectral_density(
    spectral_density_angle_averaged_dict::Dict,
    diamond_spectral_density_peak_wavenumber::Float64 = 4.680517,
    diamond_spectral_density_peak_std::Float64 = 521.88398)

    # find the wavenumber that is the closest to the diamond peak wavenumber
    diamond_spectral_density_peak_wavenumber_index = argmin(abs.(spectral_density_angle_averaged_dict["wavenumber_vec"] .- diamond_spectral_density_peak_wavenumber))

    # get the spectral density standard deviation around the diamond peak
    peak_spectral_density_std = Measurements.uncertainty(
        spectral_density_angle_averaged_dict["spectral_density_vec"][diamond_spectral_density_peak_wavenumber_index])

    # normalize this standard deviation by the spectral density standard deviation of the diamond peak
    anisotropy_metric_from_spectral_density = peak_spectral_density_std / diamond_spectral_density_peak_std

    return anisotropy_metric_from_spectral_density
end





"""
Calculate the pore size distribution following the method described in
10.1103/PhysRevE.100.053314, modified to work with periodic boundary conditions
and by using random sampling of voxels to speed up the calculation
"""
function get_pore_size_distribution_old(structure_dict::Dict;
    nr_sampled_voxels::Int = 20000,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing,
    print_progress::Bool = false,
    thread_nr::Int64 = 0,
    print_lock = Threads.ReentrantLock())

    # create list of digital spheres with increasing radius
    sphere_pixel_radius_vec = collect(0.5001:0.5:minimum(structure_dict["size_data"]./2))
    digital_sphere_list = [get_digital_sphere(radius) for radius in sphere_pixel_radius_vec]

    # create array with pore radii
    pore_pixel_radius_array = zeros(size(structure_dict["data_binary"])...)

    # initialize counter if progress is printed
    if print_progress
        i = 1
    end

    # sample given number of voxels to speed up calculation
    sampled_coords = StatsBase.sample(CartesianIndices(structure_dict["data_binary"]), 
        nr_sampled_voxels, replace=false)

    # loop through all sampled voxels using cartesian indices
    for coord in sampled_coords

        # check if voxel is in pore
        if !structure_dict["data_binary"][coord]

            # loop through all digital spheres
            for j in eachindex(digital_sphere_list)

                # get digital sphere
                digital_sphere = digital_sphere_list[j]

                # check if voxel is in digital sphere
                if all([!structure_dict["data_binary"][
                    ((coord.I .- [1,1,1] .+ minimum(structure_dict["size_data"]) .+ digital_sphere[k]
                        ) .% structure_dict["size_data"] .+ [1,1,1])...] 
                        for k in eachindex(digital_sphere)])

                    # set pore pixel radius to maximum of current and previous radius
                    for k in eachindex(digital_sphere)

                        pore_pixel_radius_array[
                            ((coord.I .- [1,1,1] .+ minimum(structure_dict["size_data"]) .+ digital_sphere[k]
                        ).% structure_dict["size_data"] .+ [1,1,1])...
                            ] = maximum([
                                pore_pixel_radius_array[
                                    ((coord.I .- [1,1,1] .+ minimum(structure_dict["size_data"]) .+ digital_sphere[k]
                        ).% structure_dict["size_data"] .+ [1,1,1])...],
                                        sphere_pixel_radius_vec[j] 
                                ])
                    end
                
                else
                    break

                end

            end
        end

        # print progress
        if print_progress
            progress_percentage = i/nr_sampled_voxels*100

            lock(print_lock) do
                Format.printfmtln("Pore size distribution calculation progress thread nr {1:d}: {2:.1f} %", 
                    thread_nr, progress_percentage)
            end

            i += 1
        end

    end

    # shape the pore pixel radius array into a vector
    pore_pixel_radius_vec = vec(pore_pixel_radius_array)

    # filter out voxels that are not in a pore
    pore_pixel_radius_filtered_vec = pore_pixel_radius_vec[pore_pixel_radius_vec .> 0.0]

    # create histogram of pore pixel radii
    pixel_radius_histogram = StatsBase.fit(
            StatsBase.Histogram, pore_pixel_radius_filtered_vec, 
            0.2501:0.5:minimum(structure_dict["size_data"]), 
            closed=:left)

    # normalize histogram
    pixel_radius_histogram = LinearAlgebra.normalize(pixel_radius_histogram, mode=:probability)

    # get pore size distribution
    pore_size_distribution = pixel_radius_histogram.weights

    # convert pixel radii to physical radii
    pore_size_vec = sphere_pixel_radius_vec .* structure_dict["voxel_edge_length"]

    

    # create dict to save
    pore_size_distribution_dict = Dict{String, Any}("pore_size_vec" => pore_size_vec,
                                "pore_size_distribution" => pore_size_distribution,
                                "nr_sampled_voxels" => nr_sampled_voxels)

    # add label to dictionary if label is not nothing
    if label !== nothing
        pore_size_distribution_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(pore_size_distribution_dict),
            save_path*"_pore_size_distribution.h5")

    end

    return pore_size_distribution_dict
end



"""
Save spatial network to a GML format file 
"""
function save_spatial_network_to_gml(spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    save_path::String 
        = raw"..\structures\random_networks\\")

    # open new file
    open(save_path*filename*".gml", "w") do opened_file

        # write header
        write(opened_file, "graph [ \n")

        # loop through vertices
        for vertex in MetaGraphsNext.labels(spatial_network)

            # write vertex
            write(opened_file, Format.format(
                "  node [\n    id {1}\n    position [ x {2} y {3} z {4} ]\n  ]\n",
                vertex,
                spatial_network[vertex]["position"][1],
                spatial_network[vertex]["position"][2],
                spatial_network[vertex]["position"][3]))

        end

        # loop through edges
        for edge in MetaGraphsNext.edge_labels(spatial_network)

            # write edge
            write(opened_file, Format.format(
                "  edge [\n    source {1}\n    target {2}\n    vector [ x {3} y {4} z {5} ]\n    distance_squared {6}\n  ]\n",
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
Get mesh from network
"""
function save_mesh_from_spatial_network(graph_dict::Dict, filename::String;
    bond_radius::Float64 = 0.3131,
    vector_out_of_supercell_length = 1/2,
    save_path::String = raw"..\structures\random_networks\\",
    duplicate_bonds_close_to_supercell_edge::Bool = true,
    format::String = "obj")

    # create graph dict to plot
    plot_dict = deepcopy(graph_dict)
    
    # cut all bonds that reach out of supercell and replace
    # them by new bonds to duplicated vertices outside of the supercell
    plot_dict = cut_bonds_out_of_supercell!(plot_dict; 
        vector_out_of_supercell_length = vector_out_of_supercell_length)

    # loop through bonds
    for bond in MetaGraphsNext.edge_labels(plot_dict["spatial_network"])

        # get bond's start and target positions and its direction vector
        start_pos = plot_dict["spatial_network"][bond[1]]["position"]
        target_pos = plot_dict["spatial_network"][bond[2]]["position"]
        # direction_vec = plot_dict["spatial_network"][bond...]["vector"]

        # create cylinder object
        current_cylinder = GeometryBasics.Cylinder(
            GeometryBasics.Point( start_pos...),
            GeometryBasics.Point( target_pos...),
            bond_radius)
        
        # mesh cylinder object
        current_cylinder_mesh = GeometryBasics.mesh(current_cylinder)

        # save mesh
        total_path = save_path*filename*"\\"*string(bond[1])*"_"*string(bond[2])*"."*format

        FileIO.save(total_path, current_cylinder_mesh)

        # if one of the two vertices is close to the supercell edge but the vertices
        # are not on opposite sides of the supercell, save another cylinder just outside the
        # supercell on the other side
        if (duplicate_bonds_close_to_supercell_edge
            && (any(start_pos .< bond_radius ) 
            || any(target_pos .< bond_radius ) 
            || any((graph_dict["supercell_edge_length"] .- start_pos) .< bond_radius )
            || any((graph_dict["supercell_edge_length"] .- target_pos) .< bond_radius ) )
            && LinearAlgebra.norm(start_pos .- target_pos) < graph_dict["supercell_edge_length"]/2
            && all(start_pos .< graph_dict["supercell_edge_length"] )
            && all(target_pos .< graph_dict["supercell_edge_length"] )
            && all(start_pos .> 0.0 )
            && all(target_pos .> 0.0 )
            )

            # check on which side of supercell the additional bond should be added and
            # calculate new start and target positions
            if any(start_pos .< bond_radius ) || any(target_pos .< bond_radius ) 
                
                new_start_pos = (start_pos .+ graph_dict["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) .|| (target_pos .< bond_radius ) ))

                new_target_pos = (target_pos .+ graph_dict["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) .|| (target_pos .< bond_radius ) ))
            else
                new_start_pos = (start_pos .- graph_dict["supercell_edge_length"]
                    .* (((graph_dict["supercell_edge_length"] .- start_pos) .< bond_radius )
                    .|| ((graph_dict["supercell_edge_length"] .- target_pos) .< bond_radius )))

                new_target_pos = (target_pos .- graph_dict["supercell_edge_length"]
                    .* (((graph_dict["supercell_edge_length"] .- start_pos) .< bond_radius )
                    .|| ((graph_dict["supercell_edge_length"] .- target_pos) .< bond_radius )))
            end

            #println(start_pos, target_pos, new_start_pos, new_target_pos)
                
            # create cylinder object
            current_cylinder = GeometryBasics.Cylinder(
                GeometryBasics.Point( new_start_pos...),
                GeometryBasics.Point( new_target_pos...),
                bond_radius)
            
            # mesh cylinder object
            current_cylinder_mesh = GeometryBasics.mesh(current_cylinder)

            # save mesh
            total_path = save_path*filename*"\\"*string(bond[1])*"_"*string(bond[2])*"_outside_supercell."*format

            FileIO.save(total_path, current_cylinder_mesh)
        end

    end

    return
end


"""
Load spatial network and its properties from a MGformat file and a h5 dictionary
"""
function load_spatial_network_from_h5_and_MGformat(dict_path_without_format::String)

    # load spatial network in MGformat
    spatial_network = MetaGraphsNext.loadgraph(
            dict_path_without_format*".mg", MetaGraphsNext.MGFormat())

    # load rest of spatial network dict
    spatial_network = GU.load_h5_dict(dict_path_without_format*".h5")

    # add spatial network key to graph dict
    spatial_network = spatial_network

    return spatial_network
end




"""
Save spatial network to a GML format file and the rest of spatial_network and
evolution_dict to an h5 file
"""
function save_graph_to_h5_and_gml(spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    evolution_dict = nothing,
    save_path::String 
        = raw"..\structures\random_networks\\")

    # save evolution dict if passed
    if evolution_dict !== nothing
        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
    end

    # create copy of spatial_network to not change the original file
    spatial_network_to_save = deepcopy(spatial_network)

    # save graph to gml format
    save_spatial_network_to_gml(spatial_network, filename; save_path=save_path)

    # remove spatial_network from spatial_network
    delete!(spatial_network_to_save, "spatial_network")

    # save graph dict
    GU.save_dict_to_h5(spatial_network_to_save, save_path*filename*".h5")

    return
end


"""
Load graph and its properties from a GML file and a h5 dictionary
"""
function load_graph_from_h5_and_gml(dict_path_without_format::String)

    # load spatial network in MGformat
    spatial_network = load_spatial_network_from_gml(
            dict_path_without_format*".gml")

    # load rest of graph dict
    spatial_network = GU.load_h5_dict(dict_path_without_format*".h5")

    # add spatial network key to graph dict
    spatial_network = spatial_network

    return spatial_network
end




"""
Save spatial network to an MGformat file and the rest of spatial_network and
evolution_dict to an h5 file
"""
function save_spatial_network_to_h5_and_MGformat(spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    evolution_dict = nothing,
    save_path::String 
        = raw"..\structures\random_networks\\")

    # save evolution dict if passed
    if evolution_dict !== nothing
        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
    end

    # create copy of spatial_network to not change the original file
    spatial_network_to_save = deepcopy(spatial_network)

    # save spatial network to MGformat
    MetaGraphsNext.savegraph(save_path*filename*".mg", spatial_network_to_save["spatial_network"])

    # remove spatial_network from spatial_network
    delete!(spatial_network_to_save, "spatial_network")

    # save spatial network dict
    GU.save_dict_to_h5(spatial_network_to_save, save_path*filename*".h5")

    return
end


"""
Convert a graph in MGformat to a gml format
"""
function convert_MGformat_to_gml(
    filename::String;
    save_path::String 
        = raw"..\structures\random_networks\\")

    # load graph dict
    spatial_network = load_graph_from_h5_and_MGformat(save_path*filename)

    # save graph to gml format
    save_spatial_network_to_gml(spatial_network, filename; save_path=save_path)

    return
end




"""
Convert all files in a directory in MGformat to gml format
"""
function convert_all_files_in_directory_MGformat_to_gml(directory_path::String)

    # get all files in directory
    filenames = readdir(directory_path)

    # loop through files
    for filename in filenames

        if endswith(filename, ".mg")

            # convert file to gml format
            convert_MGformat_to_gml(filename[1:end-3]; save_path=directory_path)
        end

    end

    return
    
end



"""
For each bond in the network, if both its vertices are on the same side of
the supercell but at least one of them lies close to the supercell edge,
duplicate the bond on the other side of the supercell just outside the edge.
This is required when cylinders are assigned to the bonds and it is plotted or
used in an optical simulation.
"""
function duplicate_bonds_close_to_supercell_edge!(
    spatial_network::MetaGraphsNext.MetaGraph;
    bond_radius::Float64 = 0.35)

    # count current vertex
    vertex_count = copy(spatial_network[]["nr_vertices"])

    # loop through bonds
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # get bond's start and target positions and its direction vector
        start_pos = spatial_network[bond[1]]["position"]
        target_pos = spatial_network[bond[2]]["position"]

        # if one of the two vertices is close to the supercell edge but the
        # vertices are not on opposite sides of the supercell, save another
        # cylinder just outside the supercell on the other side
        if ((any(start_pos .< bond_radius ) 
            || any(target_pos .< bond_radius ) 
            || any((spatial_network[]["supercell_edge_length"] .- start_pos) 
                .< bond_radius )
            || any((spatial_network[]["supercell_edge_length"] .- target_pos) 
                .< bond_radius ) )
            && LinearAlgebra.norm(start_pos .- target_pos) 
                < spatial_network[]["supercell_edge_length"]/2
            && all(start_pos .< spatial_network[]["supercell_edge_length"] )
            && all(target_pos .< spatial_network[]["supercell_edge_length"] )
            && all(start_pos .> 0.0 )
            && all(target_pos .> 0.0 ) )

            # check on which side of supercell the additional bond should be
            # added and calculate new start and target positions
            if (any(start_pos .< bond_radius ) 
                || any(target_pos .< bond_radius ) )
                
                new_start_pos = (
                    start_pos .+ spatial_network[]["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) 
                    .|| (target_pos .< bond_radius ) ))

                new_target_pos = (target_pos 
                    .+ spatial_network[]["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) 
                    .|| (target_pos .< bond_radius ) ))
            else
                new_start_pos = (start_pos 
                    .- spatial_network[]["supercell_edge_length"]
                    .* (((spatial_network[]["supercell_edge_length"] 
                        .- start_pos) .< bond_radius )
                    .|| ((spatial_network[]["supercell_edge_length"] 
                        .- target_pos) .< bond_radius )))

                new_target_pos = (target_pos 
                    .- spatial_network[]["supercell_edge_length"]
                    .* (((spatial_network[]["supercell_edge_length"] 
                        .- start_pos) .< bond_radius )
                    .|| ((spatial_network[]["supercell_edge_length"] 
                        .- target_pos) .< bond_radius )))
            end

            # add two new vertices and the bond between them to the spatial
            # network
            spatial_network[vertex_count + 1] = (
                    Dict("position" => new_start_pos) )
            spatial_network[vertex_count + 2] = (
                        Dict("position" => new_target_pos) )

            spatial_network[vertex_count + 1, vertex_count + 2] = (
                Dict("vector" => (new_target_pos .- new_start_pos), 
                    "distance_squared" => (
                LinearAlgebra.norm(new_target_pos .- new_start_pos)^2 )) )

            vertex_count += 2
        end
    end

    spatial_network[]["nr_vertices"] = vertex_count

    return spatial_network
end



"""
Calculate angle averaged spectral density from 3d array of spectral density
"""
function get_spectral_density_angle_averaged(
    spectral_density_dict::Dict;
    gaussian_filter::Bool = true,
    gaussian_filter_sigma_x::Float64 = 2*pi/25, 
    gaussian_filter_filtered_data_x_step_length::Float64 = 2*pi/25,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name")

    # create a dictionary for angle averaged spectral density
    # with the length squared of the index vector as the key,
    # because it is proportional to the square of the wavenumber
    spectral_density_angle_averaged_dict = Dict{Int64, Vector{Float64}}()

    # determine origin vector of cartesian coordinates
    index_vector_origin = Int.(floor.( (size(
        spectral_density_dict["spectral_density_array"])[1:3] .+ 1 ) ./ 2) )

    for i in CartesianIndices(spectral_density_dict["spectral_density_array"])

        # determine actual index vector
        index_vector = i.I .- index_vector_origin

        # calculate length squared of index vector
        index_vector_length_squared = sum(index_vector.^2)

        # check if key exists in dictionary
        if index_vector_length_squared in keys(
            spectral_density_angle_averaged_dict)

            # add spectral density to dictionary
            push!(spectral_density_angle_averaged_dict[index_vector_length_squared], 
                abs(spectral_density_dict["spectral_density_array"][i]))

        else

            # add key and spectral density to dictionary
            spectral_density_angle_averaged_dict[index_vector_length_squared] = 
                [abs(spectral_density_dict["spectral_density_array"][i])]

        end
    end
    
    # initialize wavenumber vector and spectral density vector
    unfiltered_wavenumber_vec = Vector{Float64}()
    unfiltered_spectral_density_vec = Vector{
        Measurements.Measurement{Float64}}()

    # get lattice constant of reciprocal lattice
    reciprocal_lattice_constant = LinearAlgebra.norm(
        spectral_density_dict["wavevector_array"][(
            index_vector_origin .+ [1,0,0])...,:])

    # get vector of wavenumbers and angle averaged spectral density including
    # their uncertainty
    for key in keys(spectral_density_angle_averaged_dict)

        # get wavenumber from wavevector
        push!(unfiltered_wavenumber_vec, reciprocal_lattice_constant*sqrt(key))

        # get angle averaged spectral density
        push!(unfiltered_spectral_density_vec, 
            Measurements.measurement(Statistics.mean(
                spectral_density_angle_averaged_dict[key]),
            Statistics.std(spectral_density_angle_averaged_dict[key])))

    end

    # sort wavenumber vector and spectral density vector
    unfiltered_spectral_density_vec = unfiltered_spectral_density_vec[
        sortperm(unfiltered_wavenumber_vec)]
    sort!(unfiltered_wavenumber_vec)

    # create dict to save
    spectral_density_angle_averaged_dict = Dict{String, Any}(
        "unfiltered_wavenumber_vec" => unfiltered_wavenumber_vec, 
        "unfiltered_spectral_density_vec" => unfiltered_spectral_density_vec,
        "label" => spectral_density_dict["label"])

    # apply gaussian filter if desired
    if gaussian_filter
        filtered_data_x, filtered_data_y = GU.gaussian_filter_1d(
            unfiltered_wavenumber_vec[2:end], 
            unfiltered_spectral_density_vec[2:end]; 
            sigma_x=gaussian_filter_sigma_x, 
            filtered_data_x_step_length
                =gaussian_filter_filtered_data_x_step_length)

        spectral_density_angle_averaged_dict["wavenumber_vec"] = (
            filtered_data_x)
        spectral_density_angle_averaged_dict["spectral_density_vec"] = (
            filtered_data_y)
        spectral_density_angle_averaged_dict["gaussian_filter_sigma_x"] = (
            gaussian_filter_sigma_x)
    end

    if save_result
        GU.save_dict_to_h5(copy(spectral_density_angle_averaged_dict),
            save_path*"_spectral_density_angle_averaged.h5")
    end
                                            
    return spectral_density_angle_averaged_dict
end


"""
Define a anisotropy metric as the ratio of uncertainty and value of the sum of
the spectral density for several small wavenumbers normalized by the same value
for the diamond lattice
"""
function get_anisotropy_metric_from_spectral_density(
    spectral_density_angle_averaged_dict::Dict;
    diamond_std_value_ratio = 1.2588849)

    # set the wavenumbers where spectral density will be checked
    wavenumbers_to_check_vec = 2*pi*collect(0.2:0.1:1.0)

    # get the wavenumbers that lie the clostest to the wavenumbers
    # to be checked
    index_vec = [argmin(abs.(
        spectral_density_angle_averaged_dict["wavenumber_vec"] 
            .- wavenumbers_to_check_vec[i]))
        for i in eachindex(wavenumbers_to_check_vec)]

    # sum the checked spectral densities
    summed_spectral_density_to_check = sum( 
        spectral_density_angle_averaged_dict["spectral_density_vec"][
            index_vec] )

    # determine ratio of uncertainty to value
    std_value_ratio = (Measurements.uncertainty(
            summed_spectral_density_to_check) 
        / Measurements.value(summed_spectral_density_to_check))

    # normalize this ratio by the corresponding ratio for the diamond
    anisotropy_metric_from_spectral_density = (std_value_ratio 
    / diamond_std_value_ratio)

    return anisotropy_metric_from_spectral_density
end


function get_cluster_metric(correlation_functions_dict::Dict;
    nr_vertices_second_neighbor_shell::Int = 17)

    # get distance where cumulative coordination nr equals the given number of 
    # vertices minus one to account for the fact that the central vertex is not
    # contained in the cumulative coordination number 
    cluster_metric = correlation_functions_dict["vertex_distance_vec"][
        findfirst(x -> x > nr_vertices_second_neighbor_shell, 
        correlation_functions_dict["cumulative_coord_nr_vec"])
        ]/second_neighbor_distance_periodic_network

    return cluster_metric
end


"""
Calculate the statistical difference between two networks by comparing their
bond length and bond angle standard deviations
"""
function get_statistical_difference(
    bond_length_std_1,
    bond_length_std_2,
    bond_angle_std_1,
    bond_angle_std_2)

    # calculate ratio in bond length standard deviation
    bond_length_std_ratio = (bond_length_std_1 / bond_length_std_2)
    
    # calculate ratio in bond angle standard deviation
    bond_angle_std_ratio = (bond_angle_std_1 / bond_angle_std_2)

    # define difference metric delta
    delta = sqrt( abs2(bond_length_std_ratio - 1) 
        + abs2(bond_angle_std_ratio - 1))

    return delta
end


"""
Create a padding of 0.6 bond length outside of the supercell for all bonds that
are close to the supercell edge. I choose 0.6 times the equilibrium bond length
to account for bonds that are slightly longer than the equilibrium bond length
"""
function duplicate_bonds_close_to_supercell_edge!(
    spatial_network::MetaGraphsNext.MetaGraph)

    # count current vertex
    vertex_count = copy(spatial_network[]["nr_vertices"])

    original_bond_vec = collect(MetaGraphsNext.edge_labels(spatial_network))

    # loop through bonds
    for bond in original_bond_vec

        # get bond's start and target positions and its direction vector
        start_pos = spatial_network[bond[1]]["position"]
        target_pos = spatial_network[bond[2]]["position"]

        identity_matrix = Matrix(LinearAlgebra.I, 3, 3)

        # loop through all dimensions
        for i in 1:3

            # copy lower face
            if (((start_pos[i] < 0.6) && (start_pos[i] > 0)) 
                || ((target_pos[i] < 0.6) && (target_pos[i] > 0)))

                new_start_pos = (start_pos 
                    .+ spatial_network[]["supercell_edge_length"] 
                    .* identity_matrix[i,:])

                new_target_pos = (target_pos 
                    .+ spatial_network[]["supercell_edge_length"] 
                    .* identity_matrix[i,:])

                # add two new vertices and the bond between them to the spatial
                # network
                spatial_network[vertex_count + 1] = (
                    Dict("position" => new_start_pos) )
                spatial_network[vertex_count + 2] = (
                    Dict("position" => new_target_pos) )

                spatial_network[vertex_count + 1, vertex_count + 2] = (
                    Dict("vector" => spatial_network[bond...]["vector"], 
                        "distance_squared" 
                            => spatial_network[bond...]["distance_squared"]))

                vertex_count += 2

                # copy edges
                for j in i+1:3

                    # copy lower edge
                    if (((start_pos[j] < 0.6) && (start_pos[j] > 0))
                        || ((target_pos[j] < 0.6) && (target_pos[j] > 0)))

                        new_start_pos = (start_pos 
                            .+ spatial_network[]["supercell_edge_length"] 
                            .* (identity_matrix[i,:] .+ identity_matrix[j,:]))

                        new_target_pos = (target_pos 
                            .+ spatial_network[]["supercell_edge_length"] 
                            .* (identity_matrix[i,:] .+ identity_matrix[j,:]))

                        # add two new vertices and the bond between them to the
                        # spatial network
                        spatial_network[vertex_count + 1] = (
                            Dict("position" => new_start_pos) )
                        spatial_network[vertex_count + 2] = (
                            Dict("position" => new_target_pos) )

                        spatial_network[vertex_count + 1, vertex_count + 2] = (
                            Dict("vector" 
                                    => spatial_network[bond...]["vector"], 
                                "distance_squared" 
                                    => spatial_network[bond...][
                                        "distance_squared"]))

                        vertex_count += 2

                        # copy corners
                        if j == 2
                            k=3

                            # copy lower corner 
                            if (((start_pos[k] < 0.6) && (start_pos[k] > 0))
                                || ((target_pos[k] < 0.6) 
                                    && (target_pos[k] > 0)))
                                new_start_pos = (start_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (identity_matrix[i,:] 
                                        .+ identity_matrix[j,:] 
                                        .+ identity_matrix[k,:]))

                                new_target_pos = (target_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (identity_matrix[i,:] 
                                        .+ identity_matrix[j,:] 
                                        .+ identity_matrix[k,:]))

                                # add two new vertices and the bond between 
                                # them to the spatial network
                                spatial_network[vertex_count + 1] = (
                                    Dict("position" => new_start_pos) )
                                spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )

                                spatial_network[vertex_count + 1, 
                                    vertex_count + 2] = (
                                    Dict("vector" 
                                            => spatial_network[bond...][
                                                "vector"], 
                                        "distance_squared" 
                                            => spatial_network[bond...][
                                                "distance_squared"]))
                                vertex_count += 2

                            # copy upper corner
                            elseif (((start_pos[k] 
                                    > spatial_network[][
                                        "supercell_edge_length"] 
                                        - 0.6) 
                                && (start_pos[k] 
                                    < spatial_network[][
                                        "supercell_edge_length"])) 
                                || ((target_pos[k] 
                                    > spatial_network[][
                                        "supercell_edge_length"] 
                                        - 0.6) 
                                && (target_pos[k] 
                                    < spatial_network[][
                                        "supercell_edge_length"])))

                                new_start_pos = (start_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (identity_matrix[i,:] 
                                        .+ identity_matrix[j,:] 
                                        .- identity_matrix[k,:]))
                                new_target_pos = (target_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (identity_matrix[i,:] 
                                        .+ identity_matrix[j,:] 
                                        .- identity_matrix[k,:]))

                                # add two new vertices and the bond between
                                # them to the spatial network
                                spatial_network[vertex_count + 1] = (
                                    Dict("position" => new_start_pos) )
                                spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )
                                spatial_network[vertex_count + 1, 
                                    vertex_count + 2] = (
                                    Dict("vector" 
                                            => spatial_network[bond...][
                                                "vector"], 
                                        "distance_squared" 
                                            => spatial_network[bond...][
                                                "distance_squared"]))
                                vertex_count += 2
                            end
                        end

                    # copy upper edge
                    elseif (((start_pos[j] 
                            > spatial_network[]["supercell_edge_length"] - 0.6) 
                        && (start_pos[j] 
                            < spatial_network[]["supercell_edge_length"])) 
                        || ((target_pos[j] 
                            > spatial_network[]["supercell_edge_length"] - 0.6) 
                        && (target_pos[j] 
                            < spatial_network[]["supercell_edge_length"])))

                        new_start_pos = (start_pos 
                            .+ spatial_network[]["supercell_edge_length"] 
                            .* (identity_matrix[i,:] .- identity_matrix[j,:]))

                        new_target_pos = (target_pos 
                            .+ spatial_network[]["supercell_edge_length"] 
                            .* (identity_matrix[i,:] .- identity_matrix[j,:]))

                        # add two new vertices and the bond between them to the
                        # spatial network
                        spatial_network[vertex_count + 1] = (
                                Dict("position" => new_start_pos) )
                        spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )

                        spatial_network[vertex_count + 1, vertex_count + 2] = (
                            Dict("vector" 
                                    => spatial_network[bond...]["vector"], 
                                "distance_squared" 
                                    => spatial_network[bond...][
                                        "distance_squared"]))

                        vertex_count += 2

                        # copy corners
                        if j == 2
                            k=3

                            # copy lower corner 
                            if (((start_pos[k] < 0.6) && (start_pos[k] > 0)) 
                                || ((target_pos[k] < 0.6) 
                                    && (target_pos[k] > 0)))

                                new_start_pos = (start_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (identity_matrix[i,:] 
                                        .- identity_matrix[j,:] 
                                        .+ identity_matrix[k,:]))

                                new_target_pos = (target_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (identity_matrix[i,:] 
                                        .- identity_matrix[j,:] 
                                        .+ identity_matrix[k,:]))

                                # add two new vertices and the bond between
                                # them to the spatial network
                                spatial_network[vertex_count + 1] = (
                                    Dict("position" => new_start_pos) )
                                spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )

                                spatial_network[vertex_count + 1, 
                                    vertex_count + 2] = (
                                    Dict("vector" 
                                            => spatial_network[bond...][
                                                "vector"], 
                                        "distance_squared" 
                                            => spatial_network[bond...][
                                                "distance_squared"]))
                                vertex_count += 2

                            # copy upper corner
                            elseif (((start_pos[k] 
                                    > spatial_network[][
                                        "supercell_edge_length"] 
                                        - 0.6) 
                                && (start_pos[k] 
                                    < spatial_network[][
                                        "supercell_edge_length"])) 
                                || ((target_pos[k] 
                                    > spatial_network[][
                                        "supercell_edge_length"] 
                                        - 0.6) 
                                && (target_pos[k] 
                                    < spatial_network[][
                                        "supercell_edge_length"])))

                                new_start_pos = (start_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (identity_matrix[i,:] 
                                        .- identity_matrix[j,:] 
                                        .- identity_matrix[k,:]))
                                new_target_pos = (target_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (identity_matrix[i,:] 
                                        .- identity_matrix[j,:] 
                                        .- identity_matrix[k,:]))

                                # add two new vertices and the bond between
                                # them to the spatial network
                                spatial_network[vertex_count + 1] = (
                                    Dict("position" => new_start_pos) )
                                spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )
                                spatial_network[vertex_count + 1, 
                                    vertex_count + 2] = (
                                    Dict("vector" 
                                            => spatial_network[bond...][
                                                "vector"], 
                                        "distance_squared" 
                                            => spatial_network[bond...][
                                                "distance_squared"]))
                                vertex_count += 2
                            end
                        end
                    end
                end

            # copy upper face
            elseif (((start_pos[i] 
                        > spatial_network[]["supercell_edge_length"] - 0.6) 
                    && (start_pos[i] 
                        < spatial_network[]["supercell_edge_length"]))
                || ((target_pos[i] 
                        > spatial_network[]["supercell_edge_length"] - 0.6) 
                    && (target_pos[i] 
                        < spatial_network[]["supercell_edge_length"])))

                new_start_pos = (start_pos 
                    .- spatial_network[]["supercell_edge_length"] 
                    .* identity_matrix[i,:])

                new_target_pos = (target_pos 
                    .- spatial_network[]["supercell_edge_length"] 
                    .* identity_matrix[i,:])

                # add two new vertices and the bond between them to the spatial
                # network
                spatial_network[vertex_count + 1] = (
                        Dict("position" => new_start_pos) )
                spatial_network[vertex_count + 2] = (
                            Dict("position" => new_target_pos) )

                spatial_network[vertex_count + 1, vertex_count + 2] = (
                    Dict("vector" => spatial_network[bond...]["vector"], 
                        "distance_squared" 
                            => spatial_network[bond...]["distance_squared"]))

                vertex_count += 2

                # copy edges
                for j in i+1:3

                    # copy lower edge
                    if (((start_pos[j] < 0.6) && (start_pos[j] > 0)) 
                        || ((target_pos[j] < 0.6) && (target_pos[j] > 0)))

                        new_start_pos = (start_pos 
                            .+ spatial_network[]["supercell_edge_length"] 
                            .* (.- identity_matrix[i,:] 
                                .+ identity_matrix[j,:]))

                        new_target_pos = (target_pos 
                            .+ spatial_network[]["supercell_edge_length"] 
                            .* (.- identity_matrix[i,:] 
                                .+ identity_matrix[j,:]))

                        # add two new vertices and the bond between them to the
                        # spatial network
                        spatial_network[vertex_count + 1] = (
                                Dict("position" => new_start_pos) )
                        spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )

                        spatial_network[vertex_count + 1, vertex_count + 2] = (
                            Dict("vector" 
                                    => spatial_network[bond...]["vector"], 
                                "distance_squared" 
                                    => spatial_network[bond...][
                                        "distance_squared"]))

                        vertex_count += 2

                        # copy corners
                        if j == 2
                            k=3

                            # copy lower corner 
                            if (((start_pos[k] < 0.6) && (start_pos[k] > 0)) 
                                || ((target_pos[k] < 0.6) 
                                    && (target_pos[k] > 0)))

                                new_start_pos = (start_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (.- identity_matrix[i,:] 
                                        .+ identity_matrix[j,:] 
                                        .+ identity_matrix[k,:]))

                                new_target_pos = (target_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (.- identity_matrix[i,:] 
                                        .+ identity_matrix[j,:] 
                                        .+ identity_matrix[k,:]))

                                # add two new vertices and the bond between
                                # them to the spatial network
                                spatial_network[vertex_count + 1] = (
                                    Dict("position" => new_start_pos) )
                                spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )

                                spatial_network[vertex_count + 1, 
                                    vertex_count + 2] = (
                                    Dict("vector" 
                                            => spatial_network[bond...][
                                                "vector"], 
                                        "distance_squared" 
                                            => spatial_network[bond...][
                                                "distance_squared"]))
                                vertex_count += 2
                                

                            # copy upper corner
                            elseif (((start_pos[k] 
                                    > spatial_network[][
                                        "supercell_edge_length"] 
                                        - 0.6) 
                                && (start_pos[k] 
                                    < spatial_network[][
                                        "supercell_edge_length"])) 
                                || ((target_pos[k] 
                                    > spatial_network[][
                                        "supercell_edge_length"] 
                                        - 0.6) 
                                && (target_pos[k] 
                                    < spatial_network[][
                                        "supercell_edge_length"])))

                                new_start_pos = (start_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (.- identity_matrix[i,:] 
                                        .+ identity_matrix[j,:] 
                                        .- identity_matrix[k,:]))
                                new_target_pos = (target_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (.- identity_matrix[i,:] 
                                        .+ identity_matrix[j,:] 
                                        .- identity_matrix[k,:]))

                                # add two new vertices and the bond between
                                # them to the spatial network
                                spatial_network[vertex_count + 1] = (
                                    Dict("position" => new_start_pos) )
                                spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )
                                spatial_network[vertex_count + 1, 
                                    vertex_count + 2] = (
                                    Dict("vector" 
                                            => spatial_network[bond...][
                                                "vector"], 
                                        "distance_squared" 
                                            => spatial_network[bond...][
                                                "distance_squared"]))
                                vertex_count += 2
                                
                            end
                        end

                    # copy upper edge
                    elseif (((start_pos[j] 
                            > spatial_network[]["supercell_edge_length"] - 0.6) 
                        && (start_pos[j] 
                            < spatial_network[]["supercell_edge_length"])) 
                        || ((target_pos[j] 
                            > spatial_network[]["supercell_edge_length"] - 0.6) 
                        && (target_pos[j] 
                            < spatial_network[]["supercell_edge_length"])))

                        new_start_pos = (start_pos 
                            .- spatial_network[]["supercell_edge_length"] 
                            .* (identity_matrix[i,:] .+ identity_matrix[j,:]))

                        new_target_pos = (target_pos 
                            .- spatial_network[]["supercell_edge_length"] 
                            .* (identity_matrix[i,:] .+ identity_matrix[j,:]))

                        # add two new vertices and the bond between them to the
                        # spatial network
                        spatial_network[vertex_count + 1] = (
                                Dict("position" => new_start_pos) )
                        spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos))

                        spatial_network[vertex_count + 1, vertex_count + 2] = (
                            Dict("vector" 
                                    => spatial_network[bond...]["vector"], 
                                "distance_squared" 
                                    => spatial_network[bond...][
                                        "distance_squared"]))

                        vertex_count += 2

                        # copy corners
                        if j == 2
                            k=3

                            # copy lower corner 
                            if (((start_pos[k] < 0.6) && (start_pos[k] > 0)) 
                                || ((target_pos[k] < 0.6) 
                                    && (target_pos[k] > 0)))

                                new_start_pos = (start_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (.- identity_matrix[i,:] 
                                        .- identity_matrix[j,:] 
                                        .+ identity_matrix[k,:]))

                                new_target_pos = (target_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (.- identity_matrix[i,:] 
                                        .- identity_matrix[j,:] 
                                        .+ identity_matrix[k,:]))

                                # add two new vertices and the bond between
                                # them to the spatial network
                                spatial_network[vertex_count + 1] = (
                                    Dict("position" => new_start_pos) )
                                spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )

                                spatial_network[vertex_count + 1, 
                                    vertex_count + 2] = (
                                    Dict("vector" 
                                            => spatial_network[bond...][
                                                "vector"], 
                                        "distance_squared" 
                                            => spatial_network[bond...][
                                                "distance_squared"]))
                                vertex_count += 2
                                

                            # copy upper corner
                            elseif (((start_pos[k] 
                                    > spatial_network[][
                                        "supercell_edge_length"] 
                                        - 0.6) 
                                && (start_pos[k] 
                                    < spatial_network[][
                                        "supercell_edge_length"])) 
                                || ((target_pos[k] 
                                    > spatial_network[][
                                        "supercell_edge_length"] 
                                        - 0.6) 
                                && (target_pos[k] 
                                    < spatial_network[][
                                        "supercell_edge_length"])))

                                new_start_pos = (start_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (.- identity_matrix[i,:] 
                                        .- identity_matrix[j,:] 
                                        .- identity_matrix[k,:]))
                                new_target_pos = (target_pos 
                                    .+ spatial_network[][
                                        "supercell_edge_length"] 
                                    .* (.- identity_matrix[i,:] 
                                        .- identity_matrix[j,:] 
                                        .- identity_matrix[k,:]))

                                # add two new vertices and the bond between
                                # them to the spatial network
                                spatial_network[vertex_count + 1] = (
                                    Dict("position" => new_start_pos) )
                                spatial_network[vertex_count + 2] = (
                                    Dict("position" => new_target_pos) )
                                spatial_network[vertex_count + 1, 
                                    vertex_count + 2] = (
                                    Dict("vector" 
                                            => spatial_network[bond...][
                                                "vector"], 
                                        "distance_squared" 
                                            => spatial_network[bond...][
                                                "distance_squared"]))

                                vertex_count += 2
                            end
                        end
                    end 
                end
            end 
        end
    end

    spatial_network[]["nr_vertices"] = vertex_count

    return spatial_network
end


"""
Measure the standard deviation of dihedral angles. This function might need to
be revisited when considering other coordination numbers than 4 because the
dihedral angle might have a peak around 0 which is not considered here
"""
function get_dihedral_angle_std(spatial_network::MetaGraphsNext.MetaGraph)

    # initialize vector of diehedral angles
    dihedral_angle_vec = Vector{Float64}()

    # loop through all bonds
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # get vector along bond
        bond_vec = spatial_network[bond...]["vector"]

        # loop through all neighbors of one vertex 
        for first_neighbor in setdiff( MetaGraphsNext.neighbor_labels(
            spatial_network, bond[1]), bond[2])

            # get vector from first neighbor to first bond vertex
            first_neighbor_to_bond_vertex_vec = ( sign(bond[1] 
                - first_neighbor)*spatial_network[first_neighbor, 
                    bond[1]]["vector"])

            # loop through all neighbors of other vertex
            for second_neighbor in setdiff( MetaGraphsNext.neighbor_labels(
                spatial_network, bond[2]), bond[1])

                # get vector from second bond vertex to second neighbor
                bond_vertex_to_second_neighbor_vec = ( sign(second_neighbor 
                    - bond[2])*spatial_network[bond[2], 
                        second_neighbor]["vector"])

                # calculate dihedral angle according to the equation given in
                # https://en.wikipedia.org/wiki/Dihedral_angle# 
                dihedral_angle = atan( LinearAlgebra.norm(bond_vec)
                        *LinearAlgebra.dot(first_neighbor_to_bond_vertex_vec,
                            LinearAlgebra.cross(bond_vec,
                                bond_vertex_to_second_neighbor_vec )),
                    LinearAlgebra.dot(
                        LinearAlgebra.cross(first_neighbor_to_bond_vertex_vec,
                            bond_vec ),
                        LinearAlgebra.cross(bond_vec,
                              bond_vertex_to_second_neighbor_vec )))

                # save dihedral angle
                push!(dihedral_angle_vec, dihedral_angle)
            end
        end
    end

    # save one peak of dihedral angle distribution
    lower_limit = 0
    upper_limit =  2 * pi / (spatial_network[]["coordination_nr"] - 1)

    dihedral_angle_one_peak_vec = dihedral_angle_vec[
        (dihedral_angle_vec .> lower_limit) .& (dihedral_angle_vec .< upper_limit)]

    # determine standard deviation
    dihedral_angle_std = Statistics.std(dihedral_angle_one_peak_vec)

    return [dihedral_angle_std, dihedral_angle_vec]
end


"""
Get hyperuniformity metric which is the structure factor at zero momentum
normalized by the height of the first peak in the structure factor as defined
in equation 251 in 10.1016/j.physrep.2018.03.001
"""
function get_hyperuniformity_metric(structure_factor_dict::Dict)

    # locate peaks of structure factor
    pks, vals = Peaks.findmaxima(structure_factor_dict["structure_factor_vec"])

    # cut structure factor data at momentum just above first peak
    structure_factor_cut_vec = structure_factor_dict[
        "structure_factor_vec"][1:pks[2]-1]
    wavenumber_cut_vec = structure_factor_dict["wavenumber_vec"][1:pks[2]-1]

    # set the order of the fitted polynomial
    polynomial_order = 5

    # fit polynomial of given order to cut data
    polynomial_fit = Polynomials.fit(wavenumber_cut_vec, 
                                    structure_factor_cut_vec,
                                    polynomial_order)

    # get extrapolated structure factor at zero momentum
    structure_factor_zero_momentum = polynomial_fit(0)

    # get first derivative of polynomial
    polynomial_derivative = Polynomials.derivative(polynomial_fit)

    # in case the structure factor is provided with uncertainty, obtain the
    # values
    polynomial_derivative_values = Polynomials.Polynomial(Measurements.value.( 
        collect(polynomial_derivative) ))
    
    # get critical momenta which is roots of first derivative of polynomial
    critical_momenta = Polynomials.roots(polynomial_derivative_values)

    # get real critical momenta
    critical_momenta_real = real.(
        critical_momenta[imag.(critical_momenta) .== 0])

    # get fitted structure factor at highest peak
    structure_factor_first_peak = maximum( 
        polynomial_fit.(critical_momenta_real) )

    # get hyperuniformity metric
    hyperuniformity_metric = (structure_factor_zero_momentum
        /structure_factor_first_peak)

    return [hyperuniformity_metric, polynomial_fit]
end



"""
Given a vector, this function returns all windows within the vector that range
over one decade, like from 1 to 10. The passed vector is the logarithm of the
vector containing the original values
"""
function get_one_decade_windows(
    log_vector::Vector{Float64})

    # get the minimum and maximum of the logarithm of the vector
    min_vector = minimum(log_vector)
    max_vector = maximum(log_vector)

    window_index_vec = Vector{Tuple{Int64, Int64}}()

    for i in eachindex(log_vector)
        if max_vector < log_vector[i]+1
            break
        end
        # get the element of the vector that is one decade larger than the
        # current element
        upper_index = argmin(abs.(log_vector .- (log_vector[i]+1)))

        push!(window_index_vec, (i, upper_index))
    end
    return window_index_vec
end


"""
Get the exponent alpha, that determines the scaling of the structure factor
or the spectral density with the wavenumber k for k -> 0 as in eq 25 of
10.1103/PhysRevE.109.064108. If this exponent is 0, the system is
nonhyperuniform, if it is >0, the system is hyperuniform. Three classes of
hyperuniformity are defined: 0 < alpha < 1, alpha = 1 and alpha > 1 as written
below eq 21 of 10.1103/PhysRevE.109.064108.
"""
function get_hyperuniformity_alpha(structure_factor_angle_averaged_dict::Dict;
    consider_spectral_density::Bool = false,
    t_range = (1e-5, 1))

    if consider_spectral_density
        data_type_string = "spectral_density"
    else
        data_type_string = "structure_factor"
    end

    # get the vector of times to sample the excess spreadability such that the
    # t values are equally spaced on a logarithmic scale
    time_vec = exp.(LinRange(log(t_range[1]), log(t_range[2]), 1000))

    # calculate the excess spreadability for each t value
    excess_spreadability_vec = [get_excess_spreadability(
        structure_factor_angle_averaged_dict, time_vec[i]; 
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

    # get logarithm of time and excess spreadability
    log_time_vec = log10.(time_vec)
    log_excess_spreadability_vec = log10.(excess_spreadability_vec)

    # get all windows within the time vector that range over one decade
    window_index_vec = get_one_decade_windows(log_time_vec)

    # perform a fit for each window and store the results of the fits
    fit_alpha_vec = Vector{Measurements.Measurement{Float64}}(undef, 
        length(window_index_vec))

    for i in eachindex(window_index_vec)
        # get the indices of the current window
        window_indices = window_index_vec[i]

        # fit a line to the data
        polynomial_fit = Polynomials.fit(
            log_time_vec[window_indices[1]:window_indices[2]], 
            log_excess_spreadability_vec[window_indices[1]:window_indices[2]], 
            1)

        # get alpha from the fit slope according to section IIIB of 
        # 10.1103/PhysRevE.109.064108
        alpha = -2 * polynomial_fit[1] - 3

        fit_alpha_vec[i] = alpha
        
        # print progress
        println("Progress: ", i/length(window_index_vec)*100, "%")
    end

    # determine the fit result with the smallest uncertainty
    min_uncertainty_index = argmin(Measurements.uncertainty.(fit_alpha_vec))
    
    # get the time window with the smallest uncertainty
    time_window = (log_time_vec[window_index_vec[min_uncertainty_index][1]], 
        log_time_vec[window_index_vec[min_uncertainty_index][2]])
    
    # get the exponent alpha and its uncertainty
    alpha = fit_alpha_vec[min_uncertainty_index]
    
    return [fit_alpha_vec, time_window, alpha]
end