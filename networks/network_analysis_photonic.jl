"""
These functions can be used to analyze the photonic response of networks.
"""


"""
Given an h5 dictionary with the keys 'transmission', 'reflection', and 'freqs',
this function returns a so called saturation_metric. For this metric, one
calculated two values. The first one is the maximal mean reflection inside a
window that has a size of 10% of its central frequency. The second value is the
smallest mean reflection inside a window that spans one octave of frequencies
and that fully contains the first window. The saturation metric is then defined 
as the ratio of these two values.
"""
function get_saturation_metric(r_t_dict::Dict{String, Any};
    consider_freq_range = [0.15, 0.35],
    freq_window_size = 0.1)

    # get the frequencies and the reflection
    freqs = r_t_dict["freqs"]
    reflection = r_t_dict["reflection"]

    # get the mean reflection in the frequency range
    mask = ((freqs .>= consider_freq_range[1]) 
        .& (freqs .<= consider_freq_range[2]))
    mean_reflection = Statistics.mean(reflection[mask])

    # search for the first frequeny f where 0.95*f is inside the frequency
    # range
    f1 = findfirst(x -> (1-freq_window_size/2)*x >= consider_freq_range[1], freqs)
    # search for the last frequeny f where 1.05*f is inside the frequency
    # range
    f2 = findlast(x -> (1+freq_window_size/2)*x <= consider_freq_range[2], freqs)
    # if no such frequency exists, return nothing
    if isnothing(f1) || isnothing(f2)
        @warn "No frequency found in the range $consider_freq_range. Returning nothing."
        return nothing
    end

    max_window_reflection = 0

    # loop through the frequencies between f1 and f2
    for i in f1:f2
        lower_index = findfirst(x -> x >= (1-freq_window_size/2)*freqs[i], freqs)
        upper_index = findlast(x -> x <= (1+freq_window_size/2)*freqs[i], freqs)

        # get the mean reflection in the window
        mean_window_reflection = Statistics.mean(
            reflection[lower_index:upper_index])

        if mean_window_reflection > max_window_reflection
            max_window_reflection = mean_window_reflection
        end
    end

    # calculate the saturation metric
    saturation_metric = max_window_reflection / mean_reflection

    return saturation_metric
end


"""
This function calculates the Pearson correlation and the Spearman correlation
between the saturation metric and the order metrics. It returns two
dictionaries with the keys being the order metric names and the values being
the correlation coefficients.
"""
function get_order_saturation_correlations(
    order_metric_dict::Dict{String, Any},    
    r_t_dict_path::String,
    analysis_data_path::String;
    consider_freq_range = [0.15, 0.35],
    freq_window_size = 0.1,
    save_results::Bool = false)
    
    # get all files in the directory that end with "r_t_only"
    files = readdir(r_t_dict_path, join=true)
    r_t_files = filter(x -> endswith(x, "r_t_only.hdf5"), files)

    # get the saturation metrics and the beta, t_max, and t_gradient values for
    # each file
    saturation_metric_vec = zeros(length(r_t_files))
    beta_t_max_t_gradient_vec = []
    
    for (i, file) in enumerate(r_t_files)
        # get the saturation metric
        r_t_dict = GU.load_h5_dict(file)
        saturation_metric = get_saturation_metric(r_t_dict, 
            consider_freq_range=consider_freq_range,
            freq_window_size=freq_window_size)
        saturation_metric_vec[i] = saturation_metric

        # get the beta, t_max, and t_gradient values from the filename
        filename = basename(file)
        parts = split(filename, '_')
        push!(beta_t_max_t_gradient_vec, (parse(Float64, parts[3]), 
        parse(Float64, parts[6]), parse(Float64, parts[9])))
    end

    # equally, get the combinations of beta, t_max, and t_gradient from the 
    # order metrics dict
    order_metric_beta_t_max_t_gradient_vec = [
        (order_metric_dict["bond_bending_const_vec"][i], 
        order_metric_dict["t_max_vec"][i], 
        order_metric_dict["t_gradient_vec"][i]) 
        for i in eachindex(order_metric_dict["bond_bending_const_vec"])]

    # find the indices where the entries of the
    # order_metric_beta_t_max_t_gradient_vec are the sames as the entries of 
    # the beta_t_max_t_gradient_vec
    for i in eachindex(beta_t_max_t_gradient_vec)
        if !(beta_t_max_t_gradient_vec[i] 
            in order_metric_beta_t_max_t_gradient_vec)
            @warn "The combination $(beta_t_max_t_gradient_vec[i]) 
            is not in the order metric dictionary."
            println(basename(r_t_files[i]))
        end
    end

    indices = [findfirst(x -> x == beta_t_max_t_gradient_vec[i],
        order_metric_beta_t_max_t_gradient_vec) 
        for i in eachindex(beta_t_max_t_gradient_vec)]

    # loop through keys of the dictionary which are all vectors and add them to
    # a list of vectors
    order_metric_vec_vec = []
    key_vec = []
    current_index = 1
    for key in keys(order_metric_dict)
        if (key == "filenames_vec" || key == "bond_bending_const_vec" 
            || key == "t_max_vec" || key == "t_gradient_vec" 
            || key == "total_keating_energy_vec")
            continue
        elseif key == "q_l_mat"
            for j in 1:13
                push!(order_metric_vec_vec, Measurements.value.(
                    order_metric_dict[key][j,indices]))
            end
            q_l_key_vec = ["q_$(i)_mean" for i in 0:12]
            append!(key_vec, q_l_key_vec)
            current_index += 13
            for j in 1:13
                push!(order_metric_vec_vec, Measurements.uncertainty.(
                    order_metric_dict[key][j,indices]))
            end
            q_l_key_vec = ["q_$(i)_std" for i in 0:12]
            append!(key_vec, q_l_key_vec)
            current_index += 13
        else
            push!(order_metric_vec_vec, order_metric_dict[key][indices])
            push!(key_vec, key)
            current_index += 1
        end
    end

    # using Pearson's correlation coefficient, calculate the correlation 
    # between the saturation metrics and the order metrics
    correlations = [Statistics.cor(saturation_metric_vec, 
       order_metric_vec_vec[i]) for i in 1:length(order_metric_vec_vec)]

    # sort keys according to the absolute values of their correlation 
    # coefficients
    sorted_indices = sortperm(abs.(correlations), rev=true)
    sorted_keys = key_vec[sorted_indices]
    sorted_correlations = correlations[sorted_indices]

    # Print the correlation coefficients
    for (i, key) in enumerate(sorted_keys)
        println("Pearson corr. for $key: ", sorted_correlations[i])
    end

    # create a dictionary with the keys being the order metric names and the
    # values being the correlation coefficients
    correlations_dict = Dict(zip(sorted_keys, sorted_correlations))

    if save_results
        GU.save_dict_to_h5(copy(correlations_dict),
            analysis_data_path*"order_saturation_pearson_correlations.h5")
    end

    spearman_correlations = [StatsBase.corspearman(saturation_metric_vec, 
        order_metric_vec_vec[i]) for i in 1:length(order_metric_vec_vec)]

    # sort keys according to the absolute values of their Spearman's rank
    # correlation coefficients
    sorted_spearman_indices = sortperm(abs.(spearman_correlations), rev=true)
    sorted_spearman_keys = key_vec[sorted_spearman_indices]
    sorted_spearman_correlations = (
        spearman_correlations[sorted_spearman_indices])

    # Print the Spearman's rank correlation coefficients
    for (i, key) in enumerate(sorted_spearman_keys)
        println("Spearman corr. for $key: ", sorted_spearman_correlations[i])
    end

    # create a dictionary with the keys being the order metric names and the
    # values being the Spearman's rank correlation coefficients
    spearman_correlations_dict = Dict(zip(sorted_spearman_keys, 
        sorted_spearman_correlations))

    if save_results
        GU.save_dict_to_h5(copy(spearman_correlations_dict),
            analysis_data_path*"order_saturation_spearman_correlations.h5")
    end

    return [correlations_dict, spearman_correlations_dict]
end



"""
Calculate the relative reflection peak to relative width ratio of a reflection 
spectrum. Here, the relative reflection peak is defined as the ratio of the 
peak height without background and the background at the frequency of the peak
center. The relative peak width is defined as the sigma of the peak over the
central peak frequency. The initial guesses for the fit parameters depend on
the network type. For a refractive index of n=1.5 and fill fractions of about
0.4, the initial guesses are as follows:
ctn: p0=[0.28, 0.7, 0.25, 0.2, 0.34, 0.05]
dia: p0=[0.28, 0.7, 0.25, 0.2, 0.4, 0.05]
lcs: p0=[0.28, 0.7, 0.25, 0.2, 0.38, 0.05]
srs: p0=[0.28, 0.7, 0.25, 0.2, 0.3, 0.05]
"""
function get_reflection_peak_height_to_width(r_t_dict::Dict{String, Any};
    p0::Vector{Float64}=[0.28, 0.7, 0.25, 0.2, 0.4, 0.05],
    upper_bounds::Vector{Float64}=[0.4, 0.8, 0.35, 0.4, 0.45, 0.1],
    lower_bounds::Vector{Float64}=[0.2, 0.55, 0.15, 0.0, 0.25, 0.01],
    save_plot::Bool = false,
    save_path::String = raw"..\..\photonics\tidy3d\plots\neural_network_networks\dia\run_1_2_r_t_low_n\some_filename")

    freqs = r_t_dict["freqs"]
    reflection = r_t_dict["reflection"]

    # define a fit model consisting of two gaussians, one broad one as the
    # background and a narrower one sitting on top of the first one
    background_peak(x, p) = (p[1] .* exp.((-1/2) .* ((x .- p[2]) ./ p[3]).^2))
    reflection_peak(x, p) = (p[4] .* exp.((-1/2) .* ((x .- p[5]) ./ p[6]).^2))

    model(x, p) = background_peak(x, p) .+ reflection_peak(x, p)

    # fit the model to the reflection data
    fit_result = LsqFit.curve_fit(model, freqs, reflection, p0, 
        lower=lower_bounds, upper=upper_bounds)
    
    covariance_matrix = LsqFit.estimate_covar(fit_result)
    standard_errors = sqrt.(LinearAlgebra.diag(covariance_matrix))

    fit_params = Measurements.measurement.(fit_result.param, standard_errors)

    # calculate the value of the background peak at the center frequency of the
    # reflection peak
    background_at_reflection_peak = background_peak(
        fit_params[5], fit_params[1:3])

    # calculate the relative peak height of the reflection peak
    relative_peak_height = fit_params[4] / background_at_reflection_peak

    # calculate the relative reflection peak width
    reflection_peak_width = fit_params[6] / fit_params[5]

    # calculate the ratio of relative peak height to relative peak width
    peak_height_to_width = relative_peak_height / reflection_peak_width

    if save_plot
        fitted_reflection = [model(freq, fit_params) for freq in freqs]

        Plots.plot(freqs, reflection, label="FDTD", xlabel="Frequency", 
            ylabel="Reflectance")
        Plots.plot!(freqs, Measurements.value.(fitted_reflection), 
            ribbon=Measurements.uncertainty.(fitted_reflection), label="Fit", 
            color=:orange)

        Plots.savefig(save_path*"_reflection_fit.png")
    end

    return peak_height_to_width
end


"""
This function calculates the Pearson correlation and the Spearman correlation
between the peak height to width ratio of a reflection spectrum and the order 
metrics. It returns two dictionaries with the keys being the order metric names 
and the values being the correlation coefficients. The initial guesses for the 
fit parameters of the reflection spectrum depend on the network type. For a 
refractive index of n=1.5 and fill fractions of about 0.4, the initial guesses 
are as follows:
ctn: p0=[0.28, 0.7, 0.25, 0.2, 0.34, 0.05]
dia: p0=[0.28, 0.7, 0.25, 0.2, 0.4, 0.05]
lcs: p0=[0.28, 0.7, 0.25, 0.2, 0.38, 0.05]
srs: p0=[0.28, 0.7, 0.25, 0.2, 0.3, 0.05]
"""
function get_order_peak_height_to_width_correlations(
    order_metric_dict::Dict{String, Any},    
    r_t_dict_path::String,
    analysis_data_path::String;
    p0::Vector{Float64}=[0.28, 0.7, 0.25, 0.2, 0.4, 0.05],
    upper_bounds::Vector{Float64}=[0.4, 0.8, 0.35, 0.4, 0.45, 0.1],
    lower_bounds::Vector{Float64}=[0.2, 0.55, 0.15, 0.0, 0.25, 0.01],
    save_results::Bool = false,
    print_progress::Bool = true)
    
    # get all files in the directory that end with "r_t_only"
    files = readdir(r_t_dict_path, join=true)
    r_t_files = filter(x -> endswith(x, "r_t_only.hdf5"), files)

    # get the peak to width metrics and the beta, t_max, and t_gradient values 
    # for each file
    peak_height_to_width_vec = Vector{Measurements.Measurement{Float64}}(undef, 
        length(r_t_files))
    beta_t_max_t_gradient_vec = []
    
    for (i, file) in enumerate(r_t_files)
        # get the peak to width metric
        r_t_dict = GU.load_h5_dict(file)
        peak_height_to_width = get_reflection_peak_height_to_width(r_t_dict, 
            p0=p0,
            upper_bounds=upper_bounds,
            lower_bounds=lower_bounds,
            save_plot=false)
        peak_height_to_width_vec[i] = peak_height_to_width

        # get the beta, t_max, and t_gradient values from the filename
        filename = basename(file)
        parts = split(filename, '_')
        push!(beta_t_max_t_gradient_vec, (parse(Float64, parts[3]), 
        parse(Float64, parts[6]), parse(Float64, parts[9])))

        if print_progress && i%10 == 0
            progress_percentage = i/length(r_t_files) * 100
            println("Processed fits: $progress_percentage%")
        end
    end

    # equally, get the combinations of beta, t_max, and t_gradient from the 
    # order metrics dict
    order_metric_beta_t_max_t_gradient_vec = [
        (order_metric_dict["bond_bending_const_vec"][i], 
        order_metric_dict["t_max_vec"][i], 
        order_metric_dict["t_gradient_vec"][i]) 
        for i in eachindex(order_metric_dict["bond_bending_const_vec"])]

    # find the indices where the entries of the
    # order_metric_beta_t_max_t_gradient_vec are the sames as the entries of 
    # the beta_t_max_t_gradient_vec
    for i in eachindex(beta_t_max_t_gradient_vec)
        if !(beta_t_max_t_gradient_vec[i] 
            in order_metric_beta_t_max_t_gradient_vec)
            @warn "The combination $(beta_t_max_t_gradient_vec[i]) 
            is not in the order metric dictionary."
            println(basename(r_t_files[i]))
        end
    end

    indices = [findfirst(x -> x == beta_t_max_t_gradient_vec[i],
        order_metric_beta_t_max_t_gradient_vec) 
        for i in eachindex(beta_t_max_t_gradient_vec)]

    # loop through keys of the dictionary which are all vectors and add them to
    # a list of vectors
    order_metric_vec_vec = []
    key_vec = []
    current_index = 1
    for key in keys(order_metric_dict)
        if (key == "filenames_vec" || key == "bond_bending_const_vec" 
            || key == "t_max_vec" || key == "t_gradient_vec" 
            || key == "total_keating_energy_vec")
            continue
        elseif key == "q_l_mat"
            for j in 1:13
                push!(order_metric_vec_vec, Measurements.value.(
                    order_metric_dict[key][j,indices]))
            end
            q_l_key_vec = ["q_$(i)_mean" for i in 0:12]
            append!(key_vec, q_l_key_vec)
            current_index += 13
            for j in 1:13
                push!(order_metric_vec_vec, Measurements.uncertainty.(
                    order_metric_dict[key][j,indices]))
            end
            q_l_key_vec = ["q_$(i)_std" for i in 0:12]
            append!(key_vec, q_l_key_vec)
            current_index += 13
        else
            push!(order_metric_vec_vec, order_metric_dict[key][indices])
            push!(key_vec, key)
            current_index += 1
        end
    end

    # using Pearson's correlation coefficient, calculate the correlation 
    # between the peak height to width metrics and the order metrics
    correlations = [Statistics.cor( 
        Measurements.value.(peak_height_to_width_vec), 
        order_metric_vec_vec[i]) for i in 1:length(order_metric_vec_vec)]

    # sort keys according to the absolute values of their correlation 
    # coefficients
    sorted_indices = sortperm(abs.(correlations), rev=true)
    sorted_keys = key_vec[sorted_indices]
    sorted_correlations = correlations[sorted_indices]

    # Print the correlation coefficients
    for (i, key) in enumerate(sorted_keys)
        println("Pearson corr. for $key: ", sorted_correlations[i])
    end

    # create a dictionary with the keys being the order metric names and the
    # values being the correlation coefficients
    correlations_dict = Dict(zip(sorted_keys, sorted_correlations))

    if save_results
        GU.save_dict_to_h5(copy(correlations_dict),
            analysis_data_path
            *"order_peak_height_to_width_pearson_correlations.h5")
    end

    spearman_correlations = [StatsBase.corspearman(peak_height_to_width_vec, 
        order_metric_vec_vec[i]) for i in 1:length(order_metric_vec_vec)]

    # sort keys according to the absolute values of their Spearman's rank
    # correlation coefficients
    sorted_spearman_indices = sortperm(abs.(spearman_correlations), rev=true)
    sorted_spearman_keys = key_vec[sorted_spearman_indices]
    sorted_spearman_correlations = (
        spearman_correlations[sorted_spearman_indices])

    # Print the Spearman's rank correlation coefficients
    for (i, key) in enumerate(sorted_spearman_keys)
        println("Spearman corr. for $key: ", sorted_spearman_correlations[i])
    end

    # create a dictionary with the keys being the order metric names and the
    # values being the Spearman's rank correlation coefficients
    spearman_correlations_dict = Dict(zip(sorted_spearman_keys, 
        sorted_spearman_correlations))

    if save_results
        GU.save_dict_to_h5(copy(spearman_correlations_dict),
            analysis_data_path*"order_peak_height_to_width_spearman_correlations.h5")
    end

    return [correlations_dict, spearman_correlations_dict]
end