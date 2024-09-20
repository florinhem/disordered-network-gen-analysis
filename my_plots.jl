
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex
import Measurements
import Polynomials

fontsize=18

Plots.gr()
Plots.default(grid=false, 
legend = true, 
dpi=250,
xtickfontsize=fontsize,
ytickfontsize=fontsize,
xguidefontsize=fontsize,
yguidefontsize=fontsize,
legendfontsize=fontsize,
bottom_margin = 3Plots.mm,
linewidth=3, 
thickness_scaling = 1,
framestyle = :box)

# functions to have pi ticks
function pitick(start, stop, denom; mode=:text)
    a = Int(cld(start, 2*π/denom))
    b = Int(fld(stop, 2*π/denom))
    tick = range(a*2*π/denom, b*2*π/denom; step=2*π/denom)
    ticklabel = piticklabel.( 2 .* (a:b) .// denom, Val(mode))
    tick, ticklabel
end

function piticklabel(x::Rational, ::Val{:text})
    iszero(x) && return "0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return S * N * "π"
    S * N * "π/" * repr(d)
end

function piticklabel(x::Rational, ::Val{:latex})
    iszero(x) && return Latex.L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return Latex.L"%$S%$N\pi"
    Latex.L"%$S\frac{%$N\pi}{%$d}"
end


"""
    get_tickslogscale(lims; skiplog=false)
Return a tuple (ticks, ticklabels) for the axis limit `lims`
where multiples of 10 are major ticks with label and minor ticks have no label
skiplog argument should be set to true if `lims` is already in log scale.
"""
function get_tickslogscale(lims::Tuple{T, T}; skiplog::Bool=false) where {T<:AbstractFloat}
    mags = if skiplog
        # if the limits are already in log scale
        floor.(lims)
    else
        floor.(log10.(lims))
    end
    rlims = if skiplog; 10 .^(lims) else lims end

    total_tickvalues = []
    total_ticknames = []

    rgs = range(mags..., step=1)
    for (i, m) in enumerate(rgs)
        if m >= 0
            tickvalues = range(Int(10^m), Int(10^(m+1)); step=Int(10^m))
            ticknames  = vcat([string(round(Int, 10^(m)))],
                              ["" for i in 2:9],
                              [string(round(Int, 10^(m+1)))])
        else
            tickvalues = range(10^m, 10^(m+1); step=10^m)
            ticknames  = vcat([string(10^(m))], ["" for i in 2:9], [string(10^(m+1))])
        end

        if i==1
            # lower bound
            indexlb = findlast(x->x<rlims[1], tickvalues)
            if isnothing(indexlb); indexlb=1 end
        else
            indexlb = 1
        end
        if i==length(rgs)
            # higher bound
            indexhb = findfirst(x->x>rlims[2], tickvalues)
            if isnothing(indexhb); indexhb=10 end
        else
            # do not take the last index if not the last magnitude
            indexhb = 9
        end

        total_tickvalues = vcat(total_tickvalues, tickvalues[indexlb:indexhb])
        total_ticknames = vcat(total_ticknames, ticknames[indexlb:indexhb])
    end
    return (total_tickvalues, total_ticknames)
end

"""
    fancylogscale!(p; forcex=false, forcey=false)
Transform the ticks to log scale for the axis with scale=:log10.
forcex and forcey can be set to true to force the transformation
if the variable is already expressed in log10 units.
"""
function fancylogscale!(p::Plots.Subplot; forcex::Bool=false, forcey::Bool=false)
    kwargs = Dict()
    for (ax, force, lims) in zip((:x, :y), (forcex, forcey), (Plots.xlims, Plots.ylims))
        axis = Symbol("$(ax)axis")
        ticks = Symbol("$(ax)ticks")

        if force || p.attr[axis][:scale] == :log10
            # Get limits of the plot and convert to Float
            ls = float.(lims(p))
            ts = if force
                (vals, labs) = get_tickslogscale(ls; skiplog=true)
                (log10.(vals), labs)
            else
                get_tickslogscale(ls)
            end
            kwargs[ticks] = ts
        end
    end

    if length(kwargs) > 0
        Plots.plot!(p; kwargs...)
    end
    p
end
fancylogscale!(p::Plots.Plot; kwargs...) = (fancylogscale!(p.subplots[1]; kwargs...); return p)
fancylogscale!(; kwargs...) = fancylogscale!(Plots.plot!(); kwargs...)


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\1000_vertices_bond_bending_0.285\run_1\\"

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.285\run_1\\"

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[1:6]

Plots.plot()

for temperature in temperature_vec
    volume_fraction_variance_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_volume_fraction_variance.h5")

    Plots.plot!(volume_fraction_variance_dict["sphere_radius_vec"], Measurements.value.(volume_fraction_variance_dict["volume_fract_variance_times_window_volume_vec"]), 
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "window radius / "*Latex.L"d^{-1}", ylabel = Latex.L"\sigma_V^2 \cdot"*"window volume")

Plots.savefig(path*"volume_fraction_variance_times_window_volume_low_t.png")

Plots.plot()

for temperature in temperature_vec
    volume_fraction_variance_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_volume_fraction_variance.h5")

    Plots.plot!(volume_fraction_variance_dict["sphere_radius_vec"], Measurements.value.(volume_fraction_variance_dict["volume_fract_variance_vec"]), 
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "window radius / "*Latex.L"d^{-1}", ylabel = Latex.L"\sigma_V^2")

Plots.savefig(path*"volume_fraction_variance_low_t.png")

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[7:end]

Plots.plot()

for temperature in temperature_vec
    volume_fraction_variance_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_volume_fraction_variance.h5")

    Plots.plot!(volume_fraction_variance_dict["sphere_radius_vec"], Measurements.value.(volume_fraction_variance_dict["volume_fract_variance_times_window_volume_vec"]), 
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "window radius / "*Latex.L"d^{-1}", ylabel = Latex.L"\sigma_V^2 \cdot"*"window volume")

Plots.savefig(path*"volume_fraction_variance_times_window_volume_high_t.png")


Plots.plot()

for temperature in temperature_vec
    volume_fraction_variance_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_volume_fraction_variance.h5")

    Plots.plot!(volume_fraction_variance_dict["sphere_radius_vec"], Measurements.value.(volume_fraction_variance_dict["volume_fract_variance_vec"]), 
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "window radius / "*Latex.L"d^{-1}", ylabel = Latex.L"\sigma_V^2")

Plots.savefig(path*"volume_fraction_variance_high_t.png")




plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\heat_cool_bond_bending_0.285\\"

diamonds_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

diamond_order_metrics_dicts = [GU.load_h5_dict(diamonds_analysis_data_path*"216_vertices_perfect_diamond_small_scale_order_metrics.h5"),
GU.load_h5_dict(diamonds_analysis_data_path*"512_vertices_perfect_diamond_small_scale_order_metrics.h5"),
GU.load_h5_dict(diamonds_analysis_data_path*"1000_vertices_perfect_diamond_small_scale_order_metrics.h5")]


analysis_data_paths = [
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285_heat_cool\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.285\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.285\\"]

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "cluster_metric_vec", 
"pore_size_distribution_second_moment_vec",
"anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec"]


order_metrics_labels = ["Bond length std", "Bond angle std", "Dihedral angle std", "Cluster metric",
"Pore size dist. 2nd m.",
"Anisotropy from s. f.", "Anisotropy from s. d."]

order_metrics_dicts = []

# loop through folders and append all order metrics to the order_metrics_dict
for analysis_data_path in analysis_data_paths 

    order_metrics_dict = Dict()

    for i in 1:5

        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"

        current_order_metrics_dict = GU.load_h5_dict(current_analysis_data_path*"all_order_metrics.h5")

        for (key, value) in current_order_metrics_dict
            if haskey(order_metrics_dict, key)
                order_metrics_dict[key] = vcat(order_metrics_dict[key], value)
            else
                order_metrics_dict[key] = value
            end
        end

    end

    # sort all vectors in order of the total keating energy
    for order_metric_name in order_metrics_names
        order_metrics_dict[order_metric_name] = order_metrics_dict[order_metric_name][sortperm(order_metrics_dict["total_keating_energy_vec"])]
    end
    order_metrics_dict["filenames_vec"] = order_metrics_dict["filenames_vec"][sortperm(order_metrics_dict["total_keating_energy_vec"])]
    sort!(order_metrics_dict["total_keating_energy_vec"])

    push!(order_metrics_dicts, order_metrics_dict)
    
end

markershapes = [:circle, :rect, :diamond]

nr_vertices = [216, 512, 1000]

for i in eachindex(order_metrics_names)

    Plots.scatter()

    for j in 1:3


        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dicts[j]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = minimum(temperatures)
        max_temp = maximum(temperatures)
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        if order_metrics_names[i] == "anisotropy_metric_from_structure_factor_vec" || order_metrics_names[i] == "anisotropy_metric_from_spectral_density_vec"
            Plots.scatter!(order_metrics_dicts[j]["total_keating_energy_vec"] ./ nr_vertices[j], order_metrics_dicts[j][order_metrics_names[i]] ./ diamond_order_metrics_dicts[j][order_metrics_names[i][1:end-4]],
            color = mapped_colors, markershape = markershapes[j], yscale = :log10)
        else
            
            Plots.scatter!([0], [diamond_order_metrics_dicts[j][order_metrics_names[i][1:end-4]]], markershape = markershapes[j], color = :black)

            Plots.scatter!(order_metrics_dicts[j]["total_keating_energy_vec"] ./ nr_vertices[j], order_metrics_dicts[j][order_metrics_names[i]],
            color = mapped_colors, markershape = markershapes[j])
        end

    end

    if order_metrics_names[i] == "pore_size_distribution_second_moment_vec"
        #Plots.scatter!( yscale = :log10)
        #fancylogscale!()
        Plots.scatter!(yscale = :lin)
    else
        Plots.scatter!(yscale = :lin)
        
    end

    Plots.scatter!(legend = false, ylabel = order_metrics_labels[i], xlabel = "Keating energy per vertex")
    Plots.savefig(plots_save_path*order_metrics_names[i]*".png")
end

