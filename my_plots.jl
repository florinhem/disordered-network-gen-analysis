
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

path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_globally_relaxed\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\\"

order_metrics_dict = Dict()

# loop through folders and append all order metrics to the order_metrics_dict
for i in 1:5
    
    current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_"*string(i)*"\\"

    current_order_metrics_dict = GU.load_h5_dict(current_analysis_data_path*"all_order_metrics.h5")

    for (key, value) in current_order_metrics_dict
        if haskey(order_metrics_dict, key)
            order_metrics_dict[key] = vcat(order_metrics_dict[key], value)
        else
            order_metrics_dict[key] = value
        end
    end
    
end

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec", "cluster_metric_vec"]

# sort all vectors in order of the total keating energy
for order_metric_name in order_metrics_names
    order_metrics_dict[order_metric_name] = order_metrics_dict[order_metric_name][sortperm(order_metrics_dict["total_keating_energy_vec"])]
end
order_metrics_dict["filenames_vec"] = order_metrics_dict["filenames_vec"][sortperm(order_metrics_dict["total_keating_energy_vec"])]
sort!(order_metrics_dict["total_keating_energy_vec"])

# now start the happy plotting

y_labels = ["bond length std", "bond angle std", "dihedral angle std", "anisotropy in structure f.", "anisotropy in spectral d.", "cluster metric"]

for i in eachindex(order_metrics_names)
    Plots.scatter(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216, order_metrics_dict[order_metrics_names[i]][1:end-10], xlabel = "Keating energy per vertex", ylabel = y_labels[i], legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
    ylims = (minimum(order_metrics_dict[order_metrics_names[i]][1:end-10]), maximum(order_metrics_dict[order_metrics_names[i]][1:end-10])))
    Plots.savefig(path*order_metrics_names[i][1:end-4]*".png")
end

mask_vec = [contains.(order_metrics_dict["filenames_vec"], "heated_for_0.1_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_0.25_steps"),
    contains.(order_metrics_dict["filenames_vec"], "heated_for_0.5_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_1.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_5.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_10.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.025"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.05"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.1"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.2"),
( contains.(order_metrics_dict["filenames_vec"], "cool_0.1")
    .& .!(contains.(order_metrics_dict["filenames_vec"], "heat")) ),
]

filename_vec = ["heated_for_0.1_steps", "heated_for_0.25_steps", "heated_for_0.5_steps", "heated_for_1.0_steps", "heated_for_5.0_steps", "heated_for_10.0_steps", "heat_cool_0.025", "heat_cool_0.05", "heat_cool_0.1", "heat_cool_0.2", "cool_0.1"]

title_vec = ["heated for 0.1 steps", "heated for 0.25 steps", "heated for 0.5 steps", "heated for 1.0 steps", "heated for 5.0 steps", "heated for 10.0 steps", "heat and cool 0.025/MC step", "heat and cool 0.05/MC step", "heat and cool 0.1/MC step", "heat and cool 0.2/MC step", "cool 0.1/MC step"]

for i in eachindex(mask_vec)
    mask = mask_vec[i]
    filtered_filenames_vec = order_metrics_dict["filenames_vec"][mask]
    filtered_total_keating_energy_vec = order_metrics_dict["total_keating_energy_vec"][mask]
    filtered_bond_length_std_vec = order_metrics_dict["bond_length_std_vec"][mask]
    filtered_bond_angle_std_vec = order_metrics_dict["bond_angle_std_vec"][mask]
    filtered_dihedral_angle_std_vec = order_metrics_dict["dihedral_angle_std_vec"][mask]
    filtered_anisotropy_metric_from_structure_factor_vec = order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][mask]
    filtered_anisotropy_metric_from_spectral_density_vec = order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][mask]
    filtered_cluster_metric_vec = order_metrics_dict["cluster_metric_vec"][mask]

    # get the temperature from the filtered filenames
    pattern = r"T_([0-9\.]+)"
    extracted_numbers = [match(pattern, s).captures[1] for s in filtered_filenames_vec]
    temperatures = parse.(Float64, extracted_numbers)

    min_temp = minimum(temperatures)
    max_temp = maximum(temperatures)
    normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
    colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
    mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]


    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_bond_length_std_vec, xlabel = "Keating energy per vertex", ylabel = "bond length std", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["bond_length_std_vec"][1:end-10]), maximum(order_metrics_dict["bond_length_std_vec"][1:end-10])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_length_std.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_bond_angle_std_vec, xlabel = "Keating energy per vertex", ylabel = "bond angle std", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["bond_angle_std_vec"][1:end-10]), maximum(order_metrics_dict["bond_angle_std_vec"][1:end-10])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_angle_std.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_dihedral_angle_std_vec, xlabel = "Keating energy per vertex", ylabel = "dihedral angle std", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["dihedral_angle_std_vec"][1:end-10]), maximum(order_metrics_dict["dihedral_angle_std_vec"][1:end-10])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_dihedral_angle_std.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_anisotropy_metric_from_structure_factor_vec, xlabel = "Keating energy per vertex", ylabel = "anisotropy in structure f.", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][1:end-10]), maximum(order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][1:end-10])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_anisotropy_in_structure_factor.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_anisotropy_metric_from_spectral_density_vec, xlabel = "Keating energy per vertex", ylabel = "anisotropy in spectral d.", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][1:end-10]), maximum(order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][1:end-10])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_anisotropy_in_spectral_density.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_cluster_metric_vec, xlabel = "Keating energy per vertex", ylabel = "cluster metric", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["cluster_metric_vec"][1:end-10]), maximum(order_metrics_dict["cluster_metric_vec"][1:end-10])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_cluster_metric.png")
end