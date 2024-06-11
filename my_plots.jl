
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex
import Measurements
import Polynomials


path = raw"..\..\presentations\material\\"

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


filenames = ["216_vertices_T_0.1_heated_for_1.0_steps_quenched",
"216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched",
"216_vertices_T_0.3_heated_for_10.0_steps_quenched",
"216_vertices_T_0.4_heated_for_0.25_steps_quenched",
"216_vertices_T_0.25_heated_for_0.05_steps_quenched",
"216_vertices_T_2.0_heated_for_0.5_steps_quenched"
]

for filename in filenames
    data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename
    save_path = raw"..\plots\random_networks\\"*filename

    autocovariance_fct_direction_dict = GU.load_h5_dict(data_path*"_autocovariance_fct_direction.h5")

    NA.plot_autocovariance_fct_heatmap(autocovariance_fct_direction_dict,
        save_path;
        save_plot = true,
        clims = nothing,
        x_y_lims = nothing,
        sampling_vector_component_to_fix = 3,
        sampling_vector_value_fixed = 0)


    spectral_density_dict = GU.load_h5_dict(data_path*"_spectral_density_array.h5")

    NA.plot_spectral_density_heatmap(spectral_density_dict,
        save_path;
        save_plot = true,
        clims = (0,0.1),
        x_y_lims = nothing,
        wavevector_component_to_fix = 3,
        wavevector_value_fixed = 0)

    spectral_density_dict = GU.load_h5_dict(data_path*"_spectral_density_array.h5")

    spectral_density_angle_averaged_dict = GU.load_h5_dict(data_path*"_spectral_density_angle_averaged.h5")

    Plots.plot(spectral_density_angle_averaged_dict["wavenumber_vec"], 
                        Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"]) , 
                        ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"]))

    Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", legend = false, xlims=(0,15), ylims=(0,minimum([1000.0, maximum(Measurements.value.(spectral_density_angle_averaged_dict["unfiltered_spectral_density_vec"]))])), xtick=pitick(0, 15, 1; mode=:latex))

    Plots.savefig(raw"..\plots\random_networks\\"*filename*"_spectral_density_angle_averaged.png")

    Plots.plot(spectral_density_angle_averaged_dict["unfiltered_wavenumber_vec"], 
    Measurements.value.(spectral_density_angle_averaged_dict["unfiltered_spectral_density_vec"]),
    ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["unfiltered_spectral_density_vec"]))

    Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", legend = false, xlims=(0,15), ylims=(0, minimum([1000.0, maximum(Measurements.value.(spectral_density_angle_averaged_dict["unfiltered_spectral_density_vec"]))])), xtick=pitick(0, 15, 1; mode=:latex))

    Plots.savefig(raw"..\plots\random_networks\\"*filename*"_spectral_density_angle_averaged_unfiltered.png")


    volume_fract_variance_dict = NA.get_volume_fract_variance(autocovariance_fct_direction_dict;
            save_result = false)

    # plot the volume fraction variance
    Plots.plot(volume_fract_variance_dict["sphere_radius_vec"], volume_fract_variance_dict["volume_fract_variance_times_window_volume_vec"], xlabel="window radius "*Latex.L"R / d", ylabel=Latex.L"\sigma_V^2(R) \cdot v_1(R)", xlims=(0, maximum(volume_fract_variance_dict["sphere_radius_vec"])), ylims=(0, maximum(volume_fract_variance_dict["volume_fract_variance_times_window_volume_vec"])), legend=false)

    Plots.savefig(raw"..\plots\random_networks\\"*filename*"_volume_fraction_variance.png")
end