include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import Random
import LaTeXStrings as Latex
import Measurements
import Polynomials
import FFTW
import Statistics


path = raw"..\..\presentations\material\\"

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

x = collect(0:0.01:10)
plot_1 = 1/2 * exp.( .- x ./2 )
plot_2 = 1/4 * exp.( .- x ./4 )

my_plot = Plots.plot(x, plot_1, label = Latex.L"T_1")
my_plot = Plots.plot!(x, plot_2, label = Latex.L"T_2>T_1")
my_plot = Plots.plot!(xlabel="energy",
ylabel = "probability")

Plots.savefig(path*"boltzmann_distribution_2.png")


x = collect(-8:0.01:8)

plot_1 = min.(1, exp.( .-  x ./ 1  ) )  
plot_2 = min.(1, exp.( .-  x ./ 4  ) )  

myplot = Plots.plot(x, plot_2, label = Latex.L"kT=4", linecolor=Plots.palette(:tab10)[2])
myplot = Plots.plot!(x, plot_1, label = Latex.L"kT=1", linecolor=Plots.palette(:tab10)[1])

Plots.plot!(grid=false, xlabel="energy difference",
ylabel = "acceptance probability")

Plots.savefig(path*"boltzmann_metropolis_3.png")


random_array = rand(Float64, (2, 100))

my_plot = Plots.plot(random_array[1,:], random_array[2,:], seriestype=:scatter, aspect_ratio=:equal, fillcolor=Plots.palette(:tab10)[5], markercolor=Plots.palette(:tab10)[5])
my_plot = Plots.plot!(xlabel=Latex.L"x", ylabel=Latex.L"y",
legend = false, dpi=400, xlims=(0,1), ylims=(0,1))

Plots.savefig(path*"poisson_process.png")


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_1 = NG.load_graph_from_h5_and_MGformat(dict_path*"1000_vertices_T_1_quenched")
graph_dict_4 = NG.load_graph_from_h5_and_MGformat(dict_path*"1000_vertices_T_4_quenched")


bond_length_std_1, bond_length_vec_1 = NA.get_bond_length_std(graph_dict_1)
bond_length_std_4, bond_length_vec_4 = NA.get_bond_length_std(graph_dict_4)

b_range = range(0.4, 1.8, length=71)

my_plot = Plots.stephist(bond_length_vec_4, bins=b_range, label = Latex.L"kT=4", normalize=:probability, linecolor=Plots.palette(:tab10)[2])
my_plot = Plots.stephist!(bond_length_vec_1, bins=b_range, label = Latex.L"kT=1", normalize=:probability, linecolor=Plots.palette(:tab10)[1])
my_plot = Plots.plot!(xlabel="bond length",
ylabel = "relative frequency")

Plots.savefig(path*"bond_length_T_1_4.png")


bond_angle_std_1, bond_angle_vec_1 = NA.get_bond_angle_std(graph_dict_1)
bond_angle_std_4, bond_angle_vec_4 = NA.get_bond_angle_std(graph_dict_4)

b_range = range(0, 180, length=61)

my_plot = Plots.stephist(bond_angle_vec_4 ./pi .* 180, bins=b_range, label = Latex.L"kT=4", normalize=:probability, linecolor=Plots.palette(:tab10)[2])
my_plot = Plots.stephist!(bond_angle_vec_1 ./pi .* 180, bins=b_range, label = Latex.L"kT=1", normalize=:probability, linecolor=Plots.palette(:tab10)[1])
my_plot = Plots.plot!(xlabel="bond angle",
ylabel = "relative frequency")

Plots.savefig(path*"bond_angle_T_1_4.png")


steinhardt_order_parameter_dict_1 = NA.get_steinhardt_order_parameter_dict(graph_dict_1, 6)
steinhardt_order_parameter_dict_4 = NA.get_steinhardt_order_parameter_dict(graph_dict_4, 6)

diamond = [0.509, 0.629]
cubic = [0.764, 0.354]
fcc = [0.191, 0.575]
hcp = [0.097, 0.485]

marker_size = 9

my_plot = Plots.plot([steinhardt_order_parameter_dict_4[4]], [steinhardt_order_parameter_dict_4[6]], seriestype=:scatter, label="kT=4", mc=Plots.palette(:tab10)[2], ms=marker_size)
my_plot = Plots.plot!([steinhardt_order_parameter_dict_1[4]], [steinhardt_order_parameter_dict_1[6]], seriestype=:scatter, label="kT=1", mc=Plots.palette(:tab10)[1], ms=marker_size)
my_plot = Plots.plot!([diamond[1]], [diamond[2]], seriestype=:scatter, label="kT=0", ms=marker_size)
my_plot = Plots.plot!([cubic[1]], [cubic[2]], seriestype=:scatter, label="cubic", ms=marker_size)
my_plot = Plots.plot!([fcc[1]], [fcc[2]], seriestype=:scatter, label="fcc", ms=marker_size)
my_plot = Plots.plot!([hcp[1]], [hcp[2]], seriestype=:scatter, label="hcp", ms=marker_size)
my_plot = Plots.plot!(grid=true, legend=false, xlabel=Latex.L"q_4",
ylabel = Latex.L"q_6")

Plots.savefig(path*"steinhardt_order_parameter_T_1_4.png")


x_vec = collect(0:0.01:4)
plot_1 = (0.8 .* sinc.(( 2 .* (x_vec .- 2.5)) .^2) .+ 0.2) .* exp.( .-(x_vec .- 2.5).^2)
Plots.plot(x_vec, plot_1, linecolor="black", ylims=(0,1), xlims=(0,4), framestyle = :box, legend=false)

Plots.savefig(path*"reflectivity_band_structure_2.png")


x_vec = collect(0:0.01:2)
y_vec = collect(0.5:0.01:1.5)
Plots.plot(x_vec, (3/16) .* (x_vec.^2 .- 1).^2   )
Plots.plot!(y_vec, (3/4) .* (y_vec .- 1).^2 , linestyle=:dash )
Plots.plot!(xlabel="bond length / "*Latex.L"d", ylabel="energy", right_margin = 3Plots.mm, ylims=(0,0.4), xlims=(0,2), legend=false)

Plots.savefig(path*"bond_stretching_energy.png")


x_vec = collect(0:0.1:180)
y_vec = collect(40:0.1:180)
Plots.plot(x_vec, (3/8 * 0.285) .* (cosd.(x_vec) .+ 1/3).^2  )
Plots.plot!(y_vec, (0.095 .* (y_vec .* pi ./ 180 .- acos(-1/3)).^2 ) , linestyle=:dash  )
Plots.plot!(xlabel="bond angle / °", ylabel="energy", right_margin = 5Plots.mm, ylims=(0,0.4), xlims=(0,180), legend=false)
Plots.savefig(path*"bond_bending_energy.png")


x_vec = collect(0:0.01:2)
y_vec = collect(0.5:0.01:1.5)
Plots.plot(x_vec, (3/16) .* (x_vec.^2 .- 1).^2   )
Plots.plot!(y_vec, (3/4) .* (y_vec .- 1).^2 , linestyle=:dash )
Plots.plot!(xlabel="bond length / "*Latex.L"d", ylabel="energy", right_margin = 4Plots.mm, bottom_margin = 8Plots.mm, left_margin = 4Plots.mm, ylims=(0,0.15), xlims=(0.7,1.3), legend=false, size = (650, 300), xticks=collect(0.7:0.1:1.3))

Plots.savefig(path*"bond_stretching_energy_zoom.png")

x_vec = collect(40:0.1:180)
y_vec = collect(40:0.1:180)
Plots.plot(x_vec, (3/8 * 0.285) .* (cosd.(x_vec) .+ 1/3).^2  )
Plots.plot!(y_vec, (0.095 .* (y_vec .* pi ./ 180 .- acos(-1/3)).^2 ) , linestyle=:dash  )
Plots.plot!(xlabel="bond angle / °", ylabel="energy",  right_margin = 5Plots.mm, bottom_margin = 8Plots.mm, left_margin = 4Plots.mm, ylims=(0,0.15), xlims=(40,180), legend=false, size = (650, 300), xticks=collect(40:20:180))
Plots.savefig(path*"bond_bending_energy_zoom.png")


dict_path_1 = raw"..\analysis_data\random_networks\1000_vertices_T_1_quenched_structure_factor_bartlett_isotrope.h5"

dict_path_4 = raw"..\analysis_data\random_networks\1000_vertices_T_4_quenched_structure_factor_bartlett_isotrope.h5"

structure_factor_dict_1 = GU.load_h5_dict(dict_path_1)
structure_factor_dict_4 = GU.load_h5_dict(dict_path_4)

my_plot = Plots.plot(structure_factor_dict_4["wavenumber_vec"], structure_factor_dict_4["structure_factor_vec"], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=4" )
my_plot = Plots.plot!(structure_factor_dict_1["wavenumber_vec"], structure_factor_dict_1["structure_factor_vec"], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=1" )
my_plot = Plots.plot!(xlabel="wavenumber",
ylabel = "structure factor", xlims=(0,32.5), ylims=(0,2.75), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_T_1_4.png")


# get structure factor dicts
dict_path_1 = raw"..\analysis_data\random_networks\1000_vertices_T_1_quenched_structure_factor_bartlett_isotrope.h5"
dict_path_4 = raw"..\analysis_data\random_networks\1000_vertices_T_4_quenched_structure_factor_bartlett_isotrope.h5"

structure_factor_dict_1 = GU.load_h5_dict(dict_path_1)
structure_factor_dict_4 = GU.load_h5_dict(dict_path_4)

# get effective hyperuniformity parameter and fit parameters for T=1 and T=4
hyperuniformity_parameter_1, polynomial_fit_1 = NA.get_hyperuniformity_metric(structure_factor_dict_1)
hyperuniformity_parameter_4, polynomial_fit_4 = NA.get_hyperuniformity_metric(structure_factor_dict_4)

# plot structure factor
x_vec = collect(0:10/200:10)
fit_1_vec = polynomial_fit_1.(x_vec)
fit_4_vec = polynomial_fit_4.(x_vec)

Plots.plot(structure_factor_dict_4["wavenumber_vec"], structure_factor_dict_4["structure_factor_vec"], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=4" )
Plots.plot!(structure_factor_dict_1["wavenumber_vec"], structure_factor_dict_1["structure_factor_vec"], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=1" )
Plots.plot!(x_vec, fit_4_vec, linecolor="sienna", ls=:dash, label = "fit "*Latex.L"kT=4")
Plots.plot!(x_vec, fit_1_vec, linecolor="cyan2", ls=:dash, label = "fit "*Latex.L"kT=1")
Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,22.5), ylims=(0,20), size = (500, 600), bottom_margin = 0Plots.mm, xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_T_1_4_fits.png")


Plots.plot(structure_factor_dict_4["wavenumber_vec"], structure_factor_dict_4["structure_factor_vec"], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=4" )
Plots.plot!(structure_factor_dict_1["wavenumber_vec"], structure_factor_dict_1["structure_factor_vec"], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=1" )
Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,22.5), ylims=(0,20), size = (500, 600), bottom_margin = 0Plots.mm, xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_T_1_4_stretched.png")


# load structure factor dictionaries
dict_path = raw"..\analysis_data\random_networks\\"

filenames = ["1000_vertices_T_1_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5",
    "512_vertices_T_1_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5",
    "216_vertices_T_1_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5"]

structure_factor_dict_1 = GU.load_h5_dict(dict_path*filenames[1])
structure_factor_dict_5 = GU.load_h5_dict(dict_path*filenames[2])
structure_factor_dict_2 = GU.load_h5_dict(dict_path*filenames[3])


# plot structure factors 
Plots.plot(
    structure_factor_dict_2["wavenumber_vec"][6:end],
    structure_factor_dict_2["structure_factor_vec"][6:end],
    label = "216 vertices"
)
Plots.plot!(
    structure_factor_dict_5["wavenumber_vec"][6:end],
    structure_factor_dict_5["structure_factor_vec"][6:end], ls=:dash,
    label = "512 vertices"
)
Plots.plot!(
    structure_factor_dict_1["wavenumber_vec"][6:end],
    structure_factor_dict_1["structure_factor_vec"][6:end], ls=:dash,
    label = "1000 vertices"
)
Plots.plot!(
    xlabel = "wavenumber",
    ylabel = "structure factor",
    legend = true,
    xlims=(0,22.5), ylims=(0,10), xtick=pitick(0, 32, 1; mode=:latex)
)

# save plot
Plots.savefig(path*"structure_factor_T_1_size_comparison.png")


# load structure factor dictionaries
dict_path = raw"..\analysis_data\random_networks\\"

filenames = ["1000_vertices_T_1_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5",
    "1000_vertices_T_4_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5",
    "1000_vertices_T_1_quenched_structure_factor_isotrope.h5",
    "1000_vertices_T_4_quenched_structure_factor_isotrope.h5"]

structure_factor_bartlett_dict_1 = GU.load_h5_dict(dict_path*filenames[1])
structure_factor_bartlett_dict_4 = GU.load_h5_dict(dict_path*filenames[2])
structure_factor_dict_1 = GU.load_h5_dict(dict_path*filenames[3])
structure_factor_dict_4 = GU.load_h5_dict(dict_path*filenames[4])


window_size = length(structure_factor_dict_1["wavenumber_vec"])/100.0

# plot structure factors 
Plots.plot(
    structure_factor_bartlett_dict_1["wavenumber_vec"],
    structure_factor_bartlett_dict_1["structure_factor_vec"],
    label = "Bartlett "*Latex.L"kT=1"
)
Plots.plot!(
    structure_factor_bartlett_dict_4["wavenumber_vec"],
    structure_factor_bartlett_dict_4["structure_factor_vec"],
    label ="Bartlett "*Latex.L"kT=4"
)
Plots.plot!(
    structure_factor_dict_1["wavenumber_vec"],
    NaNStatistics.movmean(structure_factor_dict_1["structure_factor_vec"], window_size),
    ls = :dash,
    label = "S. I. "*Latex.L"kT=1"
)
Plots.plot!(
    structure_factor_dict_4["wavenumber_vec"],
    NaNStatistics.movmean(structure_factor_dict_4["structure_factor_vec"], window_size),
    ls = :dash,
    label = "S. I. "*Latex.L"kT=4"
)
Plots.plot!(
    xlabel = "wavenumber",
    ylabel = "structure factor",
    legend = true,
    xlims=(0,22.5), ylims=(0,10), xtick=pitick(0, 32, 1; mode=:latex)
)

# save plot
Plots.savefig(path*"structure_factor_scattering_intensity_bartlett_comparison.png")


# get structure factor dicts
dict_path_1 = raw"..\analysis_data\random_networks\1000_vertices_T_1_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5"
dict_path_4 = raw"..\analysis_data\random_networks\1000_vertices_T_4_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5"

structure_factor_dict_1 = GU.load_h5_dict(dict_path_1)
structure_factor_dict_4 = GU.load_h5_dict(dict_path_4)

Plots.plot(structure_factor_dict_4["wavenumber_vec"][10:end], structure_factor_dict_4["structure_factor_vec"][10:end], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=4" )
Plots.plot!(structure_factor_dict_1["wavenumber_vec"][10:end], structure_factor_dict_1["structure_factor_vec"][10:end], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=1" )
Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,22.5), ylims=(0,10), size = (500, 600), bottom_margin = 0Plots.mm, xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_T_1_4_stretched_high_sampling_rate.png")


# get structure factor dicts
dict_path_1 = raw"..\analysis_data\random_networks\1000_vertices_T_1_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5"
dict_path_4 = raw"..\analysis_data\random_networks\1000_vertices_T_4_quenched_high_sampling_rate_structure_factor_bartlett_isotrope.h5"

structure_factor_dict_1 = GU.load_h5_dict(dict_path_1)
structure_factor_dict_4 = GU.load_h5_dict(dict_path_4)

Plots.plot(structure_factor_dict_4["wavenumber_vec"][10:end], structure_factor_dict_4["structure_factor_vec"][10:end], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=4" )
Plots.plot!(structure_factor_dict_1["wavenumber_vec"][10:end], structure_factor_dict_1["structure_factor_vec"][10:end], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=1" )
Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,22.5), ylims=(0,5), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_T_1_4_high_sampling_rate.png")


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.125, 0.25, 0.5, 1, 2, 4, 8]

b_range = range(0.4, 1.8, length=71)

my_plot = Plots.stephist()


bond_length_std_vec = zeros(length(temperatures))

for i in eachindex(temperatures)

    graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_"*string(temperatures[i])*"_quenched")

    bond_length_std_vec[i], bond_length_vec = NA.get_bond_length_std(graph_dict)

    my_plot = Plots.stephist!(bond_length_vec, bins=b_range, label = Latex.L"kT="*string(temperatures[i]), normalize=:probability)
end

my_plot = Plots.plot!(xlabel="bond length",
ylabel = "relative frequency")

Plots.savefig(path*"bond_length_216_vertices_T_0.125_8.png")


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.125, 0.25, 0.5, 1, 2, 4, 8]

b_range = range(0, 180, length=61)
my_plot = Plots.stephist()


bond_angle_std_vec = zeros(length(temperatures))

for i in eachindex(temperatures)

    graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_"*string(temperatures[i])*"_quenched")

    bond_angle_std_vec[i], bond_angle_vec = NA.get_bond_angle_std(graph_dict)


    my_plot = Plots.stephist!(bond_angle_vec./pi .* 180, bins=b_range, label = Latex.L"kT="*string(temperatures[i]), normalize=:probability)
end

my_plot = Plots.plot!(xlabel="bond angle",
ylabel = "relative frequency", xlims=(0, 180), legend = :topleft)

Plots.savefig(path*"bond_angle_216_vertices_T_0.125_8.png")




dict_path = raw"..\analysis_data\random_networks\216_vertices_T_"

temperatures = [0.125, 0.25, 0.5, 0.0625]

my_plot = Plots.plot()


bond_angle_std_vec = zeros(length(temperatures))

for i in eachindex(temperatures)

    structure_factor_dict = GU.load_h5_dict(dict_path*string(temperatures[i])*"_quenched_structure_factor_bartlett_isotrope.h5")

    my_plot = Plots.plot!(structure_factor_dict["wavenumber_vec"], structure_factor_dict["structure_factor_vec"], label = Latex.L"kT="*string(temperatures[i]) )
end

my_plot = Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,27.5), ylims=(0,5), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_bartlett_216_vertices_T_0.125_0.5.png")



temperatures = [1, 2, 4, 0.0625]

my_plot = Plots.plot()


bond_angle_std_vec = zeros(length(temperatures))

for i in eachindex(temperatures)

    structure_factor_dict = GU.load_h5_dict(dict_path*string(temperatures[i])*"_quenched_structure_factor_bartlett_isotrope.h5")

    my_plot = Plots.plot!(structure_factor_dict["wavenumber_vec"], structure_factor_dict["structure_factor_vec"], label = Latex.L"kT="*string(temperatures[i]) )
end

my_plot = Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,27.5), ylims=(0,5), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_bartlett_216_vertices_T_1_4.png")



temperatures = [4, 6, 8, 0.0625]

my_plot = Plots.plot()


bond_angle_std_vec = zeros(length(temperatures))

for i in eachindex(temperatures)

    structure_factor_dict = GU.load_h5_dict(dict_path*string(temperatures[i])*"_quenched_structure_factor_bartlett_isotrope.h5")

    my_plot = Plots.plot!(structure_factor_dict["wavenumber_vec"], structure_factor_dict["structure_factor_vec"], label = Latex.L"kT="*string(temperatures[i]) )
end

my_plot = Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,27.5), ylims=(0,5), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_bartlett_216_vertices_T_4_8.png")



dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

temperature = temperatures[1]

filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)
evolution_dict = GU.load_h5_dict(dict_path*filename*"_evolution.h5")

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"]) ) ./(graph_dict["nr_vertices"]*18), 
evolution_dict["total_energy_vec"]./graph_dict["nr_vertices"], label=Latex.L"kT_\mathrm{max}="*string(temperature) )

temperature = temperatures[3]

filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)
evolution_dict = GU.load_h5_dict(dict_path*filename*"_evolution.h5")

Plots.plot!(collect(1:length(evolution_dict["total_energy_vec"]) ) ./(graph_dict["nr_vertices"]*18), 
evolution_dict["total_energy_vec"] ./graph_dict["nr_vertices"], label=Latex.L"kT_\mathrm{max}="*string(temperature) )

Plots.plot!(xlabel="Monte Carlo step", ylabel="energy per vertex", xlims=(0, 20), right_margin = 3Plots.mm)

Plots.savefig(path*"total_energy_216_vertices_T_0.1_0.15_cool_0.1_per_mc_quenched.png")


temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

temperature = temperatures[1]

filename = "216_vertices_T_"*string(temperature)*"_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)
evolution_dict = GU.load_h5_dict(dict_path*filename*"_evolution.h5")

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"]) ) ./(graph_dict["nr_vertices"]*18), 
evolution_dict["total_energy_vec"]./graph_dict["nr_vertices"], label=Latex.L"kT_\mathrm{max}="*string(temperature) )

temperature = temperatures[3]

filename = "216_vertices_T_"*string(temperature)*"_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)
evolution_dict = GU.load_h5_dict(dict_path*filename*"_evolution.h5")

Plots.plot!(collect(1:length(evolution_dict["total_energy_vec"]) ) ./(graph_dict["nr_vertices"]*18), 
evolution_dict["total_energy_vec"]./graph_dict["nr_vertices"], label=Latex.L"kT_\mathrm{max}="*string(temperature) )

Plots.plot!(xlabel="Monte Carlo step", ylabel="energy per vertex", xlims=(0, 20), right_margin = 3Plots.mm)

Plots.savefig(path*"total_energy_216_vertices_T_0.1_0.15_quenched.png")


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

function heaviside(t)
    0.5 * (sign(t) + 1)
 end
 
temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

temperature = temperatures[1]

mc_step_vec = collect(0:0.1:20)
temperature_vec = heaviside.(- (mc_step_vec .- 2)) .* temperature

Plots.plot(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}="*string(temperature), ls= :dash)

temperature = temperatures[3]

temperature_vec = heaviside.(- (mc_step_vec .- 2)) .* temperature

Plots.plot!(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}="*string(temperature), ls= :dash )

Plots.plot!(xlabel="Monte Carlo step", ylabel=Latex.L"kT", xlims=(0, 20), right_margin = 3Plots.mm)

Plots.savefig(path*"temperature_T_0.1_0.15_quenched.png")




dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

function heaviside(t)
    0.5 * (sign(t) + 1)
 end
 
temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

temperature = temperatures[1]

mc_step_vec = collect(0:0.1:20)
temperature_vec = (heaviside.(- (mc_step_vec .- 2)) .* temperature 
.+ (0.3 .- 0.1 .* mc_step_vec ).* heaviside.(- (mc_step_vec .- 3)) .* heaviside.( (mc_step_vec .- 2)))


Plots.plot(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}="*string(temperature), ls= :dash)

temperature = temperatures[3]
temperature_vec = (heaviside.(- (mc_step_vec .- 2)) .* temperature 
.+ (0.35 .- 0.1 .* mc_step_vec) .* heaviside.(- (mc_step_vec .- 3.5)) .* heaviside.( (mc_step_vec .- 2)))

Plots.plot!(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}="*string(temperature), ls= :dash )

Plots.plot!(xlabel="Monte Carlo step", ylabel=Latex.L"kT", xlims=(0, 20), right_margin = 3Plots.mm)

Plots.savefig(path*"temperature_T_0.1_0.15_cool_0.1_per_mc_quenched.png")


x = collect(-1:0.01:1)

plot_1 = min.(1, exp.( .-  x ./ 0.1  ) )  
plot_2 = min.(1, exp.( .-  x ./ 0.4  ) )  

myplot = Plots.plot(x, plot_2, label = Latex.L"kT=0.4", linecolor=Plots.palette(:tab10)[2])
myplot = Plots.plot!(x, plot_1, label = Latex.L"kT=0.1", linecolor=Plots.palette(:tab10)[1])

Plots.plot!(grid=false, xlabel="energy difference",
ylabel = "acceptance probability")

Plots.savefig(path*"boltzmann_metropolis_temperatures_0.1_0.4.png")


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_1 = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_0.1_heated_for_0.5_steps_quenched"
)
graph_dict_4 = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_1.0_heated_for_0.5_steps_quenched"
)


bond_length_std_1, bond_length_vec_1 = NA.get_bond_length_std(graph_dict_1)
bond_length_std_4, bond_length_vec_4 = NA.get_bond_length_std(graph_dict_4)

b_range = range(0.4, 1.8, length=71)

my_plot = Plots.stephist(bond_length_vec_4, bins=b_range, label = Latex.L"kT=1.0", normalize=:probability, linecolor=Plots.palette(:tab10)[2])
my_plot = Plots.stephist!(bond_length_vec_1, bins=b_range, label = Latex.L"kT=0.1", normalize=:probability, linecolor=Plots.palette(:tab10)[1])
my_plot = Plots.plot!(xlabel="bond length",
ylabel = "relative frequency")

Plots.savefig(path*"bond_length_T_0.1_1.0.png")


bond_angle_std_1, bond_angle_vec_1 = NA.get_bond_angle_std(graph_dict_1)
bond_angle_std_4, bond_angle_vec_4 = NA.get_bond_angle_std(graph_dict_4)

b_range = range(0, 180, length=61)

my_plot = Plots.stephist(bond_angle_vec_4 ./pi .* 180, bins=b_range, label = Latex.L"kT=1.0", normalize=:probability, linecolor=Plots.palette(:tab10)[2])
my_plot = Plots.stephist!(bond_angle_vec_1 ./pi .* 180, bins=b_range, label = Latex.L"kT=0.1", normalize=:probability, linecolor=Plots.palette(:tab10)[1])
my_plot = Plots.plot!(xlabel="bond angle",
ylabel = "relative frequency")

Plots.savefig(path*"bond_angle_T_0.1_1.0.png")


dict_path_low_t = raw"..\analysis_data\random_networks\216_vertices_T_0.1_heated_for_0.5_steps_quenched_structure_factor_isotrope.h5"

dict_path_high_t = raw"..\analysis_data\random_networks\216_vertices_T_1.0_heated_for_0.5_steps_quenched_structure_factor_isotrope.h5"

structure_factor_dict_low_t = GU.load_h5_dict(dict_path_low_t)
structure_factor_dict_high_t = GU.load_h5_dict(dict_path_high_t)

my_plot = Plots.plot(structure_factor_dict_high_t["wavenumber_vec"], structure_factor_dict_high_t["structure_factor_vec"], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=0.1" )
my_plot = Plots.plot!(structure_factor_dict_low_t["wavenumber_vec"], structure_factor_dict_low_t["structure_factor_vec"], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=1.0" )
my_plot = Plots.plot!(xlabel="wavenumber",
ylabel = "structure factor", xlims=(0,32.5), ylims=(0,2.75), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_T_0.1_1.0.png")


x = collect(-2:0.01:2)

plot_1 = min.(1, exp.( .-  x ./ 0.1  ) )  
plot_2 = min.(1, exp.( .-  x ./ 1  ) )  

myplot = Plots.plot(x, plot_2, label = Latex.L"kT=1.0", linecolor=Plots.palette(:tab10)[2])
myplot = Plots.plot!(x, plot_1, label = Latex.L"kT=0.1", linecolor=Plots.palette(:tab10)[1])

Plots.plot!(grid=false, xlabel="energy difference",
ylabel = "acceptance probability")

Plots.savefig(path*"boltzmann_metropolis_temperatures_0.1_1.0.png")



dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_1 = NG.load_graph_from_h5_and_MGformat(dict_path*"1000_vertices_T_0.1_heated_for_0.5_steps_quenched"
)
graph_dict_4 = NG.load_graph_from_h5_and_MGformat(dict_path*"1000_vertices_T_1.0_heated_for_0.5_steps_quenched"
)


bond_length_std_1, bond_length_vec_1 = NA.get_bond_length_std(graph_dict_1)
bond_length_std_4, bond_length_vec_4 = NA.get_bond_length_std(graph_dict_4)

diamond_vec = ones(20)

b_range = range(0.705, 1.305, length=71)

my_plot = Plots.stephist(bond_length_vec_4, bins=b_range, label = Latex.L"kT=1.0", normalize=:probability, linecolor=Plots.palette(:tab10)[2])
my_plot = Plots.stephist!(bond_length_vec_1, bins=b_range, label = Latex.L"kT=0.1", normalize=:probability, linecolor=Plots.palette(:tab10)[1])
Plots.plot!(vcat([0.9], collect(0.995:0.001:1.005), [1.1]), vcat([0,0], ones(length(collect(0.995:0.001:1.005)) -2 ) .* 0.15, [0,0]), label = Latex.L"kT=0", linecolor=Plots.palette(:tab10)[3])
#Plots.vline!([1], label = Latex.L"kT=0", linecolor=Plots.palette(:tab10)[3])
my_plot = Plots.plot!(xlabel="bond length / "*Latex.L"d",
ylabel = "relative frequency", xlim=(0.7, 1.3), rightmargin=3Plots.mm)

Plots.savefig(path*"bond_length_1000_vertices_T_0.1_1.0.png")


bond_angle_std_1, bond_angle_vec_1 = NA.get_bond_angle_std(graph_dict_1)
bond_angle_std_4, bond_angle_vec_4 = NA.get_bond_angle_std(graph_dict_4)

b_range = range(40.6, 180.6, length=61)

my_plot = Plots.stephist(bond_angle_vec_4 ./pi .* 180, bins=b_range, label = Latex.L"kT=1.0", normalize=:probability, linecolor=Plots.palette(:tab10)[2])
my_plot = Plots.stephist!(bond_angle_vec_1 ./pi .* 180, bins=b_range, label = Latex.L"kT=0.1", normalize=:probability, linecolor=Plots.palette(:tab10)[1])
Plots.plot!(vcat([90], collect(109.5-1:0.01:109.5+1), [130]), vcat([0,0], ones(length(collect(109.5-1:0.01:109.5+1)) -2 ) .* 0.17, [0,0]), label = Latex.L"kT=0", linecolor=Plots.palette(:tab10)[3])
#Plots.vline!([109.5], label = Latex.L"kT=0", linecolor=Plots.palette(:tab10)[3])
my_plot = Plots.plot!(xlabel="bond angle / °",
ylabel = "relative frequency", xlim=(40, 180), rightmargin=5Plots.mm)

Plots.savefig(path*"bond_angle_1000_vertices_T_0.1_1.0.png")



dict_path_low_t = raw"..\analysis_data\random_networks\1000_vertices_T_0.1_heated_for_0.5_steps_quenched_structure_factor_bartlett_isotrope.h5"

dict_path_high_t = raw"..\analysis_data\random_networks\1000_vertices_T_1.0_heated_for_0.5_steps_quenched_structure_factor_bartlett_isotrope.h5"

structure_factor_dict_low_t = GU.load_h5_dict(dict_path_low_t)
structure_factor_dict_high_t = GU.load_h5_dict(dict_path_high_t)

my_plot = Plots.plot(structure_factor_dict_high_t["wavenumber_vec"], structure_factor_dict_high_t["structure_factor_vec"], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=1.0" )
my_plot = Plots.plot!(structure_factor_dict_low_t["wavenumber_vec"], structure_factor_dict_low_t["structure_factor_vec"], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0.1" )
my_plot = Plots.plot!(xlabel="wavenumber",
ylabel = "structure factor", xlims=(0,32.5), ylims=(0,2.75), xtick=pitick(0, 32, 1; mode=:latex))

Plots.savefig(path*"structure_factor_bartlett_1000_vertices_T_0.1_1.0.png")


dict_path_low_t = raw"..\analysis_data\random_networks\1000_vertices_T_0.1_heated_for_0.5_steps_quenched_structure_factor_bartlett_isotrope.h5"
dict_path_high_t = raw"..\analysis_data\random_networks\1000_vertices_T_1.0_heated_for_0.5_steps_quenched_structure_factor_bartlett_isotrope.h5"
dict_path_diamond = raw"..\analysis_data\random_networks\1000_vertices_perfect_diamond_structure_factor_bartlett_isotrope.h5"

structure_factor_dict_low_t = GU.load_h5_dict(dict_path_low_t)
structure_factor_dict_high_t = GU.load_h5_dict(dict_path_high_t)
structure_factor_dict_diamond = GU.load_h5_dict(dict_path_diamond)


my_plot = Plots.plot(structure_factor_dict_high_t["wavenumber_vec"], structure_factor_dict_high_t["structure_factor_vec"], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=1.0" )
my_plot = Plots.plot!(structure_factor_dict_low_t["wavenumber_vec"], structure_factor_dict_low_t["structure_factor_vec"], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0.1" )
my_plot = Plots.plot!(structure_factor_dict_diamond["wavenumber_vec"], structure_factor_dict_diamond["structure_factor_vec"], linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0", alpha=0.5 )
my_plot = Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}",
ylabel = "structure factor", xlims=(0,32.5), ylims=(0,4), size = (500, 600),  xtick=pitick(0, 32, 1; mode=:latex), bottommargin = 0Plots.mm)

Plots.savefig(path*"structure_factor_bartlett_1000_vertices_T_0.1_1.0_stretched.png")


dict_path_low_t = raw"..\analysis_data\random_networks\1000_vertices_T_0.1_heated_for_0.5_steps_quenched_structure_factor_bartlett_isotrope.h5"

dict_path_high_t = raw"..\analysis_data\random_networks\1000_vertices_T_1.0_heated_for_0.5_steps_quenched_structure_factor_bartlett_isotrope.h5"

structure_factor_dict_low_t = GU.load_h5_dict(dict_path_low_t)
structure_factor_dict_high_t = GU.load_h5_dict(dict_path_high_t)

first_index = 11

my_plot = Plots.plot(structure_factor_dict_high_t["wavenumber_vec"][first_index:end], structure_factor_dict_high_t["structure_factor_vec"][first_index:end], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=1.0" )
my_plot = Plots.plot!(structure_factor_dict_low_t["wavenumber_vec"][first_index:end], structure_factor_dict_low_t["structure_factor_vec"][first_index:end], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0.1" )
my_plot = Plots.plot!(structure_factor_dict_diamond["wavenumber_vec"][first_index:end], structure_factor_dict_diamond["structure_factor_vec"][first_index:end], linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0", alpha=0.5 )
my_plot = Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}",
ylabel = "structure factor", xlims=(0,14), ylims=(0,3.75), xtick=pitick(0, 14, 1; mode=:latex), legend = :topright)

Plots.savefig(path*"structure_factor_bartlett_small_k_1000_vertices_T_0.1_1.0.png")


x = collect(0:0.005:1)

plot_1 = (1 .- 2 .* x.^2 + 2 .* x.^4) .* x

myplot = Plots.plot(x, plot_1)

Plots.plot!(grid=false, xlabel=Latex.L"q", legend = false,
ylabel = Latex.L"P(q)", xlims = (0, 1), ylims = (0, 1),
right_margin = 4Plots.mm,)

Plots.savefig(path*"polynomial_multiple_scattering.png")



l_max = 12

some_vertex = 10

y_vec = collect(0:l_max)

dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_low_t = NG.load_graph_from_h5_and_MGformat(dict_path*
    "1000_vertices_T_0.1_heated_for_0.5_steps_quenched")

graph_dict_high_t = NG.load_graph_from_h5_and_MGformat(dict_path*
    "1000_vertices_T_1.0_heated_for_0.5_steps_quenched")

graph_dict_diamond = NG.load_graph_from_h5_and_MGformat(dict_path*
    "216_vertices_T_0.1_heated_for_0.01_steps_quenched")

q_l_total_network_mean_dict = NA.get_q_l_total_network_mean_dict(graph_dict_high_t, l_max)
q_l_total_network_mean_vec = NA.convert_q_l_dict_to_vec(q_l_total_network_mean_dict, l_max)

Plots.scatter(y_vec, Measurements.value.(q_l_total_network_mean_vec),
        yerr=Measurements.uncertainty.(q_l_total_network_mean_vec), linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=1.0", markerstrokecolor=Plots.palette(:tab10)[2] , markercolor=Plots.palette(:tab10)[2], markersize = 6)


q_l_total_network_mean_dict = NA.get_q_l_total_network_mean_dict(graph_dict_low_t, l_max)
q_l_total_network_mean_vec = NA.convert_q_l_dict_to_vec(q_l_total_network_mean_dict, l_max)

Plots.scatter!(y_vec, Measurements.value.(q_l_total_network_mean_vec),
    yerr=Measurements.uncertainty.(q_l_total_network_mean_vec), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0.1", markerstrokecolor=Plots.palette(:tab10)[1], markercolor=Plots.palette(:tab10)[1], markersize = 6)


q_l_total_network_mean_dict = NA.get_q_l_averaged_single_vertex_dict(graph_dict_diamond, some_vertex, l_max)
q_l_total_network_mean_vec = NA.convert_q_l_dict_to_vec(q_l_total_network_mean_dict, l_max)

Plots.scatter!(y_vec, Measurements.value.(q_l_total_network_mean_vec), linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0", markerstrokecolor=Plots.palette(:tab10)[3], markercolor=Plots.palette(:tab10)[3], markersize = 6)

Plots.scatter!([3,7,9,11], zeros(4), linecolor=:grey, label = "BCC", markerstrokecolor=:grey, markercolor=:grey, markersize = 6)
    
Plots.plot!(xlabel=Latex.L"l", 
    ylabel=Latex.L"\overline{q}_l", ylims=(0,1.3), xticks=0:2:12, yticks=0:0.25:1)


Plots.savefig(path*"q_l_total_network_mean_1000_vertices_T_0.1_1.0_diamond.png")



myplot = Plots.plot(collect(0:0.1:0.5), ones(length(collect(0:0.1:0.5))), linecolor=Plots.palette(:tab10)[2] )
myplot = Plots.plot!(collect(0.5:0.1:40), zeros(length(collect(0.5:0.1:40))), linecolor=Plots.palette(:tab10)[2], ls=:dot  )
myplot = Plots.plot!(collect(0:0.1:0.5), ones(length(collect(0:0.1:0.5))) .* 0.1, linecolor=Plots.palette(:tab10)[1]  )
myplot = Plots.plot!(collect(0.5:0.1:40) .+ 0.05, zeros(length(collect(0.5:0.1:40))), linecolor=Plots.palette(:tab10)[1], ls=:dot  )

Plots.plot!(grid=false, xlabel="step", legend = false,
ylabel = Latex.L"kT", xlims=(0, 2.5), size=(300,300), xticks=[0,1,2], yticks=[0,1])

Plots.savefig(path*"temperature_profile_T_0.1_1.0.png")


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

filename = "1000_vertices_T_1.0_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)
evolution_dict = GU.load_h5_dict(dict_path*filename*"_evolution.h5")

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"]) ) ./(graph_dict["nr_vertices"]*18), 
evolution_dict["total_energy_vec"]./graph_dict["nr_vertices"], label=Latex.L"kT_\mathrm{max}=1.0", linecolor=Plots.palette(:tab10)[2])

filename = "1000_vertices_T_0.1_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)
evolution_dict = GU.load_h5_dict(dict_path*filename*"_evolution.h5")

Plots.plot!(collect(1:length(evolution_dict["total_energy_vec"]) ) ./(graph_dict["nr_vertices"]*18), 
evolution_dict["total_energy_vec"] ./graph_dict["nr_vertices"], label=Latex.L"kT_\mathrm{max}=0.1", linecolor=Plots.palette(:tab10)[1] )

Plots.plot!(xlabel="Monte Carlo step", ylabel="energy per vertex", xlims=(0, 20), right_margin = 3Plots.mm)

Plots.savefig(path*"1000_vertices_T_0.1_1.0_heated_for_0.5_steps_quenched.png")


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

function heaviside(t)
    0.5 * (sign(t) + 1)
 end

temperature = 1.0

mc_step_vec = collect(0:0.01:20)
temperature_vec = heaviside.(- (mc_step_vec .- 0.5)) .* temperature

Plots.plot(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}=1.0", ls= :dash, linecolor=Plots.palette(:tab10)[2])

temperature = 0.1

temperature_vec = heaviside.(- (mc_step_vec .- 0.5)) .* temperature

Plots.plot!(mc_step_vec, 
temperature_vec, label=Latex.L"kT_\mathrm{max}=0.1", ls= :dash, linecolor=Plots.palette(:tab10)[1] )

Plots.plot!(xlabel="Monte Carlo step", ylabel=Latex.L"kT", xlims=(0, 20), right_margin = 3Plots.mm)

Plots.savefig(path*"temperature_T_0.1_1.0_quenched.png")



dict_path = raw"..\analysis_data\random_networks\\"
filename = "216_vertices_T_0.2_heated_for_0.5_steps_quenched"

structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

save_path = raw"..\plots\random_networks\\"*filename*"_y_z"

NA.plot_structure_factor_heatmap(structure_factor_dict,
    save_path;
    title="Structure factor",
    save_plot = true,
    clims = (0, 4 ),
    x_y_lims = nothing,
    wavevector_component_to_fix = 1,
    wavevector_value_fixed = 0)


save_path = raw"..\plots\random_networks\\"*filename*"_x_z"

NA.plot_structure_factor_heatmap(structure_factor_dict,
    save_path;
    title="Structure factor",
    save_plot = true,
    clims = (0, 4 ),
    x_y_lims = nothing,
    wavevector_component_to_fix = 2,
    wavevector_value_fixed = 0)

save_path = raw"..\plots\random_networks\\"*filename*"_x_y"

NA.plot_structure_factor_heatmap(structure_factor_dict,
    save_path;
    title="Structure factor",
    save_plot = true,
    clims = (0, 4 ),
    x_y_lims = nothing,
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)

save_path = raw"..\plots\random_networks\\"*filename*"_x_z_offset"

NA.plot_structure_factor_heatmap(structure_factor_dict,
        save_path;
        title="Structure factor",
        save_plot = true,
        clims = (0, 4 ),
        x_y_lims = nothing,
        wavevector_component_to_fix = 3,
        wavevector_value_fixed = 2*pi)


dict_path = raw"..\analysis_data\random_networks\\"
filename = "216_vertices_T_0.2_heated_for_0.5_steps_quenched"

structure_factor_angle_averaged_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_angle_averaged.h5")

Plots.plot(structure_factor_angle_averaged_dict["wavenumber_vec"], 
                    Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]) , 
                    ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict["structure_factor_vec"]))

Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "structure factor", xtick=pitick(0, 56, 1; mode=:latex), legend = false, xlims=(0, 53), ylims=(0, 2))

Plots.savefig(raw"..\plots\random_networks\\"*filename*"_structure_factor_angle_averaged.png")


dict_path = raw"..\analysis_data\random_networks\\"

filename = "1000_vertices_T_1.0_heated_for_0.5_steps_quenched"
structure_factor_angle_averaged_dict_high_t = GU.load_h5_dict(dict_path*filename*"_structure_factor_angle_averaged.h5")

filename = "1000_vertices_T_0.1_heated_for_0.5_steps_quenched"
structure_factor_angle_averaged_dict_low_t = GU.load_h5_dict(dict_path*filename*"_structure_factor_angle_averaged.h5")

filename = "1000_vertices_perfect_diamond"
structure_factor_angle_averaged_dict_diamond = GU.load_h5_dict(dict_path*filename*"_structure_factor_angle_averaged.h5")


Plots.plot(structure_factor_angle_averaged_dict_diamond["wavenumber_vec"], 
                    Measurements.value.(structure_factor_angle_averaged_dict_diamond["structure_factor_vec"]) , 
                    ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict_diamond["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[3])

Plots.plot!(structure_factor_angle_averaged_dict_low_t["wavenumber_vec"], 
                    Measurements.value.(structure_factor_angle_averaged_dict_low_t["structure_factor_vec"]) , 
                    ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict_low_t["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0.1", fillcolor = Plots.palette(:tab10)[1])

Plots.plot!(structure_factor_angle_averaged_dict_high_t["wavenumber_vec"], 
                    Measurements.value.(structure_factor_angle_averaged_dict_high_t["structure_factor_vec"]) , 
                    ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict_high_t["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=1.0", fillcolor = Plots.palette(:tab10)[2])

Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "structure factor", xlims=(0,32.5), ylims=(0,4), size = (500, 600),  xtick=pitick(0, 32, 1; mode=:latex), leftmargin = 0Plots.mm)

Plots.savefig(raw"..\plots\random_networks\structure_factor_angle_averaged_1000_vertices_T_0.1_1.0_stretched.png")

Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "structure factor", xlims=(0,32.5), ylims=(0,12), size = (500, 800),  xtick=pitick(0, 32, 1; mode=:latex), leftmargin = 4Plots.mm)

Plots.savefig(raw"..\plots\random_networks\structure_factor_angle_averaged_1000_vertices_T_0.1_1.0_very_stretched.png")


filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"
data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename
save_path = raw"..\plots\random_networks\\"*filename

plot_dict = GU.load_h5_dict(data_path*"_autocovariance_fct_direction.h5")

NA.plot_autocovariance_fct_heatmap(plot_dict,
    save_path;
    save_plot = true,
    clims = nothing,
    x_y_lims = nothing,
    sampling_vector_component_to_fix = 3,
    sampling_vector_value_fixed = 0)


plot_dict = GU.load_h5_dict(data_path*"_spectral_density_array.h5")

NA.plot_spectral_density_heatmap(plot_dict,
    save_path;
    save_plot = true,
    clims = (0,0.1),
    x_y_lims = nothing,
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)


filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"
data_path = raw"..\analysis_data\random_networks\\"*filename
save_path = raw"..\plots\random_networks\\"*filename

spectral_density_angle_averaged_dict = GU.load_h5_dict(data_path*"_spectral_density_angle_averaged.h5")

Plots.plot(spectral_density_angle_averaged_dict["wavenumber_vec"], 
                    Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"]) , 
                    ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"]))

Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", legend = false, xlims=(0,30), ylims=(0,30), xtick=pitick(0, 30, 1; mode=:latex))

Plots.savefig(raw"..\plots\random_networks\\"*filename*"_spectral_density_angle_averaged.png")



filename = "216_vertices_T_2.0_heated_for_0.5_steps_quenched"
structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\run_4\\"
structure_dict = GU.load_h5_dict(structure_dict_path*filename*"_structure.h5")

NA.plot_binary_structure(structure_dict["data_binary"])


filename = "216_vertices_T_0.1_heated_for_0.05_steps_quenched"
data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename
save_path = raw"..\plots\random_networks\\"*filename

spectral_density_dict = GU.load_h5_dict(data_path*"_spectral_density_array.h5")

spectral_density_angle_averaged_dict = GU.load_h5_dict(data_path*"_spectral_density_angle_averaged.h5")

Plots.plot(spectral_density_angle_averaged_dict["wavenumber_vec"], 
                    Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"]) , 
                    ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"]))

Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", legend = false, xlims=(0,15), ylims=(0,1000), xtick=pitick(0, 15, 1; mode=:latex))

Plots.savefig(raw"..\plots\random_networks\\"*filename*"_spectral_density_angle_averaged.png")

Plots.plot(spectral_density_angle_averaged_dict["unfiltered_wavenumber_vec"], 
                    Measurements.value.(spectral_density_angle_averaged_dict["unfiltered_spectral_density_vec"]))

Plots.plot!(xlabel="wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", legend = false, xlims=(0,15), ylims=(0,1000), xtick=pitick(0, 15, 1; mode=:latex))

Plots.savefig(raw"..\plots\random_networks\\"*filename*"_spectral_density_angle_averaged_unfiltered.png")


graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\run_1\\"

all_filenames = readdir(graph_dict_path)
filenames = filter(filename -> endswith(filename, "_evolution.h5"), all_filenames)
final_energy_vec = Float64[]

for filename in filenames
    println(filename)

    evolution_dict = GU.load_h5_dict(graph_dict_path*filename)

    push!(final_energy_vec, evolution_dict["total_energy_vec"][end])
end

filenames_sorted = filenames[sortperm(final_energy_vec)   ]
sort!(final_energy_vec)


evolution_dict_high = GU.load_h5_dict(graph_dict_path*filenames_sorted[end])
evolution_dict_middle = GU.load_h5_dict(graph_dict_path*filenames_sorted[13])
evolution_dict_low = GU.load_h5_dict(graph_dict_path*filenames_sorted[1])

Plots.plot(collect(1:length(evolution_dict_high["total_energy_vec"]))./(216*18), evolution_dict_high["total_energy_vec"] ./ 216, ls = :dot, label = "high final energy")
Plots.plot!(collect(1:length(evolution_dict_middle["total_energy_vec"]))./(216*18), evolution_dict_middle["total_energy_vec"] ./ 216, ls = :dot, label = "middle final energy")
Plots.plot!(collect(1:length(evolution_dict_low["total_energy_vec"]))./(216*18), evolution_dict_low["total_energy_vec"] ./ 216, ls = :dot, label = "low final energy", yaxis=:log10)

Plots.plot!(xlabel = "Monte Carlo step", ylabel = "energy per vertex", size=(600, 400), legend=:bottomright)
Plots.savefig(raw"..\plots\random_networks\\216_vertices_anneal_quench_total_energy.png")


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


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_globally_relaxed\run_2\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_2\\"

order_metrics_dict = GU.load_h5_dict(analysis_data_path*"all_order_metrics.h5")

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec", "cluster_metric_vec"]

y_labels = ["bond length std", "bond angle std", "dihedral angle std", "anisotropy in structure f.", "anisotropy in spectral d.", "cluster metric"]

for i in eachindex(order_metrics_names)
    Plots.scatter(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216, order_metrics_dict[order_metrics_names[i]][1:end-2], xlabel = "Keating energy per vertex", ylabel = y_labels[i], legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216)),
    ylims = (minimum(order_metrics_dict[order_metrics_names[i]][1:end-2]), maximum(order_metrics_dict[order_metrics_names[i]][1:end-2])))
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
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216)),
ylims = (minimum(order_metrics_dict["bond_length_std_vec"][1:end-2]), maximum(order_metrics_dict["bond_length_std_vec"][1:end-2])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_length_std.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_bond_angle_std_vec, xlabel = "Keating energy per vertex", ylabel = "bond angle std", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216)),
ylims = (minimum(order_metrics_dict["bond_angle_std_vec"][1:end-2]), maximum(order_metrics_dict["bond_angle_std_vec"][1:end-2])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_angle_std.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_dihedral_angle_std_vec, xlabel = "Keating energy per vertex", ylabel = "dihedral angle std", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216)),
ylims = (minimum(order_metrics_dict["dihedral_angle_std_vec"][1:end-2]), maximum(order_metrics_dict["dihedral_angle_std_vec"][1:end-2])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_dihedral_angle_std.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_anisotropy_metric_from_structure_factor_vec, xlabel = "Keating energy per vertex", ylabel = "anisotropy in structure f.", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216)),
ylims = (minimum(order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][1:end-2]), maximum(order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][1:end-2])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_anisotropy_in_structure_factor.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_anisotropy_metric_from_spectral_density_vec, xlabel = "Keating energy per vertex", ylabel = "anisotropy in spectral d.", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216)),
ylims = (minimum(order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][1:end-2]), maximum(order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][1:end-2])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_anisotropy_in_spectral_density.png")

    Plots.scatter(filtered_total_keating_energy_vec./216, filtered_cluster_metric_vec, xlabel = "Keating energy per vertex", ylabel = "cluster metric", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-2]./216)),
ylims = (minimum(order_metrics_dict["cluster_metric_vec"][1:end-2]), maximum(order_metrics_dict["cluster_metric_vec"][1:end-2])),
title = title_vec[i], color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_cluster_metric.png")
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



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\presentations\material\\"

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

i=9

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


Plots.scatter(filtered_total_keating_energy_vec./216, filtered_bond_length_std_vec, xlabel = "Keating energy per vertex", ylabel = Latex.L"\sigma_\mathrm{bond} _\mathrm{length}", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["bond_length_std_vec"][1:end-10]), maximum(order_metrics_dict["bond_length_std_vec"][1:end-10])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_length_std.png")

Plots.scatter(filtered_total_keating_energy_vec./216, filtered_bond_angle_std_vec, xlabel = "Keating energy per vertex", ylabel = Latex.L"\sigma_\mathrm{bond} _\mathrm{angle}", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["bond_angle_std_vec"][1:end-10]), maximum(order_metrics_dict["bond_angle_std_vec"][1:end-10])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_angle_std.png")

Plots.scatter(filtered_total_keating_energy_vec./216, filtered_cluster_metric_vec, xlabel = "Keating energy per vertex", ylabel = "cluster metric", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["cluster_metric_vec"][1:end-10]), maximum(order_metrics_dict["cluster_metric_vec"][1:end-10])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_cluster_metric.png")



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


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\presentations\material\\"

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

i=9

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

Plots.scatter(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216, order_metrics_dict[order_metrics_names[1]][1:end-10], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_bond_length_std_vec, xlabel = "Keating energy per vertex", ylabel = Latex.L"\sigma_\mathrm{bond} _\mathrm{length}", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["bond_length_std_vec"][1:end-10]), maximum(order_metrics_dict["bond_length_std_vec"][1:end-10])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_length_std.png")

Plots.scatter(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216, order_metrics_dict[order_metrics_names[2]][1:end-10], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_bond_angle_std_vec, xlabel = "Keating energy per vertex", ylabel = Latex.L"\sigma_\mathrm{bond} _\mathrm{angle}", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["bond_angle_std_vec"][1:end-10]), maximum(order_metrics_dict["bond_angle_std_vec"][1:end-10])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_angle_std.png")

Plots.scatter(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216, order_metrics_dict[order_metrics_names[4]][1:end-10], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_anisotropy_metric_from_structure_factor_vec, xlabel = "Keating energy per vertex", ylabel = "anisotropy", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][1:end-10]), maximum(order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][1:end-10])), color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_anisotropy_in_structure_factor.png")

Plots.scatter(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216, order_metrics_dict[order_metrics_names[5]][1:end-10], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_anisotropy_metric_from_spectral_density_vec, xlabel = "Keating energy per vertex", ylabel = "anisotropy", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][1:end-10]), maximum(order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][1:end-10])),color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_anisotropy_in_spectral_density.png")

Plots.scatter(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216, order_metrics_dict[order_metrics_names[6]][1:end-10], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_cluster_metric_vec, xlabel = "Keating energy per vertex", ylabel = "cluster metric", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216), maximum(order_metrics_dict["total_keating_energy_vec"][1:end-10]./216)),
ylims = (minimum(order_metrics_dict["cluster_metric_vec"][1:end-10]), maximum(order_metrics_dict["cluster_metric_vec"][1:end-10])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_cluster_metric.png")



function heat_cool_temperature_vec(x, max_temp, heat_cool_rate)

    if x < max_temp/heat_cool_rate 
        return heat_cool_rate * x
    elseif x < 2*max_temp/heat_cool_rate 
        return 2* max_temp - heat_cool_rate*x
    else
        return 0
    end
end

x_vec = collect(0:0.01:8)
y_vec = heat_cool_temperature_vec.(x_vec, 0.1, 0.1)

min_temp = 0.1
max_temp = 0.5
temperatures = [0.1, 0.125, 0.15,  0.2, 0.25, 0.3, 0.4, 0.5]
normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

Plots.plot(x_vec, heat_cool_temperature_vec.(x_vec, 0.25, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[5])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.2, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[4])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.15, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[3])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.1, 0.1), xlabel="Monte Carlo step", ylabel=Latex.L"kT", alpha=1.0, color=mapped_colors[1])
Plots.plot!(legend = false)

Plots.savefig(path*"heat_cool_0.1_temperature_profile.png")



analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285\\"

order_metrics_dict = Dict()

# loop through folders and append all order metrics to the order_metrics_dict
for i in 1:5
    
    current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285\run_"*string(i)*"\\"

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

i=9

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

Plots.scatter(order_metrics_dict["total_keating_energy_vec"]./216, order_metrics_dict[order_metrics_names[1]], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_bond_length_std_vec, xlabel = "Keating energy per vertex", ylabel = Latex.L"\sigma_\mathrm{bond} _\mathrm{length}", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"]./216), maximum(order_metrics_dict["total_keating_energy_vec"]./216)),
ylims = (minimum(order_metrics_dict["bond_length_std_vec"]), maximum(order_metrics_dict["bond_length_std_vec"])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_length_std.png")

Plots.scatter(order_metrics_dict["total_keating_energy_vec"]./216, order_metrics_dict[order_metrics_names[2]], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_bond_angle_std_vec, xlabel = "Keating energy per vertex", ylabel = Latex.L"\sigma_\mathrm{bond} _\mathrm{angle}", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"]./216), maximum(order_metrics_dict["total_keating_energy_vec"]./216)),
ylims = (minimum(order_metrics_dict["bond_angle_std_vec"]), maximum(order_metrics_dict["bond_angle_std_vec"])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_bond_angle_std.png")

Plots.scatter(order_metrics_dict["total_keating_energy_vec"]./216, order_metrics_dict[order_metrics_names[4]], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_anisotropy_metric_from_structure_factor_vec, xlabel = "Keating energy per vertex", ylabel = "anisotropy", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"]./216), maximum(order_metrics_dict["total_keating_energy_vec"]./216)),
ylims = (minimum(order_metrics_dict["anisotropy_metric_from_structure_factor_vec"]), maximum(order_metrics_dict["anisotropy_metric_from_structure_factor_vec"])), color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_anisotropy_in_structure_factor.png")

Plots.scatter(order_metrics_dict["total_keating_energy_vec"]./216, order_metrics_dict[order_metrics_names[5]], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_anisotropy_metric_from_spectral_density_vec, xlabel = "Keating energy per vertex", ylabel = "anisotropy", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"]./216), maximum(order_metrics_dict["total_keating_energy_vec"]./216)),
ylims = (minimum(order_metrics_dict["anisotropy_metric_from_spectral_density_vec"]), maximum(order_metrics_dict["anisotropy_metric_from_spectral_density_vec"])),color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_anisotropy_in_spectral_density.png")

Plots.scatter(order_metrics_dict["total_keating_energy_vec"]./216, order_metrics_dict[order_metrics_names[6]], xlabel = "Keating energy per vertex", alpha=0.1, color= :black)
Plots.scatter!(filtered_total_keating_energy_vec./216, filtered_cluster_metric_vec, xlabel = "Keating energy per vertex", ylabel = "cluster metric", legend = false,
    xlims = (minimum(order_metrics_dict["total_keating_energy_vec"]./216), maximum(order_metrics_dict["total_keating_energy_vec"]./216)),
ylims = (minimum(order_metrics_dict["cluster_metric_vec"]), maximum(order_metrics_dict["cluster_metric_vec"])),
color = mapped_colors)
    Plots.savefig(path*filename_vec[i]*"_cluster_metric.png")


function heat_cool_temperature_vec(x, max_temp, heat_cool_rate)

    if x < max_temp/heat_cool_rate 
        return heat_cool_rate * x
    elseif x < 2*max_temp/heat_cool_rate 
        return 2* max_temp - heat_cool_rate*x
    else
        return 0
    end
end


function heat_cool_temperature_vec(x, max_temp, heat_cool_rate)

    if x < max_temp/heat_cool_rate 
        return heat_cool_rate * x
    elseif x < 2*max_temp/heat_cool_rate 
        return 2* max_temp - heat_cool_rate*x
    else
        return 0
    end
end

x_vec = collect(0:0.01:10)
y_vec = heat_cool_temperature_vec.(x_vec, 0.1, 0.1)

min_temp = 0.1
max_temp = 0.5
temperatures = [0.1, 0.125, 0.15,  0.2, 0.25, 0.3, 0.4, 0.5]
normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]


Plots.plot(x_vec, heat_cool_temperature_vec.(x_vec, 0.4, 0.1),  alpha=1.0, color=mapped_colors[7])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.25, 0.1),  alpha=1.0, color=mapped_colors[5])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.2, 0.1),  alpha=1.0, color=mapped_colors[4])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.15, 0.1), alpha=1.0, color=mapped_colors[3])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.1, 0.1),  alpha=1.0, color=mapped_colors[1])
Plots.plot!(legend = false, xlabel="Monte Carlo move per bond chain", ylabel=Latex.L"kT", right_margin = 2Plots.mm)

Plots.savefig(path*"heat_cool_0.1_temperature_profile_2.png")



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\1000_vertices_bond_bending_0.285\run_1\\"

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.285\run_1\\"

filename = "1000_vertices_T_0.22_heat_cool_0.1_per_mc_quenched_correlation_functions.h5"

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[1:6]

Plots.plot()

for temperature in temperature_vec
    corr_fct_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_correlation_functions.h5")
    Plots.plot!(corr_fct_dict["vertex_distance_vec"], corr_fct_dict["pair_correlation_fct_vec"], label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "distance", ylabel = "pair correlation function")

Plots.savefig(path*"pair_correlation_function_low_t.png")

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[7:end]

Plots.plot()

for temperature in temperature_vec
    corr_fct_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_correlation_functions.h5")
    Plots.plot!(corr_fct_dict["vertex_distance_vec"], corr_fct_dict["pair_correlation_fct_vec"], label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "distance", ylabel = "pair correlation function")

Plots.savefig(path*"pair_correlation_function_high_t.png")



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\others\\"

network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_2\216_vertices_T_0.15_heat_cool_0.1_per_mc_quenched_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)

pore_pixel_radius_array = NA.get_pore_size_distribution(structure_dict_network)

pore_pixel_radius_vec = vec(pore_pixel_radius_array)
pore_pixel_radius_filtered_vec = pore_pixel_radius_vec[pore_pixel_radius_vec .> 0.0]
Plots.histogram(pore_pixel_radius_filtered_vec)
Plots.savefig(path*"pore_size_distribution.png")

sampling_nr = 100000
pore_pixel_radius_sampled_1 = StatsBase.sample(pore_pixel_radius_vec, sampling_nr, replace=false)
pore_pixel_radius_filtered_vec_1 = pore_pixel_radius_sampled_1[pore_pixel_radius_sampled_1 .> 0.0]
Plots.histogram(pore_pixel_radius_filtered_vec_1)
Plots.savefig(path*"pore_size_distribution_100000_samples.png")

sampling_nr = 50000
pore_pixel_radius_sampled_2 = StatsBase.sample(pore_pixel_radius_vec, sampling_nr, replace=false)
pore_pixel_radius_filtered_vec_2 = pore_pixel_radius_sampled_2[pore_pixel_radius_sampled_2 .> 0.0]
Plots.histogram(pore_pixel_radius_filtered_vec_2)
Plots.savefig(path*"pore_size_distribution_50000_samples.png")

sampling_nr = 10000
pore_pixel_radius_sampled_3 = StatsBase.sample(pore_pixel_radius_vec, sampling_nr, replace=false)
pore_pixel_radius_filtered_vec_3 = pore_pixel_radius_sampled_3[pore_pixel_radius_sampled_3 .> 0.0]
Plots.histogram(pore_pixel_radius_filtered_vec_3)
Plots.savefig(path*"pore_size_distribution_10000_samples.png")

sampling_nr = 5000
pore_pixel_radius_sampled_4 = StatsBase.sample(pore_pixel_radius_vec, sampling_nr, replace=false)
pore_pixel_radius_filtered_vec_4 = pore_pixel_radius_sampled_4[pore_pixel_radius_sampled_4 .> 0.0]
Plots.histogram(pore_pixel_radius_filtered_vec_4)
Plots.savefig(path*"pore_size_distribution_5000_samples.png")



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\1000_vertices_bond_bending_0.285\run_1\\"

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bb_0.285\run_1\\"

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[1:6]

Plots.plot()

for temperature in temperature_vec
    spectral_density_angle_averaged_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_spectral_density_angle_averaged.h5")
    Plots.plot!(spectral_density_angle_averaged_dict["wavenumber_vec"], Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"]), 
    #ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"]),
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", xlim = (0, 4*pi), ylim=(0,150) )

Plots.savefig(path*"spectral_density_angle_averaged_low_t.png")

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[7:end]

Plots.plot()

for temperature in temperature_vec
    spectral_density_angle_averaged_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_spectral_density_angle_averaged.h5")
    Plots.plot!(spectral_density_angle_averaged_dict["wavenumber_vec"], Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"]), 
    #ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"]),
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", xlim = (0, 4*pi), ylim=(0,150) )

Plots.savefig(path*"spectral_density_angle_averaged_high_t.png")



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\\"

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\\"

Plots.plot()

spectral_density_angle_averaged_dict = GU.load_h5_dict(load_path * "pachy_blue_spectral_density_angle_averaged.h5")

Plots.plot!(spectral_density_angle_averaged_dict["wavenumber_vec"], Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"])./3, 
    ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"])./3,
    label = "blue" )

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bb_0.285\run_1\\"

temperature_vec = [0.1, 0.15, 0.17]


for temperature in temperature_vec
    spectral_density_angle_averaged_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_spectral_density_angle_averaged.h5")
    Plots.plot!(spectral_density_angle_averaged_dict["wavenumber_vec"], Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"]), 
    ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"]),
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density / a. u.", xlim = (0, 3*pi), ylim=(0,300) )

Plots.savefig(path*"spectral_density_angle_averaged_weevil_random_networks.png")



load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\\"

Plots.plot()

spectral_density_angle_averaged_dict = GU.load_h5_dict(load_path * "pachy_blue_spectral_density_angle_averaged.h5")

Plots.plot!(spectral_density_angle_averaged_dict["wavenumber_vec"], Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"])./3, 
    ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"])./3,
    label = "blue" )

spectral_density_angle_averaged_dict = GU.load_h5_dict(load_path * "pachy_red_spectral_density_angle_averaged.h5")

Plots.plot!(spectral_density_angle_averaged_dict["wavenumber_vec"], Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"])./1.7 ./3, 
    ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"])./3,
    label = "red" )


Plots.plot!(xlabel = "wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density / a. u.", xlim = (0, 3*pi), ylim=(0,600) )

Plots.savefig(path*"spectral_density_angle_averaged_weevil.png")



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
        Plots.scatter!( yscale = :log10)
        fancylogscale!()
    else
        Plots.scatter!(yscale = :lin)
        
    end

    Plots.scatter!(legend = false, ylabel = order_metrics_labels[i], xlabel = "Keating energy per vertex")
    Plots.savefig(plots_save_path*order_metrics_names[i]*".png")
end



plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\heat_cool_bond_bending_0.21\\"

diamonds_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

diamond_order_metrics_dicts = [GU.load_h5_dict(diamonds_analysis_data_path*"216_vertices_perfect_diamond_small_scale_order_metrics.h5"),
GU.load_h5_dict(diamonds_analysis_data_path*"512_vertices_perfect_diamond_small_scale_order_metrics.h5"),
GU.load_h5_dict(diamonds_analysis_data_path*"1000_vertices_perfect_diamond_small_scale_order_metrics.h5")]


analysis_data_paths = [
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.21_heat_cool\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.21\\"]

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
        Plots.scatter!( yscale = :log10)
        fancylogscale!()
    else
        Plots.scatter!(yscale = :lin)
        
    end

    Plots.scatter!(legend = false, ylabel = order_metrics_labels[i], xlabel = "Keating energy per vertex")
    Plots.savefig(plots_save_path*order_metrics_names[i]*".png")
end





path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\1000_vertices_bond_bending_0.21\run_1\\"

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.21\run_1\\"

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


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\1000_vertices_bond_bending_0.285\run_1\\"

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bb_0.285\run_1\\"

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[1:6]

Plots.plot()

for temperature in temperature_vec
    spectral_density_angle_averaged_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_spectral_density_angle_averaged.h5")
    Plots.plot!(spectral_density_angle_averaged_dict["wavenumber_vec"], Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"]), 
    #ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"]),
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", xlim = (0, 4*pi), ylim=(0,150) )

Plots.savefig(path*"spectral_density_angle_averaged_low_t.png")

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[7:end]

Plots.plot()

for temperature in temperature_vec
    spectral_density_angle_averaged_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_spectral_density_angle_averaged.h5")
    Plots.plot!(spectral_density_angle_averaged_dict["wavenumber_vec"], Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"]), 
    #ribbon =  Measurements.uncertainty.(spectral_density_angle_averaged_dict["spectral_density_vec"]),
    label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "wavenumber / "*Latex.L"d^{-1}", ylabel = "spectral density", xlim = (0, 4*pi), ylim=(0,150) )

Plots.savefig(path*"spectral_density_angle_averaged_high_t.png")


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
            Plots.scatter!(order_metrics_dicts[j]["total_keating_energy_vec"] ./ nr_vertices[j], order_metrics_dicts[j][order_metrics_names[i]] ./ Statistics.mean(order_metrics_dicts[j][order_metrics_names[i]]),
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


path = raw"..\..\presentations\material\\"

# Set the size of the checkerboard
n_rows, n_cols = 8, 8  # 8x8 checkerboard

# Create an empty grid to hold the colors
checkerboard = [Random.rand() for i in 1:n_rows, j in 1:n_cols]

colormap = Plots.cgrad(:roma)

checkerboard_colors = [colormap[checkerboard[i,j]] for i in 1:n_rows, j in 1:n_cols]

# Plot the checkerboard
Plots.plot(Plots.heatmap(1:n_rows, 1:n_cols, checkerboard_colors, aspect_ratio=:equal))
Plots.savefig(path*"checkerboard.png")

colormap = Plots.cgrad(:roma, scale=:lin)
checkerboard_fft = FFTW.fft(checkerboard)
checkerboard_fft_normalized = abs.(checkerboard_fft) ./ 4
checkerboard_fft_colors = [colormap[abs(checkerboard_fft_normalized[i,j])] for i in 1:n_rows, j in 1:n_cols]
Plots.plot(Plots.heatmap(1:n_rows, 1:n_cols, checkerboard_fft_colors, aspect_ratio=:equal))
Plots.savefig(path*"checkerboard_fft.png")