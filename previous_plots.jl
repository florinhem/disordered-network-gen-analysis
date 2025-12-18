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

fontsize=16

Plots.gr()
#Plots.plotlyjs()
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
framestyle = :box,
fontfamily="DejaVu Sans")

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
    S * N * "π\\" * repr(d)
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




bond_bending_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_vec
    
    plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\heat_cool_bond_bending_"*string(bond_bending)*"\\"

    diamonds_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

    diamond_order_metrics_dicts = [GU.load_h5_dict(diamonds_analysis_data_path*"216_vertices_perfect_diamond_small_scale_order_metrics.h5"),
    GU.load_h5_dict(diamonds_analysis_data_path*"512_vertices_perfect_diamond_small_scale_order_metrics.h5"),
    GU.load_h5_dict(diamonds_analysis_data_path*"1000_vertices_perfect_diamond_small_scale_order_metrics.h5")]


    analysis_data_paths = [
    raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_"*string(bond_bending)*"\\",
    raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_"*string(bond_bending)*"\\",
    raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_"*string(bond_bending)*"\\"]

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





function heat_cool_temperature_vec(x, max_temp, heat_cool_rate)

    if x < max_temp/heat_cool_rate 
        return heat_cool_rate * x
    elseif x < 2*max_temp/heat_cool_rate 
        return 2* max_temp - heat_cool_rate*x
    else
        return 0
    end
end

path = raw"..\..\presentations\material\\"

x_vec = collect(0:0.01:5)
y_vec = heat_cool_temperature_vec.(x_vec, 0.1, 0.1)

min_temp = 0.1
max_temp = 0.5
temperatures = [0.1, 0.125, 0.15,  0.2, 0.25, 0.3, 0.4, 0.5]
normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

Plots.plot()

Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.2, 0.1),  alpha=1.0, color=mapped_colors[4])
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.1, 0.1),  alpha=1.0, color=mapped_colors[1])
Plots.plot!(legend = false, xlabel="Attempts per bond chain", ylabel=Latex.L"kT", right_margin = 2Plots.mm, size = (500, 400))

Plots.savefig(path*"heat_cool_0.1_temperature_profile_3.png")


x = collect(-0.2:0.01:0.6)

plot_1 = min.(1, exp.( .-  x ./ 0.1  ) )  
plot_2 = min.(1, exp.( .-  x ./ 0.2  ) )  

myplot = Plots.plot(x, plot_2, label = Latex.L"kT=0.2", linecolor=mapped_colors[4])
myplot = Plots.plot!(x, plot_1, label = Latex.L"kT=0.1", linecolor=mapped_colors[1])

Plots.plot!(grid=false, xlabel="Energy difference",
ylabel = "Acceptance probability")

Plots.savefig(path*"boltzmann_metropolis_temperatures_0.1_0.2.png")



plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\heat_cool_bond_bending_comparison\\"

bond_bending_vec = [0.21, 0.285, 0.36]

nr_vertices_vec = [216, 512, 1000]

markershapes = [:circle, :rect, :diamond]

markersizes = [3, 5, 7]

color_palette = Plots.palette(:tab10)

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "cluster_metric_vec", 
    "pore_size_distribution_second_moment_vec",
    "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec"]


order_metrics_labels = ["Bond length std", "Bond angle std", "Dihedral angle std", "Cluster metric",
    "Pore size dist. 2nd m.",
    "Anisotropy from s. f.", "Anisotropy from s. d."]



order_metrics_dict_arr = Array{Dict{Any, Any}, 2}(undef, length(bond_bending_vec), length(nr_vertices_vec))

for j in eachindex(bond_bending_vec) 

    # loop through folders and append all order metrics to the order_metrics_dict
    for k in eachindex(nr_vertices_vec)

        analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\\"*string(nr_vertices_vec[k])*"_vertices_bond_bending_"*string(bond_bending_vec[j])*"\\"

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

        order_metrics_dict_arr[j, k] = order_metrics_dict

    end

end

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], order_metrics_dict_arr[j, k]["bond_angle_std_vec"], 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Bond angle st. d. / rad", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_bending.png")



Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], sqrt.(order_metrics_dict_arr[j, k]["pore_size_distribution_second_moment_vec"] ), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Critical pore radius / "*Latex.L"d", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_pore_radius.png")

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"] ./ minimum(order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"]), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Anisotropy metric", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_anisotropy_spectral_density.png")

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"], order_metrics_dict_arr[j, k]["bond_length_std_vec"], 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(legend = false, xlabel = "Bond angle st. d. / rad", ylabel = "Bond length st. d. / "*Latex.L"d", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_bending_stretching.png")

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"], sqrt.( order_metrics_dict_arr[j, k]["pore_size_distribution_second_moment_vec"] ), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(legend = false,xlabel = "Bond angle st. d. / rad", ylabel = "Critical pore radius / "*Latex.L"d", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_bending_pore_radius.png")


Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"], order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"] ./ minimum(order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"]), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(legend = false,xlabel = "Bond angle st. d. / rad", ylabel = "Anisotropy metric", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_bending_anisotropy_spectral_density.png")


plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\\"

bond_bending_vec = [0.21, 0.285, 0.36]
nr_vertices_vec = [216, 512, 1000]

markershapes = [:circle, :rect, :diamond]
markersizes = [3, 5, 7]
order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "cluster_metric_vec", 
    "pore_size_distribution_second_moment_vec",
    "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec"]

order_metrics_dict_arr = Array{Dict{Any, Any}, 2}(undef, length(bond_bending_vec), length(nr_vertices_vec))

for j in eachindex(bond_bending_vec) 

    # loop through folders and append all order metrics to the order_metrics_dict
    for k in eachindex(nr_vertices_vec)

        analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\\"*string(nr_vertices_vec[k])*"_vertices_bond_bending_"*string(bond_bending_vec[j])*"\\"

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

        order_metrics_dict_arr[j, k] = order_metrics_dict

    end

end

delta_vec = []

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        # get vector that contains each temperature only once
        unique_temperatures = unique(temperatures)

        # loop through the unique temperatures get all order_metrics dicts that have the same temperature
        for temperature in unique_temperatures

            indices = findall(x->x==temperature, temperatures)

            # get the order_metrics dict for the current temperature
            order_metrics_dict = Dict()
            for order_metric_name in order_metrics_names
                order_metrics_dict[order_metric_name] = vcat([order_metrics_dict_arr[j, k][order_metric_name][i] for i in indices]...)
            end

            # get all combinations of the networks with the current temperature
            combinations = collect(
                Combinatorics.combinations(1:length(indices), 2))

            # for each combination, calculate the statistical difference metric
            # delta
            for combination in combinations
                delta = NA.get_statistical_difference(order_metrics_dict["bond_length_std_vec"][combination[1]],
                    order_metrics_dict["bond_length_std_vec"][combination[2]],
                    order_metrics_dict["bond_angle_std_vec"][combination[1]],
                    order_metrics_dict["bond_angle_std_vec"][combination[2]])

                push!(delta_vec, delta)
            end
        end
    end
end

# plot the delta values
Plots.histogram(delta_vec, xlabel="delta", ylabel="count", legend=false, right_margin = 3Plots.mm,)
Plots.xlims!(0, 1)

# save the plot
Plots.savefig(plots_save_path * "delta_histogram.png")



load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks_th_fluct\216_vertices_bond_bending_0.285\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks_th_fluct\216_vertices_bond_bending_0.285\\"

# plot the bond angle and bond length standard deviation evolution for all
# temperatures in two plots, one for low temperatures and one for high temperatures
temperature_vec = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.2, 0.22, 0.24]

# create a color map for the temperatures from blue to red
min_temp = minimum(temperature_vec)
max_temp = maximum(temperature_vec)
normalized_temperatures = (temperature_vec .- min_temp) ./ (max_temp - min_temp)
colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

split_index = 6

# low temperatures
p1 = Plots.plot()
p2 = Plots.plot()

for i in 1:length(temperature_vec[1:split_index])

    temperature = temperature_vec[i]

    # load dictionary for current temperature
    order_metric_evolution_dict = GU.load_h5_dict(load_path*"T_$(temperature_vec[i])_order_metric_evolution.h5")

    filename_integers = order_metric_evolution_dict["filename_integers"]
    bond_angle_std_vec = order_metric_evolution_dict["bond_angle_std_vec"]
    bond_length_std_vec = order_metric_evolution_dict["bond_length_std_vec"]

    Plots.plot!(p1, filename_integers ./ 100, bond_angle_std_vec, label="T = $(temperature)", color=mapped_colors[i])
    Plots.plot!(p2, filename_integers ./ 100, bond_length_std_vec, label="T = $(temperature)", color=mapped_colors[i])

end

Plots.hline!(p1, [0.3665], label="Weevil", color=:black, linestyle=:dash)
Plots.hline!(p2, [0.18], label="Weevil", color=:black, linestyle=:dash)

Plots.plot!(p1, xlabel="MC step", ylabel="Bond angle st. d.")
Plots.plot!(p2, xlabel="MC step", ylabel="Bond length st. d.")

# save both plots
Plots.savefig(p1, save_path*"bond_angle_std_evolution_low_temperatures.png")
Plots.savefig(p2, save_path*"bond_length_std_evolution_low_temperatures.png")


# high temperatures
p3 = Plots.plot()
p4 = Plots.plot()

for i in 1:length(temperature_vec[split_index+1:end])

    temperature = temperature_vec[split_index+i]

    # load dictionary for current temperature
    order_metric_evolution_dict = GU.load_h5_dict(load_path*"T_$(temperature_vec[split_index+i])_order_metric_evolution.h5")

    filename_integers = order_metric_evolution_dict["filename_integers"]
    bond_angle_std_vec = order_metric_evolution_dict["bond_angle_std_vec"]
    bond_length_std_vec = order_metric_evolution_dict["bond_length_std_vec"]

    Plots.plot!(p3, filename_integers ./ 100, bond_angle_std_vec, label="T = $(temperature)", color=mapped_colors[split_index+i])
    Plots.plot!(p4, filename_integers ./ 100, bond_length_std_vec, label="T = $(temperature)", color=mapped_colors[split_index+i])

end

Plots.hline!(p3, [0.3665], label="Weevil", color=:black, linestyle=:dash)
Plots.hline!(p4, [0.18], label="Weevil", color=:black, linestyle=:dash)

Plots.plot!(p3, xlabel="MC step", ylabel="Bond angle st. d.")
Plots.plot!(p4, xlabel="MC step", ylabel="Bond length st. d.")

# save both plots
Plots.savefig(p3, save_path*"bond_angle_std_evolution_high_temperatures.png")
Plots.savefig(p4, save_path*"bond_length_std_evolution_high_temperatures.png")#%%


using Base.Iterators

function get_bond_length_AND_keating_energy_per_vertex_after_evolve_network(
    ;maximal_temperature=1,
    nr_vertices=nr_vertices)

    evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
    spatial_network = NG.get_periodic_network(evolution_dict)
    #NG.plot_spatial_network(spatial_network)

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
        NA.get_temperature_sequence_heating_cooling_gradient(
            maximal_temperature;
            temperature_gradient = 10.0, 
            nr_monte_carlo_steps_per_temperature = 2/(18*216),
            quench = false)

    evolution_dict["temperature_vec"] = temperature_vec
    evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = nr_monte_carlo_steps_per_temperature_vec

    spatial_network, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(
        spatial_network,
        evolution_dict;
        print_progress = false,
        print_every_nr_attempted_bond_switches = 10)

    #=
    save_path = raw".\my_networks\\"
    filename = "Test3"

    NG.save_spatial_network_to_gml(
        spatial_network,
        filename;
        evolution_dict = evolution_dict,
        save_path = save_path)
    =#

    #Now do here a standard deviation of bond
    bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)

    #Get the Keating energy
    total_energy_keating=NG.get_total_energy_keating(spatial_network)
    energy_keating_per_vertex=total_energy_keating/evolution_dict["nr_vertices"]

    return bond_length_std,energy_keating_per_vertex
end


bond_bending_const=0.285
x=[]
y=[]

nr_trials_per_temperature=1
maximal_temperature_range=0.1:0.1:0.3
#temperature_color_array=[]
#maximal_temperature_array=collect(take(cycle(maximal_temperature_range),nr_trials_per_temperature*length(maximal_temperature_range)))
shape_array=[:diamond, :rect]
nr_vertices_array=[216,512]

@assert length(shape_array)===length(nr_vertices_array)

#=
P=Plots.plot(
    linez=maximal_temperature_range,
    c=:batlow,
    label=false,
    colorbar=:true,
    colorbar_title="Maximal temperature"
)
=#

P=Plots.scatter()

for j in eachindex(nr_vertices_array)

    nr_vertices=nr_vertices_array[j]
    shape=shape_array[j]

    for i in eachindex(maximal_temperature_range)

        maximal_temperature=maximal_temperature_range[i]
        #temperature_color=temperature_color_array[i]

        for k in 1:nr_trials_per_temperature

            println("start: shape="*"$shape"*", N="*"$nr_vertices"*", T_max="*"$maximal_temperature"*", k="*"$k")
            bond_length_std,energy_keating_per_vertex=get_bond_length_AND_keating_energy_per_vertex_after_evolve_network(
                ;maximal_temperature=maximal_temperature,
                nr_vertices=nr_vertices)

            P=Plots.scatter!(
                P,
                [energy_keating_per_vertex],
                [bond_length_std],
                xlabel="Keating energy per vertex",
                ylabel="Bond length std",
                
                #marker_z=maximal_temperature_range,
                #markercolor=[]
                markershape=shape,

                zcolor=maximal_temperature_range,
                c=:batlow,
                
                label=false,
                colorbar=:true,
                colorbar_title="Maximal temperature",
                
                flip_axis=false
                )
        end
    end
end

minimum_temperature=minimum(maximal_temperature_range)
maximum_temperature=maximum(maximal_temperature_range)

minimum_nr_vertices=minimum(nr_vertices_array)
maximum_nr_vertices=maximum(nr_vertices_array)

Plots.savefig(P,raw".\my_networks\KE_VS_BLSTD\KE_VS_BLSTD_16"*
    "_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"*
    "_beta="*"$bond_bending_const"*
    "_T="*"$minimum_temperature" * "-" * "$maximum_temperature"*
    "_trials="*"$nr_trials_per_temperature"*
    ".png")



function heat_cool_temperature_vec(x, max_temp, heat_cool_rate)

    if x < max_temp/heat_cool_rate 
        return heat_cool_rate * x
    elseif x < 2*max_temp/heat_cool_rate 
        return 2* max_temp - heat_cool_rate*x
    else
        return 0
    end
end

path = raw"..\..\presentations\material\\"

x_vec = collect(0:0.01:5)

color_low_T = "#6c72f5"
color_high_T = "#b9b681"

Plots.plot()

Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.17, 0.1),  alpha=1.0, color=color_high_T)
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.11, 0.1),  alpha=1.0, color=color_low_T)
Plots.plot!(legend = false, xlabel="Attempts per bond chain", ylabel=Latex.L"kT", right_margin = 2Plots.mm, size = (450, 300))

Plots.savefig(path*"heat_cool_0.1_temperature_profile_4.png")


x = collect(-0.1:0.01:0.4)

plot_1 = min.(1, exp.( .-  x ./ 0.11  ) )  
plot_2 = min.(1, exp.( .-  x ./ 0.17  ) )  

myplot = Plots.plot(x, plot_2, label = Latex.L"kT=0.17", linecolor=color_high_T)
myplot = Plots.plot!(x, plot_1, label = Latex.L"kT=0.11", linecolor=color_low_T)

Plots.plot!(grid=false, xlabel="Energy difference",
ylabel = "Acceptance prob.", size = (450, 300))

Plots.savefig(path*"boltzmann_metropolis_temperatures_0.11_0.17.png")



plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\heat_cool_bond_bending_comparison\\"

bond_bending_vec = [0.21, 0.285, 0.36]

nr_vertices_vec = [216, 512, 1000]

markershapes = [:circle, :rect, :diamond]

markersizes = [3, 5, 7]

color_palette = Plots.palette(:tab10)

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "cluster_metric_vec", 
    "pore_size_distribution_second_moment_vec",
    "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec"]


order_metrics_labels = ["Bond length std", "Bond angle std", "Dihedral angle std", "Cluster metric",
    "Pore size dist. 2nd m.",
    "Anisotropy from s. f.", "Anisotropy from s. d."]



order_metrics_dict_arr = Array{Dict{Any, Any}, 2}(undef, length(bond_bending_vec), length(nr_vertices_vec))

for j in eachindex(bond_bending_vec) 

    # loop through folders and append all order metrics to the order_metrics_dict
    for k in eachindex(nr_vertices_vec)

        analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\\"*string(nr_vertices_vec[k])*"_vertices_bond_bending_"*string(bond_bending_vec[j])*"\\"

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

        order_metrics_dict_arr[j, k] = order_metrics_dict

    end

end

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Bond angle st. d. / °", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_bending.png")



Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], sqrt.(order_metrics_dict_arr[j, k]["pore_size_distribution_second_moment_vec"] ), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Critical pore radius / "*Latex.L"d", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_pore_radius.png")


Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, sqrt.(order_metrics_dict_arr[j, k]["pore_size_distribution_second_moment_vec"] ), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Critical pore radius / "*Latex.L"d", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_bending_pore_radius.png")

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"] ./ minimum(order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"]), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Anisotropy metric", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_anisotropy_spectral_density.png")

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"] ./ minimum(order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"]), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Anisotropy metric", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_bending_anisotropy_spectral_density.png")

Plots.scatter()
#Plots.plot()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["bond_length_std_vec"], markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
        #Plots.plot!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["bond_length_std_vec"], seriestype = :scatter, markershape = markershapes[j],  color=mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
        #Plots.plot!(x, y, zcolor = z, seriestype = :scatter, markersize=5, label = "points")zcolor = temperatures,
    end
end

min_temp = 0.08
max_temp = 0.26
Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Bond length st. d. / "*Latex.L"d", rightmargin=5Plots.mm, colorbar=true, clim=(min_temp, max_temp), color=Plots.cgrad(:roma, rev = true, scale = :exp))
#Plots.plot!(size=(550,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Bond length st. d. / "*Latex.L"d", rightmargin=5Plots.mm, clim=(min_temp, max_temp), color_palette=Plots.cgrad(:roma, rev = true, scale = :exp))

Plots.savefig(plots_save_path*"bond_bending_stretching.png")





plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_heat_cool_bond_bending_comparison\\"

bond_bending_vec = [0.06, 0.135, 0.21, 0.285, 0.36, 0.435, 0.51]

nr_vertices_vec = [216]

markershapes = [:utriangle, :pentagon, :diamond,  :circle, :star5, :rect, :dtriangle]

markersizes = [3]

color_palette = Plots.palette(:tab10)

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "cluster_metric_vec", 
    "pore_size_distribution_second_moment_vec",
    "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec"]


order_metrics_labels = ["Bond length std", "Bond angle std", "Dihedral angle std", "Cluster metric",
    "Pore size dist. 2nd m.",
    "Anisotropy from s. f.", "Anisotropy from s. d."]


order_metrics_dict_arr = Array{Dict{Any, Any}, 2}(undef, length(bond_bending_vec), length(nr_vertices_vec))

for j in eachindex(bond_bending_vec) 

    # loop through folders and append all order metrics to the order_metrics_dict
    for k in eachindex(nr_vertices_vec)

        analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\\"*string(nr_vertices_vec[k])*"_vertices_bond_bending_"*string(bond_bending_vec[j])*"\\"

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

        order_metrics_dict_arr[j, k] = order_metrics_dict

    end

end

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Bond angle st. d. / °", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_bending.png")



Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], sqrt.(order_metrics_dict_arr[j, k]["pore_size_distribution_second_moment_vec"] ), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), ylims=(0, 4.8), legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Critical pore radius / "*Latex.L"d", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_pore_radius.png")


Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, sqrt.(order_metrics_dict_arr[j, k]["pore_size_distribution_second_moment_vec"] ), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Critical pore radius / "*Latex.L"d", rightmargin=5Plots.mm, ylims=(0, 4.8))

Plots.savefig(plots_save_path*"bond_bending_pore_radius.png")

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_length_std_vec"], order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"] ./ minimum(order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"]), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond length st. d. / "*Latex.L"d", ylabel = "Anisotropy metric", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_stretching_anisotropy_spectral_density.png")

Plots.scatter()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)

        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"] ./ minimum(order_metrics_dict_arr[j, k]["anisotropy_metric_from_spectral_density_vec"]), 
        markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
    end
end

Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Anisotropy metric", rightmargin=5Plots.mm)

Plots.savefig(plots_save_path*"bond_bending_anisotropy_spectral_density.png")

Plots.scatter()
#Plots.plot()

for j in eachindex(bond_bending_vec) 
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["bond_length_std_vec"], markershape = markershapes[j], color = mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
        #Plots.plot!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["bond_length_std_vec"], seriestype = :scatter, markershape = markershapes[j],  color=mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
        #Plots.plot!(x, y, zcolor = z, seriestype = :scatter, markersize=5, label = "points")zcolor = temperatures,
    end
end

min_temp = 0.08
max_temp = 0.26
Plots.scatter!(size=(470,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Bond length st. d. / "*Latex.L"d", rightmargin=5Plots.mm, colorbar=true, clim=(min_temp, max_temp), color=Plots.cgrad(:roma, rev = true, scale = :exp))
#Plots.plot!(size=(550,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Bond length st. d. / "*Latex.L"d", rightmargin=5Plots.mm, clim=(min_temp, max_temp), color_palette=Plots.cgrad(:roma, rev = true, scale = :exp))

Plots.savefig(plots_save_path*"bond_bending_stretching.png")



path_1 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.135\run_1\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_pore_size_distribution.h5"

path_2 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.21\run_1\216_vertices_T_0.11_heat_cool_0.1_per_mc_quenched_pore_size_distribution.h5"

dict_1 = GU.load_h5_dict(path_1)
dict_2 = GU.load_h5_dict(path_2)

second_moment_1 = NA.get_pore_size_distribution_second_moment(dict_1)
second_moment_2 = NA.get_pore_size_distribution_second_moment(dict_2)

println("second_moment_1 = ", second_moment_1)
println("second_moment_2 = ", second_moment_2)

# print square roots of the second moments
println("sqrt(second_moment_1) = ", sqrt(second_moment_1))
println("sqrt(second_moment_2) = ", sqrt(second_moment_2))

Plots.plot(dict_1["pore_size_vec"], dict_1["pore_size_distribution"])
Plots.plot!(dict_2["pore_size_vec"], dict_2["pore_size_distribution"])



path = raw"..\..\presentations\material\\"

x_vec = collect(0:0.01:2)
y_vec = collect(0.5:0.01:1.5)
Plots.plot(x_vec, (3/16) .* (x_vec.^2 .- 1).^2   )
Plots.plot!(xlabel="Bond length / "*Latex.L"d", ylabel="Energy", right_margin = 3Plots.mm, ylims=(0,0.3), xlims=(0,2), legend=false)

Plots.savefig(path*"bond_stretching_energy_2.png")


x_vec = collect(0:0.1:180)
y_vec = collect(40:0.1:180)
Plots.plot(x_vec, (3/8 * 0.285) .* (cosd.(x_vec) .+ 1/3).^2  )
Plots.plot!(xlabel="Bond angle / °", ylabel="Energy", right_margin = 5Plots.mm, ylims=(0,0.3), xlims=(0,180), legend=false)
Plots.savefig(path*"bond_bending_energy_2.png")


plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_heat_cool_bond_bending_comparison\\"

bond_bending_vec = [0.06, 0.135, 0.21, 0.285, 0.36, 0.435, 0.51]

nr_vertices_vec = [216]

markershapes = [:utriangle, :pentagon, :diamond,  :circle, :star5, :rect, :dtriangle]

markersizes = [4]

color_palette = Plots.palette(:tab10)

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "cluster_metric_vec", 
    "pore_size_distribution_second_moment_vec",
    "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec"]


order_metrics_labels = ["Bond length std", "Bond angle std", "Dihedral angle std", "Cluster metric",
    "Pore size dist. 2nd m.",
    "Anisotropy from s. f.", "Anisotropy from s. d."]


order_metrics_dict_arr = Array{Dict{Any, Any}, 2}(undef, length(bond_bending_vec), length(nr_vertices_vec))

for j in eachindex(bond_bending_vec) 

    # loop through folders and append all order metrics to the order_metrics_dict
    for k in eachindex(nr_vertices_vec)

        analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\old_order_metrics\\"*string(nr_vertices_vec[k])*"_vertices_bond_bending_"*string(bond_bending_vec[j])*"\\"

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

        order_metrics_dict_arr[j, k] = order_metrics_dict

    end

end


Plots.scatter()
#Plots.plot()

for j in reverse(eachindex(bond_bending_vec) )
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["bond_length_std_vec"], markershape = markershapes[j], color = mapped_colors, markersize=markersizes[k], alpha = 0.5, label = false)
        #Plots.plot!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["bond_length_std_vec"], seriestype = :scatter, markershape = markershapes[j],  color=mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
        #Plots.plot!(x, y, zcolor = z, seriestype = :scatter, markersize=5, label = "points")zcolor = temperatures,

        
        Plots.scatter!([], [], color = "black", markershape = markershapes[j], label = Latex.L"\beta = "*Latex.latexstring(bond_bending_vec[j]), markersize=markersizes[k], alpha = 0.5, legend = true)

    end

end

Plots.scatter!([15.3, 13.6, 13.5], [0.045, 0.047, 0.066], color = ["#7BA8D6","#4549A9","#8C66EC"], markershape = [:diamond, :circle, :square], label = false, markersize=10, alpha = 1)

min_temp = 0.08
max_temp = 0.26
Plots.scatter!(size=(400,400), legend = true, xlabel = Latex.L" \sigma_\mathrm{angle} / °", ylabel = Latex.L" \sigma_\mathrm{length} / d", rightmargin=5Plots.mm, colorbar=true, clim=(min_temp, max_temp), color=Plots.cgrad(:roma, rev = true, scale = :exp), xlims=(-5, 23), ylims=(0.01, 0.1))
#Plots.plot!(size=(550,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Bond length st. d. / "*Latex.L"d", rightmargin=5Plots.mm, clim=(min_temp, max_temp), color_palette=Plots.cgrad(:roma, rev = true, scale = :exp))

Plots.savefig(plots_save_path*"bond_bending_stretching_highlighted.png")


plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_bond_bending_0.21\run_1\\"


diamond_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\216_vertices_perfect_diamond.gml"

disorder_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched.gml"

disorder_path_2 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\216_vertices_T_0.16_heat_cool_0.1_per_mc_quenched.gml"

disorder_path_3 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\216_vertices_T_0.24_heat_cool_0.1_per_mc_quenched.gml"


diamond = NG.load_spatial_network_from_gml(diamond_path)
disorder = NG.load_spatial_network_from_gml(disorder_path)
disorder_2 = NG.load_spatial_network_from_gml(disorder_path_2)
disorder_3 = NG.load_spatial_network_from_gml(disorder_path_3)


@time pore_size_distribution_dict_diamond = NA.get_pore_size_distribution(diamond, sampling_grid_size = 0.2)
@time pore_size_distribution_dict_disorder = NA.get_pore_size_distribution(disorder, sampling_grid_size = 0.2)
@time pore_size_distribution_dict_disorder_2 = NA.get_pore_size_distribution(disorder_2, sampling_grid_size = 0.2)
@time pore_size_distribution_dict_disorder_3 = NA.get_pore_size_distribution(disorder_3, sampling_grid_size = 0.2)

# plot the pore size distribution for the diamond network
Plots.plot(pore_size_distribution_dict_diamond["pore_size_vec"], pore_size_distribution_dict_diamond["pore_size_distribution"], label = "T=0")
Plots.plot!(pore_size_distribution_dict_disorder["pore_size_vec"], pore_size_distribution_dict_disorder["pore_size_distribution"], label = "T=0.1")
Plots.plot!(pore_size_distribution_dict_disorder_2["pore_size_vec"], pore_size_distribution_dict_disorder_2["pore_size_distribution"], label = "T=0.16")
Plots.plot!(pore_size_distribution_dict_disorder_3["pore_size_vec"], pore_size_distribution_dict_disorder_3["pore_size_distribution"], label = "T=0.24")
Plots.xlabel!("Pore size/d")
Plots.ylabel!("Pore size distribution")

Plots.savefig(plots_save_path*"pore_size_distribution.png")


max_ring_size_to_check = 10

ring_size_distribution_diamond = NA.get_ring_size_distribution(diamond, max_ring_size_to_check=max_ring_size_to_check)
ring_size_distribution_disorder = NA.get_ring_size_distribution(disorder, max_ring_size_to_check=max_ring_size_to_check)
ring_size_distribution_disorder_2 = NA.get_ring_size_distribution(disorder_2, max_ring_size_to_check=max_ring_size_to_check)
ring_size_distribution_disorder_3 = NA.get_ring_size_distribution(disorder_3, max_ring_size_to_check=max_ring_size_to_check)

# plot the ring size distribution for all networks
Plots.plot(ring_size_distribution_diamond["ring_size_vec"], ring_size_distribution_diamond["ring_size_distribution"], label = "T=0")
Plots.plot!(ring_size_distribution_disorder["ring_size_vec"], ring_size_distribution_disorder["ring_size_distribution"], label = "T=0.1")
Plots.plot!(ring_size_distribution_disorder_2["ring_size_vec"], ring_size_distribution_disorder_2["ring_size_distribution"], label = "T=0.16")
Plots.plot!(ring_size_distribution_disorder_3["ring_size_vec"], ring_size_distribution_disorder_3["ring_size_distribution"], label = "T=0.24")
Plots.xlabel!("Ring size")
Plots.ylabel!("Ring size distribution")

Plots.savefig(plots_save_path*"ring_size_distribution.png")

ring_radius_distribution_diamond = NA.get_ring_radius_distribution(diamond, ring_size_distribution_diamond)
ring_radius_distribution_disorder = NA.get_ring_radius_distribution(disorder, ring_size_distribution_disorder)
ring_radius_distribution_disorder_2 = NA.get_ring_radius_distribution(disorder_2, ring_size_distribution_disorder_2)
ring_radius_distribution_disorder_3 = NA.get_ring_radius_distribution(disorder_3, ring_size_distribution_disorder_3)

# plot the ring size distribution for all networks
Plots.plot(ring_radius_distribution_diamond["ring_radius_vec"], ring_radius_distribution_diamond["ring_radius_distribution"], label = "T=0")
Plots.plot!(ring_radius_distribution_disorder["ring_radius_vec"], ring_radius_distribution_disorder["ring_radius_distribution"], label = "T=0.1")
Plots.plot!(ring_radius_distribution_disorder_2["ring_radius_vec"], ring_radius_distribution_disorder_2["ring_radius_distribution"], label = "T=0.16")
Plots.plot!(ring_radius_distribution_disorder_3["ring_radius_vec"], ring_radius_distribution_disorder_3["ring_radius_distribution"], label = "T=0.24")
Plots.xlabel!("Ring radius / d")
Plots.ylabel!("Ring radius distribution")

Plots.savefig(plots_save_path*"ring_radius_distribution.png")



plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_bond_bending_0.21\run_1\\"


diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_structure_factor_array.h5")
diamond_vertices = GU.load_h5_dict(plots_save_path*"diamond_vertices_structure_factor_array.h5")
disorder_bonds = GU.load_h5_dict(plots_save_path*"disorder_bonds_structure_factor_array.h5")
disorder_vertices = GU.load_h5_dict(plots_save_path*"disorder_vertices_structure_factor_array.h5")
disorder_2_bonds = GU.load_h5_dict(plots_save_path*"disorder_2_bonds_structure_factor_array.h5")
disorder_2_vertices = GU.load_h5_dict(plots_save_path*"disorder_2_vertices_structure_factor_array.h5")
disorder_3_bonds = GU.load_h5_dict(plots_save_path*"disorder_3_bonds_structure_factor_array.h5")
disorder_3_vertices = GU.load_h5_dict(plots_save_path*"disorder_3_vertices_structure_factor_array.h5")

upper_clim_bonds = 1.0
upper_clim_vertices = 5.0

NA.plot_structure_factor_heatmap(
    diamond_bonds,
    plots_save_path*"diamond_bonds";
    title="Diamond bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    diamond_vertices,
    plots_save_path*"diamond_vertices";
    title="Diamond vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_bonds,
    plots_save_path*"disorder_bonds";
    title="Disorder bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_vertices,
    plots_save_path*"disorder_vertices";
    title="Disorder vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_2_bonds,
    plots_save_path*"disorder_2_bonds";
    title="Disorder 2 bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_2_vertices,
    plots_save_path*"disorder_2_vertices";
    title="Disorder 2 vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_3_bonds,
    plots_save_path*"disorder_3_bonds";
    title="Disorder 3 bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_3_vertices,
    plots_save_path*"disorder_3_vertices";
    title="Disorder 3 vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    x_y_lims = nothing)


plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_bond_bending_0.21\run_1\\"


diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_structure_factor_angle_averaged.h5")
diamond_vertices = GU.load_h5_dict(plots_save_path*"diamond_vertices_structure_factor_angle_averaged.h5")
disorder_bonds = GU.load_h5_dict(plots_save_path*"disorder_bonds_structure_factor_angle_averaged.h5")
disorder_vertices = GU.load_h5_dict(plots_save_path*"disorder_vertices_structure_factor_angle_averaged.h5")
disorder_2_bonds = GU.load_h5_dict(plots_save_path*"disorder_2_bonds_structure_factor_angle_averaged.h5")
disorder_2_vertices = GU.load_h5_dict(plots_save_path*"disorder_2_vertices_structure_factor_angle_averaged.h5")
disorder_3_bonds = GU.load_h5_dict(plots_save_path*"disorder_3_bonds_structure_factor_angle_averaged.h5")
disorder_3_vertices = GU.load_h5_dict(plots_save_path*"disorder_3_vertices_structure_factor_angle_averaged.h5")

# Plot angle averaged structure factors considering bonds

Plots.plot(diamond_bonds["wavenumber_vec"], 
    Measurements.value.(diamond_bonds["structure_factor_vec"]) , 
    ribbon =  Measurements.uncertainty.(diamond_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[1])
Plots.plot!(disorder_bonds["wavenumber_vec"],
    Measurements.value.(disorder_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=0.1", fillcolor = Plots.palette(:tab10)[2])
Plots.plot!(disorder_2_bonds["wavenumber_vec"],
    Measurements.value.(disorder_2_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_2_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0.16", fillcolor = Plots.palette(:tab10)[3])

Plots.plot!(xlabel="Wavenumber / "*Latex.L"d^{-1}", ylabel = "Structure factor (bonds)",  xtick=pitick(0, 32, 1; mode=:latex),
    xlim = (0, 25), ylims = (0, 8)) #, leftmargin = 0Plots.mm
    Plots.savefig(plots_save_path*"angle_averaged_structure_factor_bonds_no_T_0.24.png")

Plots.plot!(disorder_3_bonds["wavenumber_vec"],
    Measurements.value.(disorder_3_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_3_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[4], label = Latex.L"kT=0.24", fillcolor = Plots.palette(:tab10)[4])

Plots.savefig(plots_save_path*"angle_averaged_structure_factor_bonds.png")


Plots.plot(diamond_vertices["wavenumber_vec"], 
    Measurements.value.(diamond_vertices["structure_factor_vec"]) , 
    ribbon =  Measurements.uncertainty.(diamond_vertices["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[1])
Plots.plot!(disorder_vertices["wavenumber_vec"],
    Measurements.value.(disorder_vertices["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_vertices["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=0.1", fillcolor = Plots.palette(:tab10)[2])
Plots.plot!(disorder_2_vertices["wavenumber_vec"],
    Measurements.value.(disorder_2_vertices["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_2_vertices["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0.16", fillcolor = Plots.palette(:tab10)[3])


Plots.plot!(xlabel="Wavenumber / "*Latex.L"d^{-1}", ylabel = "Structure factor (vertices)",  xtick=pitick(0, 32, 1; mode=:latex), xlim = (0, 25), ylims = (0, 40)) #, leftmargin = 0Plots.mm
Plots.savefig(plots_save_path*"angle_averaged_structure_factor_vertices_no_T_0.24.png")

Plots.plot!(disorder_3_vertices["wavenumber_vec"],
    Measurements.value.(disorder_3_vertices["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_3_vertices["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[4], label = Latex.L"kT=0.24", fillcolor = Plots.palette(:tab10)[4])

Plots.savefig(plots_save_path*"angle_averaged_structure_factor_vertices.png")


diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_structure_factor_array.h5")
diamond_vertices = GU.load_h5_dict(plots_save_path*"diamond_vertices_structure_factor_array.h5")
disorder_bonds = GU.load_h5_dict(plots_save_path*"disorder_bonds_structure_factor_array.h5")
disorder_vertices = GU.load_h5_dict(plots_save_path*"disorder_vertices_structure_factor_array.h5")
disorder_2_bonds = GU.load_h5_dict(plots_save_path*"disorder_2_bonds_structure_factor_array.h5")
disorder_2_vertices = GU.load_h5_dict(plots_save_path*"disorder_2_vertices_structure_factor_array.h5")
disorder_3_bonds = GU.load_h5_dict(plots_save_path*"disorder_3_bonds_structure_factor_array.h5")
disorder_3_vertices = GU.load_h5_dict(plots_save_path*"disorder_3_vertices_structure_factor_array.h5")

upper_clim_bonds = 1.0
upper_clim_vertices = 5.0

NA.plot_structure_factor_heatmap(
    diamond_bonds,
    plots_save_path*"diamond_bonds_x_z";
    title="Diamond bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    diamond_vertices,
    plots_save_path*"diamond_vertices_x_z";
    title="Diamond vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_bonds,
    plots_save_path*"disorder_bonds_x_z";
    title="Disorder bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_vertices,
    plots_save_path*"disorder_vertices_x_z";
    title="Disorder vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_2_bonds,
    plots_save_path*"disorder_2_bonds_x_z";
    title="Disorder 2 bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_2_vertices,
    plots_save_path*"disorder_2_vertices_x_z";
    title="Disorder 2 vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_3_bonds,
    plots_save_path*"disorder_3_bonds_x_z";
    title="Disorder 3 bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_3_vertices,
    plots_save_path*"disorder_3_vertices_x_z";
    title="Disorder 3 vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)


plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\diamonds\\"

diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_structure_factor_array.h5")
diamond_vertices = GU.load_h5_dict(plots_save_path*"diamond_vertices_structure_factor_array.h5")

NA.plot_structure_factor_heatmap(
    diamond_bonds,
    plots_save_path*"diamond_bonds_x_z";
    title="Diamond bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    diamond_vertices,
    plots_save_path*"diamond_vertices_x_z";
    title="Diamond vertices",
    save_plot = true,
    clims = (0, upper_clim_vertices ),
    wavevector_component_to_fix = 2,
    x_y_lims = nothing)



diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_structure_factor_angle_averaged.h5")
diamond_vertices = GU.load_h5_dict(plots_save_path*"diamond_vertices_structure_factor_angle_averaged.h5")

Plots.plot(diamond_bonds["wavenumber_vec"], 
    Measurements.value.(diamond_bonds["structure_factor_vec"]) , 
    ribbon =  Measurements.uncertainty.(diamond_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[1])

Plots.plot!(xlabel="Wavenumber / "*Latex.L"d^{-1}", ylabel = "Structure factor (bonds)",  xtick=pitick(0, 32, 1; mode=:latex),
    xlim = (0, 25), ylims = (0, 8)) 

Plots.savefig(plots_save_path*"angle_averaged_structure_factor_bonds.png")

Plots.plot(diamond_vertices["wavenumber_vec"], 
    Measurements.value.(diamond_vertices["structure_factor_vec"]) , 
    ribbon =  Measurements.uncertainty.(diamond_vertices["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[1])

Plots.plot!(xlabel="Wavenumber / "*Latex.L"d^{-1}", ylabel = "Structure factor (vertices)",  xtick=pitick(0, 32, 1; mode=:latex), xlim = (0, 25), ylims = (0, 10))

Plots.savefig(plots_save_path*"angle_averaged_structure_factor_vertices.png")


plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_bond_bending_0.21\run_1\\"

diamond = GU.load_h5_dict(plots_save_path*"diamond_spectral_density_array.h5")
disorder = GU.load_h5_dict(plots_save_path*"disorder_spectral_density_array.h5")
disorder_2 = GU.load_h5_dict(plots_save_path*"disorder_2_spectral_density_array.h5")
disorder_3 = GU.load_h5_dict(plots_save_path*"disorder_3_spectral_density_array.h5")

upper_clim = 0.1

NA.plot_spectral_density_heatmap(
    diamond,
    save_path*"diamond";
    title="Diamond",
    save_plot = true,
    clims = (0, upper_clim),
    wavevector_component_to_fix = 2,
    wavevector_value_fixed = 0,
    plot_im_re = false)

NA.plot_spectral_density_heatmap(
    disorder,
    save_path*"disorder";
    title="Disorder",
    save_plot = true,
    clims = (0, upper_clim),
    wavevector_component_to_fix = 2,
    wavevector_value_fixed = 0,
    plot_im_re = false)

NA.plot_spectral_density_heatmap(
    disorder_2,
    save_path*"disorder_2";
    title="Disorder 2",
    save_plot = true,
    clims = (0, upper_clim),
    wavevector_component_to_fix = 2,
    wavevector_value_fixed = 0,
    plot_im_re = false)

NA.plot_spectral_density_heatmap(
    disorder_3,
    save_path*"disorder_3";
    title="Disorder 3",
    save_plot = true,
    clims = (0, upper_clim),
    wavevector_component_to_fix = 2,
    wavevector_value_fixed = 0,
    plot_im_re = false)


diamond = GU.load_h5_dict(plots_save_path*"diamond_spectral_density_angle_averaged.h5")
disorder = GU.load_h5_dict(plots_save_path*"disorder_spectral_density_angle_averaged.h5")
disorder_2 = GU.load_h5_dict(plots_save_path*"disorder_2_spectral_density_angle_averaged.h5")
disorder_3 = GU.load_h5_dict(plots_save_path*"disorder_3_spectral_density_angle_averaged.h5")

Plots.plot(diamond["wavenumber_vec"], 
    Measurements.value.(diamond["spectral_density_vec"]) , 
    ribbon =  Measurements.uncertainty.(diamond["spectral_density_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[1])
Plots.plot!(disorder["wavenumber_vec"],
    Measurements.value.(disorder["spectral_density_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder["spectral_density_vec"]), linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=0.1", fillcolor = Plots.palette(:tab10)[2])
Plots.plot!(disorder_2["wavenumber_vec"],
    Measurements.value.(disorder_2["spectral_density_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_2["spectral_density_vec"]), linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0.16", fillcolor = Plots.palette(:tab10)[3])
Plots.plot!(disorder_3["wavenumber_vec"],
    Measurements.value.(disorder_3["spectral_density_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_3["spectral_density_vec"]), linecolor=Plots.palette(:tab10)[4], label = Latex.L"kT=0.24", fillcolor = Plots.palette(:tab10)[4])

Plots.plot!(xlabel="Wavenumber / "*Latex.L"d^{-1}", ylabel = "Spectral density",  xtick=pitick(0, 32, 1; mode=:latex), xlim = (0, 25), ylims = (0, 1)) 
Plots.savefig(plots_save_path*"angle_averaged_spectral_density.png")

diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_3_structure_factor_angle_averaged.h5")
disorder_bonds = GU.load_h5_dict(plots_save_path*"disorder_bonds_3_structure_factor_angle_averaged.h5")
disorder_2_bonds = GU.load_h5_dict(plots_save_path*"disorder_2_bonds_3_structure_factor_angle_averaged.h5")
disorder_3_bonds = GU.load_h5_dict(plots_save_path*"disorder_3_bonds_3_structure_factor_angle_averaged.h5")

# Plot angle averaged structure factors considering bonds

Plots.plot(diamond_bonds["wavenumber_vec"], 
    Measurements.value.(diamond_bonds["structure_factor_vec"]) , 
    ribbon =  Measurements.uncertainty.(diamond_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[1])
Plots.plot!(disorder_bonds["wavenumber_vec"],
    Measurements.value.(disorder_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=0.1", fillcolor = Plots.palette(:tab10)[2])
Plots.plot!(disorder_2_bonds["wavenumber_vec"],
    Measurements.value.(disorder_2_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_2_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0.16", fillcolor = Plots.palette(:tab10)[3])
Plots.plot!(disorder_3_bonds["wavenumber_vec"],
    Measurements.value.(disorder_3_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_3_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[4], label = Latex.L"kT=0.24", fillcolor = Plots.palette(:tab10)[4])
Plots.plot!(xlabel="Wavenumber / "*Latex.L"d^{-1}", ylabel = "Structure factor (bonds)",  xtick=pitick(0, 32, 1; mode=:latex),
    xlim = (0, 25), ylims = (0, 6))
Plots.savefig(plots_save_path*"angle_averaged_structure_factor_bonds_3.png")


diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_rot_structure_factor_array.h5")

upper_clim_bonds = 1.0

NA.plot_structure_factor_heatmap(
    diamond_bonds,
    plots_save_path*"diamond_bonds_rot";
    title="Diamond bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    wavevector_component_to_fix = 3,
    x_y_lims = nothing)

diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_rot_structure_factor_angle_averaged.h5")


Plots.plot(diamond_bonds["wavenumber_vec"], 
    Measurements.value.(diamond_bonds["structure_factor_vec"]) , 
    ribbon =  Measurements.uncertainty.(diamond_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[1])
Plots.plot!(xlabel="Wavenumber / "*Latex.L"d^{-1}", ylabel = "Structure factor (bonds)",  xtick=pitick(0, 32, 1; mode=:latex),
    xlim = (0, 25), ylims = (0, 6))
Plots.savefig(plots_save_path*"angle_averaged_structure_factor_bonds_rot.png")


plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_bond_bending_0.21\run_1\\"

diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_structure_factor_array.h5")
disorder_bonds = GU.load_h5_dict(plots_save_path*"disorder_bonds_structure_factor_array.h5")
disorder_2_bonds = GU.load_h5_dict(plots_save_path*"disorder_2_bonds_structure_factor_array.h5")
disorder_3_bonds = GU.load_h5_dict(plots_save_path*"disorder_3_bonds_structure_factor_array.h5")

upper_clim_bonds = 1.0

NA.plot_structure_factor_heatmap(
    diamond_bonds,
    plots_save_path*"diamond_bonds";
    title="Diamond bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)


NA.plot_structure_factor_heatmap(
    disorder_bonds,
    plots_save_path*"disorder_bonds";
    title="Disorder bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)


NA.plot_structure_factor_heatmap(
    disorder_2_bonds,
    plots_save_path*"disorder_2_bonds";
    title="Disorder 2 bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)

NA.plot_structure_factor_heatmap(
    disorder_3_bonds,
    plots_save_path*"disorder_3_bonds";
    title="Disorder 3 bonds",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)


diamond_bonds = GU.load_h5_dict(plots_save_path*"diamond_bonds_structure_factor_angle_averaged.h5")
disorder_bonds = GU.load_h5_dict(plots_save_path*"disorder_bonds_structure_factor_angle_averaged.h5")
disorder_2_bonds = GU.load_h5_dict(plots_save_path*"disorder_2_bonds_structure_factor_angle_averaged.h5")
disorder_3_bonds = GU.load_h5_dict(plots_save_path*"disorder_3_bonds_structure_factor_angle_averaged.h5")

# Plot angle averaged structure factors considering bonds

Plots.plot(diamond_bonds["wavenumber_vec"], 
    Measurements.value.(diamond_bonds["structure_factor_vec"]) , 
    ribbon =  Measurements.uncertainty.(diamond_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=0", fillcolor = Plots.palette(:tab10)[1])
Plots.plot!(disorder_bonds["wavenumber_vec"],
    Measurements.value.(disorder_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=0.1", fillcolor = Plots.palette(:tab10)[2])
Plots.plot!(disorder_2_bonds["wavenumber_vec"],
    Measurements.value.(disorder_2_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_2_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[3], label = Latex.L"kT=0.16", fillcolor = Plots.palette(:tab10)[3])

Plots.plot!(disorder_3_bonds["wavenumber_vec"],
    Measurements.value.(disorder_3_bonds["structure_factor_vec"]), 
    ribbon =  Measurements.uncertainty.(disorder_3_bonds["structure_factor_vec"]), linecolor=Plots.palette(:tab10)[4], label = Latex.L"kT=0.24", fillcolor = Plots.palette(:tab10)[4])
Plots.plot!(xlabel="Wavenumber / "*Latex.L"d^{-1}", ylabel = "Structure factor (bonds)",  xtick=pitick(0, 32, 1; mode=:latex),
    xlim = (0, 25), ylims = (0, 8)) 

Plots.savefig(plots_save_path*"angle_averaged_structure_factor_bonds.png")


save_path =  raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_bond_bending_0.21\run_1\\"
consider_spectral_density = false

diamond_bonds = GU.load_h5_dict(save_path*"diamond_bonds_structure_factor_angle_averaged.h5")
disorder_bonds = GU.load_h5_dict(save_path*"disorder_bonds_structure_factor_angle_averaged.h5")
disorder_2_bonds = GU.load_h5_dict(save_path*"disorder_2_bonds_structure_factor_angle_averaged.h5")
disorder_3_bonds = GU.load_h5_dict(save_path*"disorder_3_bonds_structure_factor_angle_averaged.h5")

t_range = (1e-4, 1e1)

function comparison_func(x, a) 
    return a * x ^(-3/2)
end

# get the vector of times to sample the excess spreadability such that the
# t values are equally spaced on a logarithmic scale
time_vec = exp.(LinRange(log(t_range[1]), log(t_range[2]), 100))
# calculate the excess spreadability for each t value
excess_spreadability_vec_diamond = [NA.get_excess_spreadability(
    diamond_bonds, time_vec[i]; 
    consider_spectral_density = consider_spectral_density) for i in 
    eachindex(time_vec)]

excess_spreadability_vec_disorder = [NA.get_excess_spreadability(
    disorder_bonds, time_vec[i]; 
    consider_spectral_density = consider_spectral_density) for i in 
    eachindex(time_vec)]

excess_spreadability_vec_disorder_2 = [NA.get_excess_spreadability(
    disorder_2_bonds, time_vec[i]; 
    consider_spectral_density = consider_spectral_density) for i in 
    eachindex(time_vec)]

excess_spreadability_vec_disorder_3 = [NA.get_excess_spreadability(
    disorder_3_bonds, time_vec[i]; 
    consider_spectral_density = consider_spectral_density) for i in 
    eachindex(time_vec)]

# plot the excess spreadability as a function of time in a log-log plot
Plots.plot(time_vec,  Measurements.value.(excess_spreadability_vec_diamond),
    ribbon =  Measurements.uncertainty.(excess_spreadability_vec_diamond),
    xlabel = "t", ylabel = "Excess spreadability", xscale = :log10,
    yscale = :log10, label = Latex.L"kT=0", legend_position = :bottomleft)
Plots.plot!(time_vec,  Measurements.value.(excess_spreadability_vec_disorder),
    ribbon =  Measurements.uncertainty.(excess_spreadability_vec_disorder),
    label = Latex.L"kT=0.1", grid = true)
Plots.plot!(time_vec,  Measurements.value.(excess_spreadability_vec_disorder_2),
    ribbon =  Measurements.uncertainty.(excess_spreadability_vec_disorder_2),
    label = Latex.L"kT=0.16",)
Plots.plot!(time_vec,  Measurements.value.(excess_spreadability_vec_disorder_3),
    ribbon =  Measurements.uncertainty.(excess_spreadability_vec_disorder_3),
    label = Latex.L"kT=0.24",)
Plots.plot!(time_vec, comparison_func.(time_vec, 1), label = Latex.L"t^{-3/2}", color = :black, linestyle = :dash)
Plots.savefig(save_path*"excess_spreadability.png")



function my_func(x, t)
    return x.^2 .* exp.( .- x.^2 .* t)
end

x_vec = collect(0:0.1:10)
Plots.plot(x_vec, my_func.(x_vec, 0.1), label=Latex.L"t=0.1")
Plots.plot!(x_vec, my_func.(x_vec,1), label=Latex.L"t=1")
Plots.plot!(xlabel=Latex.L"k", ylabel=Latex.L"k^2 \mathrm{exp}(-k^2 t)", xlim=(0, 10), ylim=(0, 4))
Plots.savefig(path*"spreadability_func.png")




order_metrics = [
    "Bond length std. deviation",
    "Bond angle std. deviation",
    "Dihedral angle entropy",
    "Bond orientation entropy",
    "Anisotropy from structure factor",
    "Vertex homogeneity metric",
    "Critical pore radius",
    "Ring radius mean",
    "Ring radius std. deviation",
    "Hyperuniformity alpha ",
]

network_predictions = [
    0.18573039770126343,  # bond_length_std_vec
    0.39029574394226074,  # bond_angle_std_vec
    1.013838768005371,     # dihedral_angle_entropy_vec
    0.9829861521720886,   # bond_orientation_entropy_vec
    0.37492454051971436,  # anisotropy_metric_from_structure_factor_vec
    0.7437007427215576,   # vertex_homogeneity_metric_vec
    0.5594326257705688,   # critical_pore_radius_vec
    0.7401111125946045,   # ring_radius_mean_vec
    0.06156785786151886,  # ring_radius_std_vec
    -0.8236029148101807,  # hyperuniformity_alpha_vec_values
]

network_means = [0.203240492418587,  # bond_length_std_vec
                 0.3976747915609694,  # bond_angle_std_vec
                 0.9965101271239464,  # dihedral_angle_entropy_vec
                 0.9861363855433791,  # bond_orientation_entropy_vec
                 0.35689439329461065, # anisotropy_metric_from_structure_factor_vec
                 0.7447058823529412,  # vertex_homogeneity_metric_vec
                 0.48674841804849805, # critical_pore_radius_vec
                 0.7303691876648305,  # ring_radius_mean_vec
                 0.057474754688713,   # ring_radius_std_vec
                 -0.223,              # hyperuniformity_alpha_vec_values
                ]

network_stds = [0.053649381996599285,  # bond_length_std_vec
                0.08306661928092976,  # bond_angle_std_vec
                0.0012159451211699987, # dihedral_angle_entropy_vec
                0.003499055250514032,  # bond_orientation_entropy_vec
                0.01483614201294338,   # anisotropy_metric_from_structure_factor_vec
                0.06556215159621544,   # vertex_homogeneity_metric_vec
                0.1006467141361276,    # critical_pore_radius_vec
                0.043704267868199316,  # ring_radius_mean_vec
                0.017404312256294262,  # ring_radius_std_vec
                0.302,                 # hyperuniformity_alpha_vec_values
                ]

reversed_names = reverse(order_metrics)
reversed_means = reverse(network_means)
reversed_stds = reverse(network_stds)
reversed_predictions = reverse(network_predictions)


Plots.plot(
    [0.98, 0.37, 0.185],
    length(reversed_predictions)-2:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 5,
    markercolor = Plots.palette(:tab10)[3],
    xlims = (-0.9, 1.1),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
)

Plots.savefig(path*"pachy_metrics.png")

Plots.plot(
    reversed_predictions,
    1:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 7,
    markercolor = Plots.palette(:tab10)[2],
    xlims = (-0.9, 1.1),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
)

Plots.plot!(
    [0.98, 0.37, 0.185],
    length(reversed_predictions)-2:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 5,
    markercolor = Plots.palette(:tab10)[3],
    xlims = (-0.9, 1.1),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
)

Plots.savefig(path*"network_predictions.png")

Plots.plot(
    reversed_means,
    1:length(reversed_means),
    xerror = reversed_stds,
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 10,
    xlims = (-0.9, 1.1),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
)

Plots.plot!(
    reversed_predictions,
    1:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 7,
)

Plots.plot!(
    [0.98, 0.37, 0.185],
    length(reversed_predictions)-2:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 5,
    markercolor = Plots.palette(:tab10)[3],
    xlims = (-0.9, 1.1),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
)

Plots.savefig(path*"network_predictions_results.png")

Plots.plot(
    reversed_means,
    1:length(reversed_means),
    xerror = reversed_stds,
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 10,
    xlims = (-0.9, 1.1),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
)

Plots.plot!(
    [0.98, 0.37, 0.185],
    length(reversed_predictions)-2:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 5,
    markercolor = Plots.palette(:tab10)[3],
    xlims = (-0.9, 1.1),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
)

Plots.savefig(path*"weevil_results.png")



path = raw"..\..\presentations\material\\"

x_vec = collect(0:0.01:2)
y_vec = collect(0.5:0.01:1.5)
Plots.plot(x_vec, (3/16) .* (x_vec.^2 .- 1).^2   )
Plots.plot!(xlabel="Bond length / "*Latex.L"d", ylabel="Energy", right_margin = 3Plots.mm, ylims=(0,0.3), xlims=(0,2), legend=false)

Plots.savefig(path*"bond_stretching_energy_3.png")


x_vec = collect(0:0.1:180)
y_vec = collect(40:0.1:180)
Plots.plot(x_vec, (3/8 * 0.2) .* (cosd.(x_vec) .+ 1).^2  )
Plots.plot!(xlabel="Bond angle / °", ylabel="Energy", right_margin = 5Plots.mm, ylims=(0,0.3), xlims=(0,180), legend=false)
Plots.savefig(path*"bond_bending_energy_3.png")


analysis_data_path = raw"..\analysis_data\neural_network_networks\srs\\"

save_path = raw"..\plots\neural_network_networks\srs\\"

all_order_metrics = GU.load_h5_dict(analysis_data_path * "all_order_metrics.h5")

beta_vec = all_order_metrics["bond_bending_const_vec"]

t_max_vec = all_order_metrics["t_max_vec"]
t_gradient_vec = all_order_metrics["t_gradient_vec"]

t_melt_vec = [NA.get_melting_temperature("srs", beta) for beta in beta_vec]

t_max_over_t_melt_vec = t_max_vec ./ t_melt_vec
t_gradient_over_t_melt_vec = t_gradient_vec ./ t_melt_vec

critical_pore_radius_vec = all_order_metrics["critical_pore_radius_vec"]

# Filter the three vectors to only the data with pore radii below 0.45
filter_indices = findall(critical_pore_radius_vec .< 0.45)
t_max_over_t_melt_vec = t_max_over_t_melt_vec[filter_indices]
t_gradient_over_t_melt_vec = t_gradient_over_t_melt_vec[filter_indices]
critical_pore_radius_vec = critical_pore_radius_vec[filter_indices]
beta_vec = beta_vec[filter_indices]

# Create an interactive 3d scatter plot with Makie 
fig = GLMakie.Figure()
ax = GLMakie.Axis3(fig[1, 1], xlabel = Latex.L"T_\mathrm{max} /T_\mathrm{melt}", ylabel = Latex.L"\Delta T /T_\mathrm{melt}", zlabel = Latex.L"\beta")

# Plot the scatter points
#GLMakie.scatter!(ax, t_max_over_t_melt_vec, t_gradient_over_t_melt_vec, critical_pore_radius_vec)
GLMakie.scatter!(ax, t_max_over_t_melt_vec, t_gradient_over_t_melt_vec, beta_vec)

fig

Plots.scatter(
    t_max_over_t_melt_vec, 
    t_gradient_over_t_melt_vec, 
    xlabel = Latex.L"T_\mathrm{max} /T_\mathrm{melt}",
    ylabel = Latex.L"\Delta T /T_\mathrm{melt}",
    label = "Critical pore radius < 0.45",
)
Plots.savefig(save_path*"t_gradient_vs_t_max_critical_pore_radius_below_0.45.png")




path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\neural_network_targeted\test_networks\\"


function load_table_as_dict(filename::String)
    df = CSV.File(filename; delim='\t') |> DataFrames.DataFrame
    return Dict(col => collect(df[!, col]) for col in names(df))
end

analysis_data_path = raw"..\analysis_data\neural_network_targeted\test_networks\\"

measured_means_dict = load_table_as_dict(analysis_data_path*"measured_order_metric_means.txt")
measured_stds_dict = load_table_as_dict(analysis_data_path*"measured_order_metric_stds.txt")
predicted_dict = load_table_as_dict(analysis_data_path*"predicted_order_metrics.txt")


order_metrics_labels = [
    "Bond length std. deviation",
    "Bond angle std. deviation",
    "Dihedral angle entropy",
    "Bond orientation entropy",
    "Anisotropy from structure factor",
    "Vertex homogeneity metric",
    "Critical pore radius",
    "Ring radius mean",
    "Ring radius std. deviation",
    "Hyperuniformity alpha / 10",
]

order_metrics_keys = [
    "bond_length_std_vec",
    "bond_angle_std_vec",
    "dihedral_angle_entropy_vec",
    "bond_orientation_entropy_vec",
    "anisotropy_metric_from_structure_factor_vec",
    "vertex_homogeneity_metric_vec",
    "critical_pore_radius_vec",
    "ring_radius_mean_vec",
    "ring_radius_std_vec",
    "hyperuniformity_alpha_vec_values",
]

reversed_names = reverse(order_metrics_labels)

network_types_vec = predicted_dict["network_type"]

# loop through all network predictions
for i in eachindex(network_types_vec)
    # get the vector of predicitons for the current network type
    predictions = [predicted_dict[key][i] for key in order_metrics_keys]

    # get the vector of measured means for the current network type
    means = [measured_means_dict[key][i] for key in order_metrics_keys]
    # get the vector of measured stds for the current network type
    stds = [measured_stds_dict[key][i] for key in order_metrics_keys]

    # divide the last element of the three vectors by 10
    predictions[end] /= 10
    means[end] /= 10
    stds[end] /= 10

    # calculate the reversed predictions
    reversed_predictions = reverse(predictions)
    reversed_means = reverse(means)
    reversed_stds = reverse(stds)

    # create a label for the current network type
    Plots.plot(
    reversed_predictions,
    1:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 7,
    markercolor = Plots.palette(:tab10)[1],
    xlims = (-0.1, 1.4),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
    )

    Plots.plot!(
    reversed_means,
    1:length(reversed_means),
    xerror = reversed_stds,
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    xticks = collect(0.0:0.5:1.0),
    bottom_margin = 0Plots.mm,
    markersize = 5,
    markercolor = Plots.palette(:tab10)[2],
    xlims = (-0.1, 1.4),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
    )

    Plots.xlabel!("Order metric")

    network_type = network_types_vec[i]
    beta = predicted_dict["bond_bending_const_vec"][i]
    t_max = predicted_dict["t_max_vec"][i]
    t_gradient = predicted_dict["t_gradient_vec"][i]

    Plots.title!("$(network_type), "*Latex.L"\beta ="*"$(beta), "*Latex.L"T_\mathrm{max} ="*"$(t_max), "*Latex.L"\Delta T ="*"$(t_gradient)                         ")

    Plots.savefig(path*"$(network_type)_beta_$(beta)_t_max_$(t_max)_t_gradient_$(t_gradient).png")
end



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\neural_network_targeted\test_networks_pca_5_10_epochs\\"


function load_table_as_dict(filename::String)
    df = CSV.File(filename; delim='\t') |> DataFrames.DataFrame
    return Dict(col => collect(df[!, col]) for col in names(df))
end

analysis_data_path = raw"..\analysis_data\neural_network_targeted\test_networks\\"

measured_means_dict = load_table_as_dict(analysis_data_path*"measured_order_metric_means.txt")
measured_stds_dict = load_table_as_dict(analysis_data_path*"measured_order_metric_stds.txt")
predicted_dict = load_table_as_dict(analysis_data_path*"predicted_order_metrics_pca_5.txt")


order_metrics_labels = [
    "Bond length std. deviation",
    "Bond angle std. deviation",
    "Dihedral angle entropy",
    "Bond orientation entropy",
    "Anisotropy from structure factor",
    "Vertex homogeneity metric",
    "Critical pore radius",
    "Ring radius mean",
    "Ring radius std. deviation",
    "Hyperuniformity alpha / 10",
]

order_metrics_keys = [
    "bond_length_std_vec",
    "bond_angle_std_vec",
    "dihedral_angle_entropy_vec",
    "bond_orientation_entropy_vec",
    "anisotropy_metric_from_structure_factor_vec",
    "vertex_homogeneity_metric_vec",
    "critical_pore_radius_vec",
    "ring_radius_mean_vec",
    "ring_radius_std_vec",
    "hyperuniformity_alpha_vec_values",
]

reversed_names = reverse(order_metrics_labels)

network_types_vec = predicted_dict["network_type"]

# loop through all network predictions
for i in eachindex(network_types_vec)
    # get the vector of predicitons for the current network type
    predictions = [predicted_dict[key][i] for key in order_metrics_keys]

    # get the vector of measured means for the current network type
    means = [measured_means_dict[key][i] for key in order_metrics_keys]
    # get the vector of measured stds for the current network type
    stds = [measured_stds_dict[key][i] for key in order_metrics_keys]

    # divide the last element of the three vectors by 10
    predictions[end] /= 10
    means[end] /= 10
    stds[end] /= 10

    # calculate the reversed predictions
    reversed_predictions = reverse(predictions)
    reversed_means = reverse(means)
    reversed_stds = reverse(stds)

    # create a label for the current network type
    Plots.plot(
    reversed_predictions,
    1:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 7,
    markercolor = Plots.palette(:tab10)[1],
    xlims = (-0.1, 1.4),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
    )

    Plots.plot!(
    reversed_means,
    1:length(reversed_means),
    xerror = reversed_stds,
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    xticks = collect(0.0:0.5:1.0),
    bottom_margin = 0Plots.mm,
    markersize = 5,
    markercolor = Plots.palette(:tab10)[2],
    xlims = (-0.1, 1.4),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
    )

    Plots.xlabel!("Order metric")

    network_type = network_types_vec[i]
    beta = predicted_dict["bond_bending_const_vec"][i]
    t_max = predicted_dict["t_max_vec"][i]
    t_gradient = predicted_dict["t_gradient_vec"][i]

    Plots.title!("$(network_type), "*Latex.L"\beta ="*"$(beta), "*Latex.L"T_\mathrm{max} ="*"$(t_max), "*Latex.L"\Delta T ="*"$(t_gradient)                         ")

    Plots.savefig(path*"$(network_type)_beta_$(beta)_t_max_$(t_max)_t_gradient_$(t_gradient).png")
end


path = path = raw"..\..\presentations\material\\"

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\photonics\others\\"
csv_filename = "cie_2006_2deg_xyz_color_matching_functions.csv"

# Load the data (no header)
df = CSV.read(load_path * csv_filename, DataFrames.DataFrame; header=false)

# Extract columns into separate vectors
wavelength_nm = df[:, 1]             # Column 1: wavelength in nm
f1 = df[:, 2]                        # Column 2
f2 = df[:, 3]                        # Column 3
f3 = df[:, 4]                        # Column 4

# Convert wavelength [nm] to frequency [THz]
# λ [nm] → f [THz] = c / λ
# c = 299792458 m/s = 299792.458 nm/ps = 299792.458 THz·nm
c_THz_nm = 299792.458
frequency_THz = c_THz_nm ./ wavelength_nm

# Plot the three functions vs frequency
plt = Plots.plot(frequency_THz, f1, label="x(ν)", xlabel="Frequency / THz", ylabel="Color Matching Function", c = :red)
Plots.plot!(plt, frequency_THz, f2, label="y(ν)", c = :green)
Plots.plot!(plt, frequency_THz, f3, label="z(ν)", c = :blue)

# Set Y limits to [0, 2]
Plots.ylims!(plt, (0, 2))

# Create a grid
Plots.plot!(plt, grid=true)

# Optionally display the plot (if running in script or REPL)
Plots.savefig(path * "cie_2006_2deg_xyz_color_matching_functions.png")

# Now plot the same data as a function of wavelength
plt2 = Plots.plot(wavelength_nm, f1, label="x(λ)", xlabel="Wavelength / nm", ylabel="Color Matching Function", c = :red)
Plots.plot!(plt2, wavelength_nm, f2, label="y(λ)", c = :green)
Plots.plot!(plt2, wavelength_nm, f3, label="z(λ)", c = :blue)       
# Set Y limits to [0, 2]
Plots.ylims!(plt2, (0, 2))
# Create a grid
Plots.plot!(plt2, grid=true)
# Optionally display the plot (if running in script or REPL)
Plots.savefig(path * "cie_2006_2deg_xyz_color_matching_functions_wavelength.png")


p0_ctn = [0.28, 0.7, 0.25, 0.2, 0.34, 0.05]
p0_dia = [0.28, 0.7, 0.25, 0.2, 0.4, 0.05]
p0_lcs = [0.28, 0.7, 0.25, 0.2, 0.38, 0.05]
p0_srs = [0.28, 0.7, 0.25, 0.2, 0.3, 0.05]


network_type_vec = ["ctn", "dia", "lcs", "srs"]
p0_list = [p0_ctn, p0_dia, p0_lcs, p0_srs]

i = 2
network_type = network_type_vec[i]

upper_bounds::Vector{Float64}=[0.4, 0.8, 0.35, 0.4, 0.45, 0.1]
lower_bounds::Vector{Float64}=[0.2, 0.55, 0.15, 0.0, 0.25, 0.01]


# load the order metrics dict
analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_networks\\" * network_type * raw"\\"
plot_path = raw"..\..\photonics\tidy3d\plots\neural_network_networks\dia\run_1_2_r_t_low_n\\"

# get the r_t dict path
r_t_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\photonics\tidy3d\simulation_data\neural_network_networks\\" * network_type * raw"\run_1_2_r_t_low_n\\"

filename_vec = ["dia_beta_0.1369_t_max_0.1993_t_gradient_0.9613_for_sim", # large peak
                "dia_beta_0.1905_t_max_0.4376_t_gradient_0.6543_for_sim", # no peak
]

p0 = p0_dia

colors = ["#1C5EDF", "#0DAEDB", "#FF7307", "#FFA907"]
labels = ["Network 1 FDTD", "Network 1 fit", "Network 2 FDTD", "Network 2 fit"]

Plots.plot()

#file_nr = 1
for (i, filename) in enumerate(filename_vec)

    r_t_dict = GU.load_h5_dict(r_t_dict_path * filename * "_n_1.5_r_t_only.hdf5")

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
    fitted_reflection = [model(freq, fit_params) for freq in freqs]

    Plots.plot!(freqs, reflection, label=labels[(i-1)*2 + 1], color=colors[(i-1)*2 + 1])
    Plots.plot!(freqs, Measurements.value.(fitted_reflection), 
        ribbon=Measurements.uncertainty.(fitted_reflection), label=labels[(i-1)*2 + 2], 
        color=colors[(i-1)*2 + 2])

end

Plots.plot!(xlabel=Latex.L"Frequency $\omega d / (2 \pi c)$",  ylabel="Reflectance")
Plots.xlims!(0.2, 0.7)
Plots.ylims!(0, 0.65)


Plots.savefig(plot_path*"peak_no_peak_reflection_fit.png")


i = 1
filename = filename_vec[i]

r_t_dict = GU.load_h5_dict(r_t_dict_path * filename * "_n_1.5_r_t_only.hdf5")
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
fitted_reflection = [model(freq, fit_params) for freq in freqs]
Plots.plot!(freqs, reflection, label=labels[(i-1)*2 + 1], color=colors[(i-1)*2 + 1])
Plots.plot!(freqs, Measurements.value.(fitted_reflection), 
    ribbon=Measurements.uncertainty.(fitted_reflection), label=labels[(i-1)*2 + 2], 
    color=colors[(i-1)*2 + 2])
fitted_background = [background_peak(freq, fit_params) for freq in freqs]
# plot only the fitted background peak
Plots.plot!(freqs, Measurements.value.(fitted_background), label="Network 1 background", color="#3D7C8E", linestyle=:dot)

Plots.plot!(xlabel=Latex.L"Frequency $\omega d / (2 \pi c)$",  ylabel="Reflectance")
Plots.xlims!(0.2, 0.7)
Plots.ylims!(0, 0.65)


Plots.savefig(plot_path*"peak_reflection_fit_background.png")


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\neural_network_networks\network_comparison\\"

dicts = [] 
network_type_vec = ["dia", "ctn", "lcs", "srs"]

for (i, network_type) in enumerate(network_type_vec)

    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_networks\\" * network_type * raw"\\"
    pearson_corr_dict = GU.load_h5_dict(analysis_data_path * "order_relative_peak_height_pearson_correlations.h5")
    push!(dicts, pearson_corr_dict)
end

labels = reverse([
"Vertex homogeneity",
Latex.L"Hyperuniformity $\alpha$",
"Ring radius st. d.",
    "Bond length st. d.",
    "Bond angle st. d.",
    "Dihedral angle entropy",
    "Critical pore radius",
    "Bond orientation entropy",
])

key_list = reverse([
"vertex_homogeneity_metric_vec",
"hyperuniformity_alpha_vec",
"ring_radius_std_vec",
"bond_length_std_vec",
"bond_angle_std_vec",
"dihedral_angle_entropy_vec",
"critical_pore_radius_vec",
"bond_orientation_entropy_vec",
])

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Initialize plot
plt = Plots.plot(legend=false, xlabel="Pearson corr.", yticks=(yvals, labels),
size=(600, 600), xticks=[-0.3, 0, 0.3],  grid=true,)

# Add one scatter series per dictionary
for (dict, label) in zip(dicts, labels)
    Plots.scatter!([dict[k] for k in key_list], yvals; label=label, markersize=10)
end

Plots.savefig(save_path * "order_relative_peak_height_pearson_correlations_all_networks.png")


path = path = raw"..\..\presentations\material\\"

x = collect(-0.2:0.01:1.0)

color_low_T = "#1C5EDF"
color_high_T = "#FF7307"

plot_1 = min.(1, exp.( .-  x ./ 0.2  ) )  
plot_2 = min.(1, exp.( .-  x ./ 0.4  ) )  

myplot = Plots.plot(x, plot_2, label = Latex.L"kT=0.4", linecolor=color_high_T)
myplot = Plots.plot!(x, plot_1, label = Latex.L"kT=0.2", linecolor=color_low_T)

Plots.plot!(grid=false, xlabel="Energy difference",
ylabel = "Acceptance prob.", size = (450, 300), xlims=(-0.2,1.0))

Plots.savefig(path*"boltzmann_metropolis_temperatures_0.2_0.4.png")


function heat_cool_temperature_vec(x, max_temp, heat_cool_rate)

    if x < max_temp/heat_cool_rate 
        return heat_cool_rate * x
    elseif x < 2*max_temp/heat_cool_rate 
        return 2* max_temp - heat_cool_rate*x
    else
        return 0
    end
end


x_vec = collect(0:0.01:5)


Plots.plot()

Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.4, 0.7),  alpha=1.0, color=color_high_T)
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.2, 1.0),  alpha=1.0, color=color_low_T)
Plots.plot!(legend = false, xlabel="Attempts per bond chain", ylabel=Latex.L"kT", right_margin = 4Plots.mm, size = (450, 300), xlims=(0,2))

Plots.savefig(path*"heat_cool_temperature_profile_5.png")


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\neural_network_networks\network_comparison\\"

dicts = [] 
network_type_vec = ["dia", "ctn", "lcs", "srs"] #


key_list = [
    "dihedral_angle_entropy_vec",
    "bond_angle_std_vec",
"bond_length_std_vec",
"critical_pore_radius_vec",
"bond_orientation_entropy_vec",
"hyperuniformity_alpha_vec",
"vertex_homogeneity_metric_vec",
"ring_radius_std_vec",
]

labels = [
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    "Critical pore radius",
    "Isotropy",
    Latex.L"Hyperuniformity $\alpha$",
    "Vertex homogeneity",
    Latex.L"Ring radii $\sigma$",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)



for (i, network_type) in enumerate(network_type_vec)

    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_networks\\" * network_type * raw"\\"
    pearson_corr_dict = GU.load_h5_dict(analysis_data_path * "order_relative_peak_height_pearson_correlations.h5")
    push!(dicts, pearson_corr_dict)
end

# Example split indices
n1 = 3
n2 = 7


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:end]
key_lists3 = key_list[n2+1:end]
yvals3 = 1:length(labels3)

# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels, length(labels3)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

# Initialize the three subplots
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = ([-0.3, 0, 0.3], ["", "", ""]),
    grid = true
)

p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = ([-0.3, 0, 0.3], ["", "", ""]),
    grid = true
)

p3 = Plots.plot(
    legend = false,
    yticks = (yvals3, labels3),
    xlabel = "Correlation",
    #ylabel = "Topological",
    xticks = [-0.3, 0, 0.3],
    grid = true
)

# Plot the scatter series for each subplot
for dict in dicts
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=10, ylims=(0.5, n1+0.5), xlims=(-0.35, 0.35))

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=10, ylims=(0.5, length(labels2)+0.5), xlims=(-0.35, 0.35))

    Plots.scatter!(p3, [Measurements.value.(dict[k]) for k in key_lists3], yvals3; label=labels3, markersize=10, ylims=(0.5, length(labels3)+0.5), xlims=(-0.35, 0.35))
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, layout = (3, 1), size = (600, 600), link = :x,
bottom_margin = 3Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 1Plots.mm,
    right_margin = 1Plots.mm)

## Add group labels using annotations
#Plots.annotate!(p1, -0.5, Statistics.mean(yvals1), Plots.text("Local", :center, 12, :bold, :rotation => 90))
#Plots.annotate!(p2, -0.5, Statistics.mean(yvals2), Plots.text("Global", :center, 12, :bold, :rotation => 90))
#Plots.annotate!(p3, -0.5, Statistics.mean(yvals3), Plots.text("Topological", :center, 12, :bold, :rotation => 90))
#

Plots.savefig(save_path * "order_relative_peak_height_pearson_correlations_all_networks_grouped.png")



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\neural_network_networks\network_comparison\\"

dicts = [] 
network_type_vec = ["dia"] #, "ctn", "lcs", "srs"


key_list = [
    "dihedral_angle_entropy_vec",
    "bond_angle_std_vec",
"bond_length_std_vec",
"critical_pore_radius_vec",
"bond_orientation_entropy_vec",
"hyperuniformity_alpha_vec",
"vertex_homogeneity_metric_vec",
"ring_radius_std_vec",
]

labels = [
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    "Critical pore radius",
    "Isotropy",
    Latex.L"Hyperuniformity $\alpha$",
    "Vertex homogeneity",
    Latex.L"Ring radii $\sigma$",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)



for (i, network_type) in enumerate(network_type_vec)

    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_networks\\" * network_type * raw"\\"
    pearson_corr_dict = GU.load_h5_dict(analysis_data_path * "order_relative_peak_height_pearson_correlations.h5")
    push!(dicts, pearson_corr_dict)
end

# Example split indices
n1 = 3
n2 = 7


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:end]
key_lists3 = key_list[n2+1:end]
yvals3 = 1:length(labels3)

# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels, length(labels3)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

# Initialize the three subplots
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = ([-0.3, 0, 0.3], ["", "", ""]),
    grid = true
)

p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = ([-0.3, 0, 0.3], ["", "", ""]),
    grid = true
)

p3 = Plots.plot(
    legend = false,
    yticks = (yvals3, labels3),
    xlabel = "Correlation",
    #ylabel = "Topological",
    xticks = [-0.3, 0, 0.3],
    grid = true
)

# Plot the scatter series for each subplot
for dict in dicts
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=10, ylims=(0.5, n1+0.5), xlims=(-0.35, 0.35))

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=10, ylims=(0.5, length(labels2)+0.5), xlims=(-0.35, 0.35))

    Plots.scatter!(p3, [Measurements.value.(dict[k]) for k in key_lists3], yvals3; label=labels3, markersize=10, ylims=(0.5, length(labels3)+0.5), xlims=(-0.35, 0.35))
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, layout = (3, 1), size = (600, 600), link = :x,
bottom_margin = 3Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 1Plots.mm,
    right_margin = 1Plots.mm)

Plots.savefig(save_path * "order_relative_peak_height_pearson_correlations_dia_grouped.png")


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\neural_network_targeted\test_networks\\"


function load_table_as_dict(filename::String)
    df = CSV.File(filename; delim='\t') |> DataFrames.DataFrame
    return Dict(col => collect(df[!, col]) for col in names(df))
end

function reorder_dict(dict1::Dict, dict2::Dict, keys_for_order::Vector; digits::Int=4)
    # Helper: round floats, leave others unchanged
    process(x) = x isa AbstractFloat ? round(x; digits=digits) : x

    # Build the reference ordering from dict1
    tuples1 = collect(zip((process.(dict1[k]) for k in keys_for_order)...))

    # Build tuples from dict2
    tuples2 = collect(zip((process.(dict2[k]) for k in keys_for_order)...))

    # Find the mapping that sorts tuples2 such that it equals tuples1
    perm = indexin(tuples1, tuples2)
    reordered = Dict()
    for (k, v) in dict2
        reordered[k] = v[perm]
    end

    return reordered
end

analysis_data_path = raw"..\analysis_data\neural_network_targeted\test_networks\\"

measured_means_dict = load_table_as_dict(analysis_data_path*"measured_order_metric_means.txt")
measured_stds_dict = load_table_as_dict(analysis_data_path*"measured_order_metric_stds.txt")
predicted_dict = load_table_as_dict(analysis_data_path*"predicted_order_metrics_10_epochs_no_regularization.txt")

predicted_dict = reorder_dict(measured_means_dict, predicted_dict, ["network_type", "bond_bending_const_vec", "t_max_vec", "t_gradient_vec"])

order_metrics_labels = [
    "Bond length std. deviation",
    "Bond angle std. deviation",
    "Dihedral angle entropy",
    "Bond orientation entropy",
    "Anisotropy from structure factor",
    "Vertex homogeneity metric",
    "Critical pore radius",
    "Ring radius mean",
    "Ring radius std. deviation",
    "Hyperuniformity alpha / 10",
]

order_metrics_keys = [
    "bond_length_std_vec",
    "bond_angle_std_vec",
    "dihedral_angle_entropy_vec",
    "bond_orientation_entropy_vec",
    "anisotropy_metric_from_structure_factor_vec",
    "vertex_homogeneity_metric_vec",
    "critical_pore_radius_vec",
    "ring_radius_mean_vec",
    "ring_radius_std_vec",
    "hyperuniformity_alpha_vec_values",
]

reversed_names = reverse(order_metrics_labels)

network_types_vec = predicted_dict["network_type"]

# loop through all network predictions
for i in eachindex(network_types_vec)
    # get the vector of predicitons for the current network type
    predictions = [predicted_dict[key][i] for key in order_metrics_keys]

    # get the vector of measured means for the current network type
    means = [measured_means_dict[key][i] for key in order_metrics_keys]
    # get the vector of measured stds for the current network type
    stds = [measured_stds_dict[key][i] for key in order_metrics_keys]

    # divide the last element of the three vectors by 10
    predictions[end] /= 10
    means[end] /= 10
    stds[end] /= 10

    # calculate the reversed predictions
    reversed_predictions = reverse(predictions)
    reversed_means = reverse(means)
    reversed_stds = reverse(stds)

    # create a label for the current network type
    Plots.plot(
    reversed_predictions,
    1:length(reversed_predictions),
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    bottom_margin = 0Plots.mm,
    markersize = 7,
    markercolor = Plots.palette(:tab10)[1],
    xlims = (-0.1, 1.4),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
    )

    Plots.plot!(
    reversed_means,
    1:length(reversed_means),
    xerror = reversed_stds,
    seriestype = :scatter,
    markershape = :circle,
    legend = false,
    size = (700, 600),
    yticks = (1:length(reversed_names), reversed_names),
    xticks = collect(0.0:0.5:1.0),
    bottom_margin = 0Plots.mm,
    markersize = 5,
    markercolor = Plots.palette(:tab10)[2],
    xlims = (-0.1, 1.4),
    ylims = (0.5, length(reversed_names) + 0.5),
    grid = true,
    )

    Plots.xlabel!("Order metric")

    network_type = network_types_vec[i]
    beta = predicted_dict["bond_bending_const_vec"][i]
    t_max = predicted_dict["t_max_vec"][i]
    t_gradient = predicted_dict["t_gradient_vec"][i]

    Plots.title!("$(network_type), "*Latex.L"\beta ="*"$(beta), "*Latex.L"T_\mathrm{max} ="*"$(t_max), "*Latex.L"\Delta T ="*"$(t_gradient)                         ")

    Plots.savefig(path*"$(network_type)_beta_$(beta)_t_max_$(t_max)_t_gradient_$(t_gradient).png")
end


function load_table_as_dict(filename::String)
    df = CSV.File(filename; delim='\t') |> DataFrames.DataFrame
    return Dict(col => collect(df[!, col]) for col in names(df))
end


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\neural_network_targeted\comparison_t_melt\\"

analysis_data_path_1 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_targeted\test_networks\\"

analysis_data_path_2 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_targeted\dia\1000_vertices\\"

filename_1 = "dia_beta_0.2500_t_max_0.4592_t_gradient_0.1224"

filename_2 = "dia_beta_0.2500_t_max_2.8443_t_gradient_1.5169"

measured_means_dict_1 = load_table_as_dict(analysis_data_path_1*"measured_order_metric_means.txt")
measured_stds_dict_1 = load_table_as_dict(analysis_data_path_1*"measured_order_metric_stds.txt")

measured_means_dict_2 = load_table_as_dict(analysis_data_path_2*"measured_order_metric_means.txt")
measured_stds_dict_2 = load_table_as_dict(analysis_data_path_2*"measured_order_metric_stds.txt")

order_metrics_labels = [
    "Bond length std. deviation",
    "Bond angle std. deviation",
    "Dihedral angle entropy",
    "Bond orientation entropy",
    "Anisotropy from structure factor",
    "Vertex homogeneity metric",
    "Critical pore radius",
    "Ring radius mean",
    "Ring radius std. deviation",
    "Hyperuniformity alpha / 10",
]

order_metrics_keys = [
    "bond_length_std_vec",
    "bond_angle_std_vec",
    "dihedral_angle_entropy_vec",
    "bond_orientation_entropy_vec",
    "anisotropy_metric_from_structure_factor_vec",
    "vertex_homogeneity_metric_vec",
    "critical_pore_radius_vec",
    "ring_radius_mean_vec",
    "ring_radius_std_vec",
    "hyperuniformity_alpha_vec_values",
]

reversed_names = reverse(order_metrics_labels)

# in the measured_means_dict_1 find the entry of the vectors where the network_type is dia
# bond_bending_const_vec is 0.25, t_max_vec is 0.4592, t_gradient_vec is 0.1224
index_means = 0 

for i in 1:length(measured_means_dict_1["network_type"])

    if measured_means_dict_1["network_type"][i] == "dia" &&
        isapprox(measured_means_dict_1["bond_bending_const_vec"][i], 0.25; atol=1e-3) &&
        isapprox(measured_means_dict_1["t_max_vec"][i], 0.4592; atol=1e-3) &&
        isapprox(measured_means_dict_1["t_gradient_vec"][i], 0.1224; atol=1e-3)
        global index_means = i
        println("Found index_means: $index_means")
        break
    end
end

index_stds = 0

for j in 1:length(measured_stds_dict_1["network_type"])
    if measured_stds_dict_1["network_type"][j] == "dia" &&
        isapprox(measured_stds_dict_1["bond_bending_const_vec"][j], 0.25; atol=1e-3) &&
        isapprox(measured_stds_dict_1["t_max_vec"][j], 0.4592; atol=1e-3) &&
        isapprox(measured_stds_dict_1["t_gradient_vec"][j], 0.1224; atol=1e-3)
        global index_stds = j
        println("Found index_stds: $index_stds")    
        break
    end
end

# get the vector of measured means for the current network type
means_1 = [measured_means_dict_1[key][index_means] for key in order_metrics_keys]
# get the vector of measured stds for the current network type
stds_1 = [measured_stds_dict_1[key][index_stds] for key in order_metrics_keys]

# get the vector of measured means for the current network type
means_2 = [measured_means_dict_2[key][1] for key in order_metrics_keys]
# get the vector of measured stds for the current network type
stds_2 = [measured_stds_dict_2[key][1] for key in order_metrics_keys]

# divide the last element of the three vectors by 10
means_1[end] /= 10
stds_1[end] /= 10
means_2[end] /= 10
stds_2[end] /= 10

# calculate the reversed means and stds
reversed_means_1 = reverse(means_1)
reversed_stds_1 = reverse(stds_1)
reversed_means_2 = reverse(means_2)
reversed_stds_2 = reverse(stds_2)

# create a label for the current network type
Plots.plot(
reversed_means_1,
1:length(reversed_means_1),
xerror = reversed_stds_1,
seriestype = :scatter,
markershape = :circle,
legend = false,
size = (700, 600),
yticks = (1:length(reversed_names), reversed_names),
xticks = collect(0.0:0.5:1.0),
bottom_margin = 0Plots.mm,
markersize = 5,
markercolor = Plots.palette(:tab10)[1],
xlims = (-0.1, 1.4),
ylims = (0.5, length(reversed_names) + 0.5),
grid = true,
)
Plots.plot!(
reversed_means_2,
1:length(reversed_means_2),
xerror = reversed_stds_2,
seriestype = :scatter,
markershape = :circle,
legend = false,
size = (700, 600),
yticks = (1:length(reversed_names), reversed_names),
xticks = collect(0.0:0.5:1.0),
bottom_margin = 0Plots.mm,
markersize = 5,
markercolor = Plots.palette(:tab10)[2],
xlims = (-0.1, 1.4),
ylims = (0.5, length(reversed_names) + 0.5),
grid = true,
)
Plots.xlabel!("Order metric")

Plots.savefig(save_path*"dia_comparison_t_melt_216_1000.png")



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\comparison_t_melt\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\targeted\shell_nr_4\\"

measured_means_dict = load_table_as_dict(analysis_data_path*"measured_order_metric_means.txt")
measured_stds_dict = load_table_as_dict(analysis_data_path*"measured_order_metric_stds.txt")

order_metrics_labels = [
    "Bond length std. deviation",
    "Bond angle std. deviation",
    "Dihedral angle entropy",
    "Bond orientation entropy",
    "Anisotropy from structure factor",
    "Vertex homogeneity metric",
    "Critical pore radius",
    "Ring radius mean",
    "Ring radius std. deviation",
    "Hyperuniformity alpha / 10",
]

order_metrics_keys = [
    "bond_length_std_vec",
    "bond_angle_std_vec",
    "dihedral_angle_entropy_vec",
    "bond_orientation_entropy_vec",
    "anisotropy_metric_from_structure_factor_vec",
    "vertex_homogeneity_metric_vec",
    "critical_pore_radius_vec",
    "ring_radius_mean_vec",
    "ring_radius_std_vec",
    "hyperuniformity_alpha_vec_values",
]

reversed_names = reverse(order_metrics_labels)

# in the measured_means_dict find the entry of the vectors where the network_type is dia
# bond_bending_const_vec is 0.25, t_max_vec is 0.4592, t_gradient_vec is 0.1224
network_type_vec = ["ctn", "dia", "lcs", "srs"]
nr_vertices_small = [224, 216, 192, 216]
nr_vertices_large = [756, 1000, 648, 1000]

means_identifier = zip(measured_means_dict["network_type"],
                        measured_means_dict["bond_bending_const_vec"],
                        measured_means_dict["t_max_vec"],
                        measured_means_dict["t_gradient_vec"])

stds_identifier = zip(measured_stds_dict["network_type"],
                        measured_stds_dict["bond_bending_const_vec"],
                        measured_stds_dict["t_max_vec"],
                        measured_stds_dict["t_gradient_vec"])

for (i, mean_id) in enumerate(means_identifier)

    j = 0

    if measured_means_dict["nr_vertices_vec"][i] < 300
        
        # Find the first j that matches the same tuple and has measured_means_dict["nr_vertices_vec"][j] > 600
        j = findfirst(j -> 
            measured_stds_dict["network_type"][j] == mean_id[1]
            && isapprox(measured_stds_dict["bond_bending_const_vec"][j], mean_id[2]) 
            && isapprox(measured_stds_dict["t_max_vec"][j], mean_id[3]) 
            && isapprox(measured_stds_dict["t_gradient_vec"][j], mean_id[4]) 
            && measured_means_dict["nr_vertices_vec"][j] > 600, 1:length(measured_stds_dict["network_type"]))

        # get the vector of measured means for the current network type
        means_small = [measured_means_dict[key][i] for key in order_metrics_keys]
        means_large = [measured_means_dict[key][j] for key in order_metrics_keys]
        stds_small = [measured_stds_dict[key][i] for key in order_metrics_keys]
        stds_large = [measured_stds_dict[key][j] for key in order_metrics_keys]

        means_small[end] /= 10
        means_large[end] /= 10
        stds_small[end] /= 10
        stds_large[end] /= 10

        reversed_means_small = reverse(means_small)
        reversed_means_large = reverse(means_large)
        reversed_stds_small = reverse(stds_small)
        reversed_stds_large = reverse(stds_large)

        Plots.plot(
        reversed_means_large,
        1:length(reversed_means_large),
        xerror = reversed_stds_large,
        seriestype = :scatter,
        markershape = :circle,
        legend = false,
        size = (700, 600),
        yticks = (1:length(reversed_names), reversed_names),
        xticks = collect(0.0:0.5:1.0),
        bottom_margin = 0Plots.mm,
        markersize = 7,
        markercolor = Plots.palette(:tab10)[1],
        xlims = (-0.1, 1.4),
        ylims = (0.5, length(reversed_names) + 0.5),
        grid = true,
        )

        Plots.plot!(
        reversed_means_small,
        1:length(reversed_means_small),
        xerror = reversed_stds_small,
        seriestype = :scatter,
        markershape = :circle,
        markersize = 5,
        markercolor = Plots.palette(:tab10)[2],
        )

        Plots.xlabel!("Order metric")

        network_type = mean_id[1]
        beta = mean_id[2]
        t_max = mean_id[3]
        t_gradient = mean_id[4]

        # round these values to 4 decimal places
        beta = round(beta, digits=4)
        t_max = round(t_max, digits=4)
        t_gradient = round(t_gradient, digits=4)

        Plots.title!("$(network_type), "*Latex.L"\beta ="*"$(beta), "*Latex.L"T_\mathrm{max} ="*"$(t_max), "*Latex.L"\Delta T ="*"$(t_gradient)                         ")

        Plots.savefig(save_path*"$(network_type)_beta_$(beta)_t_max_$(t_max)_t_gradient_$(t_gradient).png")
    end
end



save_filename = "dia_64_vertices"
network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\crystals\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\\"
digital_sphere_mask_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\digital_sphere_masks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\crystals\pore_size_distribution_comparison\\"

spatial_network = NG.load_spatial_network_from_gml(
    network_path*save_filename*".gml")

println(spatial_network[]["supercell_edge_length"]/2)

pore_size_distribution_dict = NA.get_pore_size_distribution(
    spatial_network;
    sampling_grid_size = 0.2,
    max_pore_radius = 1.5,
    periodic_boundary_conditions = true,
    save_result = true,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    digital_sphere_mask_path 
        = digital_sphere_mask_path,
    print_progress = true,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

# convert to a spatial_network_without periodic boundaries
spatial_network_no_pbc = NA.convert_periodic_to_non_periodic(spatial_network)

pore_size_distribution_dict_no_pbc = NA.get_pore_size_distribution(
    spatial_network_no_pbc;
    sampling_grid_size = 0.2,
    max_pore_radius = 1.5,
    periodic_boundary_conditions = false,
    save_result = true,
    save_path = analysis_data_path*save_filename*"_no_pbc",
    label = nothing,
    digital_sphere_mask_path 
        = digital_sphere_mask_path,
    print_progress = true,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

Plots.plot(pore_size_distribution_dict["pore_size_vec"], pore_size_distribution_dict["pore_size_distribution"], label="with PBC")
Plots.plot!(pore_size_distribution_dict_no_pbc["pore_size_vec"], pore_size_distribution_dict_no_pbc["pore_size_distribution"], label="without PBC")
Plots.xlabel!("Pore radius / d")
Plots.ylabel!("Pore size distribution")
Plots.savefig(plot_path*save_filename*"_pore_size_distribution_comparison.png")


network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\targeted\shell_nr_4\run_1\\"

save_filename = "lcs_nr_vertices_192_beta_0.2500_t_max_0.2138_t_gradient_0.1140"

analysis_data_path_1 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\targeted\shell_nr_4\run_1\\"
analysis_data_path_2 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\targeted\new_pore_size_distribution\run_1\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\pore_size_distribution_comparison\\"

dict_1 = GU.load_h5_dict(analysis_data_path_1*save_filename*"_pore_size_distribution.h5")
dict_2 = GU.load_h5_dict(analysis_data_path_2*save_filename*"_pore_size_distribution.h5")

Plots.plot(dict_1["pore_size_vec"], dict_1["pore_size_distribution"], label="with error")
Plots.plot!(dict_2["pore_size_vec"], dict_2["pore_size_distribution"], label="corrected")
Plots.xlabel!("Pore radius / d")
Plots.ylabel!("Pore size distribution")
Plots.savefig(plot_path*save_filename*"_with_and_without_error.png")


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"
filename_1 = "pachy_blue_pore_size_distribution.h5"
filename_2 = "pachy_red_pore_size_distribution.h5"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\\"

dict_1 = GU.load_h5_dict(analysis_data_path*filename_1)
dict_2 = GU.load_h5_dict(analysis_data_path*filename_2)

Plots.plot(dict_1["pore_size_vec"], dict_1["pore_size_distribution"], label="Blue")
Plots.plot!(dict_2["pore_size_vec"], dict_2["pore_size_distribution"], label="Red")
Plots.xlabel!("Pore radius / d")
Plots.ylabel!("Pore size distribution")
Plots.savefig(plot_path*"pachy_blue_red_pore_size_distribution.png")


save_filename = "dia_64_vertices"
network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\crystals\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\\"
digital_sphere_mask_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\digital_sphere_masks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\crystals\structure_factor_comparison\\"

maximal_wavevector_int = 5

spatial_network = NG.load_spatial_network_from_gml(
    network_path*save_filename*".gml")

println(spatial_network[]["supercell_edge_length"]/2)

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(
    spatial_network;
    consider_bonds = false,
    maximal_wavevector_int = maximal_wavevector_int,
    periodic_boundary_conditions = true,
    wavevector_array_positive_z = NA.get_wavevector_array_positive_z(spatial_network; 
        maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=periodic_boundary_conditions),
    save_result = false,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    print_progress = true,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(
        structure_factor_dict;
        consider_bonds = false,
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = false,
        save_path = analysis_data_path*save_filename,
        label = nothing)

# convert to a spatial_network_without periodic boundaries
spatial_network_no_pbc = NA.convert_periodic_to_non_periodic(spatial_network)

structure_factor_no_pbc_dict = NA.get_structure_factor_by_wavevector_array(
    spatial_network_no_pbc;
    consider_bonds = false,
    maximal_wavevector_int = maximal_wavevector_int,
    periodic_boundary_conditions = false,
    wavevector_array_positive_z = NA.get_wavevector_array_positive_z(spatial_network_no_pbc; 
        maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=periodic_boundary_conditions),
    save_result = false,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    print_progress = true,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

structure_factor_no_pbc_angle_averaged_dict = NA.get_structure_factor_angle_averaged(
        structure_factor_no_pbc_dict;
        consider_bonds = false,
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = false,
        save_path = analysis_data_path*save_filename,
        label = nothing)

Plots.plot(structure_factor_angle_averaged_dict["wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]), ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict["structure_factor_vec"]), label="with PBC")
Plots.plot!(structure_factor_no_pbc_angle_averaged_dict["wavenumber_vec"], Measurements.value.(structure_factor_no_pbc_angle_averaged_dict["structure_factor_vec"]), ribbon =  Measurements.uncertainty.(structure_factor_no_pbc_angle_averaged_dict["structure_factor_vec"]), label="without PBC")
Plots.xlabel!("Wavenumber / (1/d)")
Plots.ylabel!("Structure factor")
Plots.savefig(plot_path*save_filename*"_structure_factor_comparison.png")

NA.plot_structure_factor_heatmap(
    structure_factor_no_pbc_dict,
    plot_path*save_filename;
    title="No PBC",
    save_plot = false,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)



save_filename = "dia_nr_vertices_216_beta_0.2500_t_max_0.1787_t_gradient_0.0477"
network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\targeted\shell_nr_4\run_1\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\\"
digital_sphere_mask_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\digital_sphere_masks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\crystals\structure_factor_comparison\\"

maximal_wavevector_int = 5

spatial_network = NG.load_spatial_network_from_gml(
    network_path*save_filename*".gml")

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(
    spatial_network;
    consider_bonds = true,
    maximal_wavevector_int = maximal_wavevector_int,
    periodic_boundary_conditions = true,
    wavevector_array_positive_z = NA.get_wavevector_array_positive_z(spatial_network; 
        maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=periodic_boundary_conditions),
    save_result = false,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    print_progress = true,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(
        structure_factor_dict;
        consider_bonds = true,
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = false,
        save_path = analysis_data_path*save_filename,
        label = nothing)

# convert to a spatial_network_without periodic boundaries
spatial_network_no_pbc = NA.convert_periodic_to_non_periodic(spatial_network)

structure_factor_no_pbc_dict = NA.get_structure_factor_by_wavevector_array(
    spatial_network_no_pbc;
    consider_bonds = true,
    maximal_wavevector_int = maximal_wavevector_int,
    periodic_boundary_conditions = false,
    wavevector_array_positive_z = NA.get_wavevector_array_positive_z(spatial_network_no_pbc; 
        maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=periodic_boundary_conditions),
    save_result = false,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    print_progress = true,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

structure_factor_no_pbc_angle_averaged_dict = NA.get_structure_factor_angle_averaged(
        structure_factor_no_pbc_dict;
        consider_bonds = true,
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = false,
        save_path = analysis_data_path*save_filename,
        label = nothing)

Plots.plot(structure_factor_angle_averaged_dict["wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]), ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict["structure_factor_vec"]), label="with PBC")
Plots.plot!(structure_factor_no_pbc_angle_averaged_dict["wavenumber_vec"], Measurements.value.(structure_factor_no_pbc_angle_averaged_dict["structure_factor_vec"]), ribbon =  Measurements.uncertainty.(structure_factor_no_pbc_angle_averaged_dict["structure_factor_vec"]), label="without PBC")
Plots.xlabel!("Wavenumber / (1/d)")
Plots.ylabel!("Structure factor")
Plots.savefig(plot_path*save_filename*"_structure_factor_comparison.png")



analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\\"

filename = "pachy_blue"

structure_factor_blue = GU.load_h5_dict(
    analysis_data_path*filename*"_structure_factor_array.h5")

upper_clim_bonds = 3.0

NA.plot_structure_factor_heatmap(
    structure_factor_blue,
    plot_path*filename;
    title="Pachy blue",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)

filename = "pachy_red"

structure_factor_blue = GU.load_h5_dict(
    analysis_data_path*filename*"_structure_factor_array.h5")

upper_clim_bonds = 3.0

NA.plot_structure_factor_heatmap(
    structure_factor_blue,
    plot_path*filename;
    title="Pachy blue",
    save_plot = true,
    clims = (0, upper_clim_bonds ),
    x_y_lims = nothing)


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\\"

filename_blue = "pachy_blue"

structure_factor_blue = GU.load_h5_dict(
    analysis_data_path*filename_blue*"_structure_factor_angle_averaged.h5")


filename_red = "pachy_red"

structure_factor_red = GU.load_h5_dict(
    analysis_data_path*filename_red*"_structure_factor_angle_averaged.h5")

Plots.plot(structure_factor_blue["wavenumber_vec"], Measurements.value.(structure_factor_blue["structure_factor_vec"]),
ribbon =  Measurements.uncertainty.(structure_factor_blue["structure_factor_vec"]),
    label = "Blue"
)
Plots.plot!(structure_factor_red["wavenumber_vec"], Measurements.value.(structure_factor_red["structure_factor_vec"]),
    ribbon =  Measurements.uncertainty.(structure_factor_red["structure_factor_vec"]),
    label = "Red"
)

Plots.xlabel!("Wavenumber / (1/d)")
Plots.ylabel!("Structure factor")
Plots.savefig(plot_path*"pachy_structure_factor_comparison.png")

structure_factor_red = GU.load_h5_dict(
    analysis_data_path*filename_red*"_structure_factor_angle_averaged.h5")

Plots.plot(structure_factor_blue["wavenumber_vec"], Measurements.value.(structure_factor_blue["structure_factor_vec"]),
ribbon =  Measurements.uncertainty.(structure_factor_blue["structure_factor_vec"]),
    label = "Blue"
)
Plots.plot!(structure_factor_red["wavenumber_vec"], Measurements.value.(structure_factor_red["structure_factor_vec"]),
    ribbon =  Measurements.uncertainty.(structure_factor_red["structure_factor_vec"]),
    label = "Red"
)

Plots.plot!(right_margin = 3Plots.mm)
Plots.xlims!(0, 20)
Plots.ylims!(0, 3)
Plots.xlabel!("Wavenumber / (1/d)")
Plots.ylabel!("Structure factor")
Plots.savefig(plot_path*"pachy_structure_factor_comparison_zoom.png")


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\apodization_test\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\apodization_test\\"

filename = "pachy_red"

apodizations = ["gauss_0.2", "gauss_0.3", "gauss_0.5", "hamming", "hann", "tukey_0.25", "tukey_0.5", "tukey_0.75"]

for apodization in apodizations
    structure_factor = GU.load_h5_dict(
        analysis_data_path*filename*"_"*apodization*"_structure_factor_array.h5")

    upper_clim = 6.0

    NA.plot_structure_factor_heatmap(
        structure_factor,
        plot_path*filename*"_large_upper_clim_"*apodization;
        title="Pachy red, "*apodization,
        save_plot = true,
        clims = (0, upper_clim ),
        x_y_lims = nothing,
        clims_from_mean = true,)
end


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

dict_names = ["pachy_red_structure_factor_array", "pachy_red_structure_factor_bonds_array",
"pachy_blue_structure_factor_array",
"pachy_blue_structure_factor_bonds_array"]


for dict_name in dict_names
    structure_factor = GU.load_h5_dict(
        analysis_data_path*dict_name*".h5")

    upper_clim = 6.0

    NA.plot_structure_factor_heatmap(
        structure_factor,
        plot_path*dict_name*"_large_upper_clim";
        title=dict_name,
        save_plot = true,
        clims = (0, upper_clim ),
        x_y_lims = nothing,
        clims_from_mean = false,)
end



analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\dia\run_1\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\structure_factor_comparison\\"

dict_names = ["dia_beta_0.0852_t_max_0.1036_t_gradient_0.0142_structure_factor_angle_averaged", "dia_beta_0.0852_t_max_0.1036_t_gradient_0.0142_structure_factor_bonds_angle_averaged"]

labels = ["vertices", "bonds"]

Plots.plot(legend = true)

for (dict_name, label) in zip(dict_names, labels)
    structure_factor_angle_averaged_dict = GU.load_h5_dict(
        analysis_data_path*dict_name*".h5")

    Plots.plot!(structure_factor_angle_averaged_dict["wavenumber_vec"], 
                    Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]) , 
                    ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict["structure_factor_vec"]),
                    label=label)
end

Plots.xlabel!("Wavenumber / (1/d)")
Plots.ylabel!("Structure factor")

Plots.xlims!(0, 20)
Plots.ylims!(0, 10)

# add a verticle line at x=pi/2
Plots.vline!([pi/2], linestyle=:dash, color=:black)

Plots.savefig(plot_path*"dia_beta_0.0852_t_max_0.1036_t_gradient_0.0142_structure_factor_angle_averaged.png")


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

dict_names = ["pachy_red_structure_factor_angle_averaged", "pachy_red_structure_factor_bonds_angle_averaged",
"pachy_blue_structure_factor_angle_averaged",
"pachy_blue_structure_factor_bonds_angle_averaged"]

labels = ["pachy red no bonds", "pachy red with bonds", "pachy blue no bonds", "pachy blue with bonds"]

Plots.plot(legend = true)


for (dict_name, label) in zip(dict_names, labels)
    structure_factor_angle_averaged_dict = GU.load_h5_dict(
        analysis_data_path*dict_name*".h5")

    Plots.plot!(structure_factor_angle_averaged_dict["wavenumber_vec"], 
                    Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]), 
                    ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict["structure_factor_vec"]),
                    label=label)
end

Plots.xlabel!("Wavenumber / (1/d)")
Plots.ylabel!("Structure factor")

Plots.xlims!(0, 20)
Plots.ylims!(0, 5)

# add a verticle line at x=pi/2
Plots.vline!([pi/2], linestyle=:dash, color=:black, label=Latex.L"\pi/2")

# add a horizontal line at y=1
Plots.hline!([1], linestyle=:dash, color=:black, label="1")

Plots.savefig(plot_path*"pachy_structure_factor_angle_averaged.png")



analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\predictions\local_relaxation\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

filename = "ctn_predictions_nr_layers_3_nr_neurons_57_full_pca_5.h5"

data_dict = GU.load_h5_dict(analysis_data_path*filename)


predictions_array = data_dict["predictions_array"]
loss_array = data_dict["loss_array"]

# permute dims of the arrays to have the shape (bond_bending_const, t_max, t_gradient )
predictions_array = permutedims(predictions_array, (4, 3, 2, 1))
loss_array = permutedims(loss_array, (3, 2, 1))

bond_bending_const_vec = data_dict["bond_bending_const_vec"]
t_max_vec = data_dict["t_max_vec"]
t_gradient_vec = data_dict["t_gradient_vec"]

# get the index, where the predictions array has the minimal value 
min_positive_value = minimum(loss_array)
min_positive_index = findfirst(x -> x == min_positive_value, loss_array)

# print the values of the parameters at this index
println("Minimum positive loss: $min_positive_value")
println("At bond_bending_const = $(bond_bending_const_vec[min_positive_index[1]])")
println("At t_max = $(t_max_vec[min_positive_index[2]])")
println("At t_gradient = $(t_gradient_vec[min_positive_index[3]])")

max_positive_value = maximum(loss_array[loss_array .< 999.0])
max_positive_index = findfirst(x -> x == max_positive_value, loss_array)

# print the values of the parameters at this index
println("Maximum loss: $max_positive_value")
println("At bond_bending_const = $(bond_bending_const_vec[max_positive_index[3]])")
println("At t_max = $(t_max_vec[max_positive_index[2]])")
println("At t_gradient = $(t_gradient_vec[max_positive_index[1]])")

# plot a heatmap of the predictions for fixed bond_bending_const = 5.0 with a
# logarithmic color scale
fixed_bond_bending_const = bond_bending_const_vec[min_positive_index[1]]
bond_bending_const_index = min_positive_index[1]

# First fixed part
part1 = [
    "bond_length_std",
    "bond_angle_std",
    "dihedral_angle_entropy",
    "bond_orientation_entropy",
    "coordination_nr_mean",
    "coordination_nr_std"
]

# q_l_value_0 ... q_l_value_12
part2 = ["q_l_value_$(i)" for i in 0:12]

# q_l_uncertainty_0 ... q_l_uncertainty_12
part3 = ["q_l_uncertainty_$(i)" for i in 0:12]

# Last fixed part
part4 = [
    "vertex_homogeneity_metric",
    "ring_size_mean",
    "ring_size_std",
    "ring_radius_mean",
    "ring_radius_std",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor",
    "anisotropy_metric_from_structure_factor_bonds",
    "hyperuniformity_alpha_value",
    "hyperuniformity_alpha_uncertainty"
]

# Combine all parts
order_metrics = vcat(part1, part2, part3, part4)

for (index, metric) in enumerate(order_metrics)
    println("Order metric: $metric ", predictions_array[min_positive_index[1], min_positive_index[2], min_positive_index[3], index])
end


Z  = loss_array[bond_bending_const_index, :, :]'
Zlog = log10.(Z)                               
exps = -4:1                                     
ticks_vals = collect(exps)                       
ticks_lbls = ["1e$(p)" for p in exps]

Plots.heatmap(t_max_vec, t_gradient_vec, Zlog;
    xlabel=Latex.L"T_\mathrm{max}",
    ylabel=Latex.L"\Delta T",
    title=Latex.L"Loss at $\beta = $"*string(fixed_bond_bending_const),
    colorbar_title=Latex.L"L",
    c=:viridis,
    clim=(minimum(exps), maximum(exps)),        
    colorbar_ticks=(ticks_vals, ticks_lbls),   
    aspect_ratio=:equal,
)

# set xlims and ylims
Plots.heatmap!(; xlims=(minimum(t_max_vec), maximum(t_max_vec)),
                      ylims=(minimum(t_gradient_vec), maximum(t_gradient_vec)))
Plots.savefig(plot_path*"ctn_loss_heatmap_fixed_bond_bending_const_$(fixed_bond_bending_const).png")


plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

key_list = [
    "ring_radius_std",
    "ring_radius_mean",
    "dihedral_angle_entropy",
    "bond_angle_std",
"bond_length_std",
"anisotropy_metric_from_structure_factor_bonds",
"bond_orientation_entropy",
"hyperuniformity_alpha",
"critical_pore_radius",
"vertex_homogeneity_metric",
]

labels = [
    Latex.L"Ring radii $\sigma$",
    "Mean ring radius",
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    "Anisotropy in str. fact.",
    "Isotropy in bonds",
    Latex.L"Hyperuniformity $\alpha / 10$",
    "Critical pore radius",
    "Vertex homogeneity",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"

filenames = ["pachy_blue_order_metrics.h5", "pachy_red_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_path * filename)
    push!(dicts, order_metrics_dict)
end

# Example split indices
n1 = 5


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1
labels2 = labels[n1+1:end]
key_lists2 = key_list[n1+1:end]
yvals2 = 1:length(labels2)


# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

xlims = (-0.1, 1.1)
xticks = [0, 0.5, 1.0]

# Initialize the three subplots
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = xticks,
    xlabel = "Order metric value",
    grid = true
)

# divide the value of hyperuniformity_alpha by 10
for dict in dicts
    dict["hyperuniformity_alpha"] = dict["hyperuniformity_alpha"] / 10
end

# Plot the scatter series for each subplot
for dict in dicts
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=10, ylims=(0.5, n1+0.5), xlims=xlims)

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=10, ylims=(0.5, length(labels2)+0.5), xlims=xlims)
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, layout = (2, 1), size = (600, 600), link = :x,
bottom_margin = 3Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 1Plots.mm,
    right_margin = 1Plots.mm)

Plots.savefig(plot_path * "pachy_order_metrics_red_blue_grouped.png")



plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

key_list = [
    "ring_radius_std",
    "ring_radius_mean",
    "dihedral_angle_entropy",
    "bond_angle_std",
"bond_length_std",
"anisotropy_metric_from_structure_factor_bonds",
"bond_orientation_entropy",
"hyperuniformity_alpha",
"critical_pore_radius",
"vertex_homogeneity_metric",
]

labels = [
    Latex.L"Ring radii $\sigma$",
    "Mean ring radius",
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    "Anisotropy in str. fact.",
    "Isotropy in bonds",
    Latex.L"Hyperuniformity $\alpha / 10$",
    "Critical pore radius",
    "Vertex homogeneity",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

analysis_data_paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\predictions\local_relaxation\\"]

filenames = ["pachy_blue_order_metrics.h5", "ctn_predictions_nr_layers_3_nr_neurons_57_full_pca_5_minimal_loss_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    push!(dicts, order_metrics_dict)
end

# Example split indices
n1 = 5


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1
labels2 = labels[n1+1:end]
key_lists2 = key_list[n1+1:end]
yvals2 = 1:length(labels2)


# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

xlims = (-0.1, 1.1)
xticks = [0, 0.5, 1.0]

# Initialize the three subplots
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = xticks,
    xlabel = "Order metric value",
    grid = true
)

# divide the value of hyperuniformity_alpha by 10
for dict in dicts
    dict["hyperuniformity_alpha"] = dict["hyperuniformity_alpha"] / 10
end

# Plot the scatter series for each subplot
markersizes = [10, 8]
for (i, dict) in enumerate(dicts)
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=markersizes[i], ylims=(0.5, n1+0.5), xlims=xlims)

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=markersizes[i], ylims=(0.5, length(labels2)+0.5), xlims=xlims)
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, layout = (2, 1), size = (600, 600), link = :x,
bottom_margin = 3Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 1Plots.mm,
    right_margin = 1Plots.mm)

Plots.savefig(plot_path * "pachy_order_metrics_blue_predicted_grouped.png")


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\predictions\local_relaxation\ctn\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\loss_prediction_ctn_nr_layers_4_nr_neurons_67_full_pca_10\\"

filename = "ctn_predictions_nr_layers_4_nr_neurons_67_full_pca_10.h5"

data_dict = GU.load_h5_dict(analysis_data_path*filename)


predictions_array = data_dict["predictions_array"]
loss_array = data_dict["loss_array"]

# permute dims of the arrays to have the shape (bond_bending_const, t_max, t_gradient )
predictions_array = permutedims(predictions_array, (4, 3, 2, 1))
loss_array = permutedims(loss_array, (3, 2, 1))

bond_bending_const_vec = data_dict["bond_bending_const_vec"]
t_max_vec = data_dict["t_max_vec"]
t_gradient_vec = data_dict["t_gradient_vec"]

# find the 3d window of size (3,3,3) in the loss array where the average loss
# is the smallest
function get_window_smallest_average(loss_array; window_size=3)
    min_loss_value = Inf
    min_i = 0
    min_j = 0
    min_k = 0

    half_window = div(window_size, 2)

    for i in 1+half_window:(size(loss_array, 1)-half_window)
        for j in 1+half_window:(size(loss_array, 2)-half_window)
            for k in 1+half_window:(size(loss_array, 3)-half_window)
                window = loss_array[(i-half_window):(i+half_window), (j-half_window):(j+half_window), (k-half_window):(k+half_window)]
                window_mean = Statistics.mean(window)
                if window_mean < min_loss_value
                    min_loss_value = window_mean
                    min_i = i
                    min_j = j
                    min_k = k
                end
            end
        end
    end
    return (min_i, min_j, min_k), min_loss_value
end

min_loss_index, min_loss_value = get_window_smallest_average(loss_array, window_size=7)

# print the values of the parameters at this index
println("Minimum positive loss: $min_loss_value")
println("At bond_bending_const = $(bond_bending_const_vec[min_loss_index[1]])")
println("At t_max = $(t_max_vec[min_loss_index[2]])")
println("At t_gradient = $(t_gradient_vec[min_loss_index[3]])")

max_loss_value = maximum(loss_array[loss_array .< 999.0])
max_loss_index = findfirst(x -> x == max_loss_value, loss_array)

# print the values of the parameters at this index
println("Maximum loss: $max_loss_value")
println("At bond_bending_const = $(bond_bending_const_vec[max_loss_index[1]])")
println("At t_max = $(t_max_vec[max_loss_index[2]])")
println("At t_gradient = $(t_gradient_vec[max_loss_index[3]])")

# plot a heatmap of the predictions for fixed bond_bending_const = 5.0 with a
# logarithmic color scale
fixed_bond_bending_const = bond_bending_const_vec[min_loss_index[1]]
bond_bending_const_index = min_loss_index[1]

# First fixed part
part1 = [
    "bond_length_std",
    "bond_angle_std",
    "dihedral_angle_entropy",
    "bond_orientation_entropy",
    "coordination_nr_mean",
    "coordination_nr_std"
]

# q_l_value_0 ... q_l_value_12
part2 = ["q_l_value_$(i)" for i in 0:12]

# q_l_uncertainty_0 ... q_l_uncertainty_12
part3 = ["q_l_uncertainty_$(i)" for i in 0:12]

# Last fixed part
part4 = [
    "vertex_homogeneity_metric",
    "ring_size_mean",
    "ring_size_std",
    "ring_radius_mean",
    "ring_radius_std",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor",
    "anisotropy_metric_from_structure_factor_bonds",
    "hyperuniformity_alpha_value",
    "hyperuniformity_alpha_uncertainty"
]

# Combine all parts
order_metrics = vcat(part1, part2, part3, part4)

for (index, metric) in enumerate(order_metrics)
    println("Order metric: $metric ", predictions_array[min_loss_index[1], min_loss_index[2], min_loss_index[3], index])
end


Z  = loss_array[bond_bending_const_index, :, :]'


Zlog = log10.(Z)

exps       = -4:-1
ticks_vals = collect(exps)
ticks_lbls = ["1e$(p)" for p in exps]

# --- key trick for GR: NBSP + newline to create horizontal offset from ticks
# use more "\n" if you need a larger offset
cb_title = "\u00A0\n" * Latex.L"\log(L)"  # NBSP (U+00A0) + newline + LaTeX title

p = Plots.heatmap(
    t_max_vec, t_gradient_vec, Zlog;
    xlabel = Latex.L"T_\mathrm{max}",
    ylabel = Latex.L"\Delta T",
    title  = Latex.L"Loss at $\beta = $" * string(fixed_bond_bending_const),

    # colorbar controls
    colorbar_title       = cb_title,
    colorbar_titlefontsize = 16,            # ↑ label font size
    colorbar_tickfontsize  = 9,             # (optional) shrink tick labels to reduce crowding

    c = :viridis,
    clim = (minimum(exps), maximum(exps)),
    colorbar_ticks = (ticks_vals, ticks_lbls),
    aspect_ratio = :equal,

    # this only affects outer spacing—keep if your figure needs extra room
    right_margin = 1Plots.mm
)



# set xlims and ylims
Plots.heatmap!(; xlims=(minimum(t_max_vec), maximum(t_max_vec)),
                      ylims=(minimum(t_gradient_vec), maximum(t_gradient_vec)))
Plots.savefig(plot_path*"ctn_loss_heatmap_fixed_bond_bending_const_$(fixed_bond_bending_const).png")


# do the same plot again but now for t_gradient fixed at the value of minimal
# loss
fixed_t_gradient = t_gradient_vec[min_loss_index[3]]
t_gradient_index = min_loss_index[3]

Z  = loss_array[:, :, t_gradient_index]'
Zlog = log10.(Z)

Plots.heatmap( t_max_vec, bond_bending_const_vec, Zlog';
    ylabel=Latex.L"\beta",
    xlabel=Latex.L"T_\mathrm{max}",
    title=Latex.L"Loss at $\Delta T = $"*string(fixed_t_gradient),
    # colorbar controls
    colorbar_title       = cb_title,
    colorbar_titlefontsize = 16,            # ↑ label font size
    colorbar_tickfontsize  = 9,             # (optional) shrink tick labels to reduce crowding

    c = :viridis,
    clim = (minimum(exps), maximum(exps)),
    colorbar_ticks = (ticks_vals, ticks_lbls),
    #aspect_ratio = :equal,

    # this only affects outer spacing—keep if your figure needs extra room
    right_margin = 10Plots.mm
)

# set xlims and ylims
Plots.heatmap!(; ylims=(5, maximum(bond_bending_const_vec)),
                      xlims=(minimum(t_max_vec), maximum(t_max_vec)))
Plots.savefig(plot_path*"ctn_loss_heatmap_fixed_t_gradient_$(fixed_t_gradient).png")



load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

filename_red = "pachy_red_structure_factor_bonds_array.h5"
filename_blue = "pachy_blue_structure_factor_bonds_array.h5"

structure_factor_blue_dict = GU.load_h5_dict(load_path*filename_red)
structure_factor_red_dict = GU.load_h5_dict(load_path*filename_blue)

t_range = (1e-2, 5)
min_wavenumber_to_consider = pi/4
consider_spectral_density = false

time_vec = exp.(LinRange(log(t_range[1]), log(t_range[2]), 40))

function comparison_func(x, a) 
    return a * x ^(-3/2)
end

# calculate the excess spreadability for each t value
excess_spreadability_vec_red = [NA.get_excess_spreadability(
        structure_factor_red_dict, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

excess_spreadability_vec_blue = [NA.get_excess_spreadability(
        structure_factor_blue_dict, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

Plots.plot(
    time_vec,
    excess_spreadability_vec_blue;
    xscale = :log10,
    yscale = :log10,
    xlabel = Latex.L"Time $t$",
    ylabel = "Excess spreadability",
    label = "Pachy blue",
    bottom_margin = 2Plots.mm
)
Plots.plot!(
    time_vec,
    excess_spreadability_vec_red;
    label = "Pachy red"
)

Plots.plot!(time_vec, comparison_func.(time_vec, 1000), label = Latex.L"t^{-3/2}", color = :black, linestyle = :dash)


Plots.savefig(plot_path*"pachy_excess_spreadability_comparison.png")



filename_red = "pachy_red_structure_factor_angle_averaged.h5"
filename_blue = "pachy_blue_structure_factor_angle_averaged.h5"

structure_factor_blue_dict = GU.load_h5_dict(load_path*filename_red)
structure_factor_red_dict = GU.load_h5_dict(load_path*filename_blue)

t_range = (1e-2, 5)
min_wavenumber_to_consider = pi/2
consider_spectral_density = false

time_vec = exp.(LinRange(log(t_range[1]), log(t_range[2]), 40))

function comparison_func(x, a) 
    return a * x ^(-3/2)
end

# calculate the excess spreadability for each t value
excess_spreadability_vec_red = [NA.get_excess_spreadability_old(
        structure_factor_red_dict, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

excess_spreadability_vec_blue = [NA.get_excess_spreadability_old(
        structure_factor_blue_dict, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

Plots.plot(
    time_vec,
    Measurements.value.(excess_spreadability_vec_blue);
    xscale = :log10,
    yscale = :log10,
    xlabel = Latex.L"Time $t$",
    ylabel = "Excess spreadability",
    label = "Pachy blue",
    bottom_margin = 2Plots.mm
)
Plots.plot!(
    time_vec,
    Measurements.value.(excess_spreadability_vec_red);
    label = "Pachy red"
)

Plots.plot!(time_vec, comparison_func.(time_vec, 1), label = Latex.L"t^{-3/2}", color = :black, linestyle = :dash)


Plots.savefig(plot_path*"pachy_excess_spreadability_comparison_old_pi_over_2_excluded.png")


load_path_crystal = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\ctn\\"
load_path_disorder = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\targeted\ctn_pachy\target_6\run_2\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\crystals\ctn\\"

filename_crystal = "ctn_1728_vertices_structure_factor_bonds_array.h5"
filename_disorder = "ctn_beta_9.4000_t_max_15.0000_t_gradient_14.0000_nr_vertices_1792_structure_factor_bonds_array.h5"

structure_factor_disorder_dict = GU.load_h5_dict(load_path_disorder*filename_disorder)
structure_factor_crystal_dict = GU.load_h5_dict(load_path_crystal*filename_crystal)

t_range = (1e-2, 5)
min_wavenumber_to_consider = pi/4
consider_spectral_density = false

time_vec = exp.(LinRange(log(t_range[1]), log(t_range[2]), 40))

function comparison_func(x, a) 
    return a * x ^(-3/2)
end

# calculate the excess spreadability for each t value
excess_spreadability_vec_crystal = [NA.get_excess_spreadability(
        structure_factor_crystal_dict, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

excess_spreadability_vec_disorder = [NA.get_excess_spreadability(
        structure_factor_disorder_dict, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

Plots.plot(
    time_vec,
    excess_spreadability_vec_disorder;
    xscale = :log10,
    yscale = :log10,
    xlabel = Latex.L"Time $t$",
    ylabel = "Excess spreadability",
    label = "Disorder",
    bottom_margin = 2Plots.mm
)
Plots.plot!(
    time_vec,
    excess_spreadability_vec_crystal;
    label = "Crystal"
)

Plots.plot!(time_vec, comparison_func.(time_vec, 100), label = Latex.L"t^{-3/2}", color = :black, linestyle = :dash)

# place legend in the bottom left
Plots.plot!(legend = :bottomleft)
Plots.savefig(plot_path*"ctn_1728_vertices_crystal_disorder_excess_spreadability_comparison.png")


load_path_crystal = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\ctn\\"
load_path_disorder = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\targeted\ctn_pachy\target_6\run_2\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\crystals\ctn\\"


filename_crystal = "ctn_1728_vertices_structure_factor_bonds_array.h5"
filename_disorder = "ctn_beta_9.4000_t_max_15.0000_t_gradient_14.0000_nr_vertices_1792_structure_factor_bonds_array.h5"

structure_factor_disorder_dict = GU.load_h5_dict(load_path_disorder*filename_disorder)
structure_factor_crystal_dict = GU.load_h5_dict(load_path_crystal*filename_crystal)

NA.plot_structure_factor_heatmap(
        structure_factor_disorder_dict,
        plot_path*"ctn_disorder_1728_vertices";
        title="disorder",
        save_plot = true,
        clims = (0, 2 ),
        x_y_lims = nothing,
        clims_from_mean = false,)

NA.plot_structure_factor_heatmap(
        structure_factor_crystal_dict,
        plot_path*"ctn_crystal_1728_vertices";
        title="crystal",
        save_plot = true,
        clims = (0, 0.1 ),
        x_y_lims = nothing,
        clims_from_mean = false,)


# place legend in the bottom left
Plots.plot!(legend = :bottomleft)
Plots.savefig(plot_path*"ctn_1728_vertices_crystal_disorder_excess_spreadability_comparison.png")


filename_crystal = "ctn_beta_0.0005_t_max_0.0005_t_gradient_0.0005_structure_factor_bonds_array.h5"
filename_disorder = "ctn_beta_9.4000_t_max_15.0000_t_gradient_14.0000_nr_vertices_224_structure_factor_bonds_array.h5"

structure_factor_disorder_dict = GU.load_h5_dict(load_path_disorder*filename_disorder)
structure_factor_crystal_dict = GU.load_h5_dict(load_path_crystal*filename_crystal)

t_range = (1e-2, 5)
min_wavenumber_to_consider = 0.0 # pi/4
consider_spectral_density = false

time_vec = exp.(LinRange(log(t_range[1]), log(t_range[2]), 40))

function comparison_func(x, a) 
    return a * x ^(-3/2)
end

# calculate the excess spreadability for each t value
excess_spreadability_vec_crystal = [NA.get_excess_spreadability(
        structure_factor_crystal_dict, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

excess_spreadability_vec_disorder = [NA.get_excess_spreadability(
        structure_factor_disorder_dict, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

Plots.plot(
    time_vec,
    excess_spreadability_vec_disorder;
    xscale = :log10,
    yscale = :log10,
    xlabel = Latex.L"Time $t$",
    ylabel = "Excess spreadability",
    label = "Disorder",
    ylim = (1e-4, 1e5),
    bottom_margin = 2Plots.mm
)
Plots.plot!(
    time_vec,
    excess_spreadability_vec_crystal;
    label = "Crystal"
)

Plots.plot!(time_vec, comparison_func.(time_vec, 1), label = Latex.L"t^{-3/2}", color = :black, linestyle = :dash)

# place legend in the bottom left
Plots.plot!(legend = :bottomleft)
Plots.savefig(plot_path*"ctn_224_vertices_crystal_disorder_excess_spreadability_comparison.png")


load_path_crystal = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\ctn\\"
load_path_disorder = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\targeted\ctn_pachy\target_6\run_2\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\crystals\ctn\\"


filename_crystal = "ctn_beta_0.0005_t_max_0.0005_t_gradient_0.0005_structure_factor_bonds_array.h5"
filename_disorder = "ctn_beta_9.4000_t_max_15.0000_t_gradient_14.0000_nr_vertices_224_structure_factor_bonds_array.h5"

structure_factor_disorder_dict = GU.load_h5_dict(load_path_disorder*filename_disorder)
structure_factor_crystal_dict = GU.load_h5_dict(load_path_crystal*filename_crystal)

NA.plot_structure_factor_heatmap(
        structure_factor_disorder_dict,
        plot_path*"ctn_disorder_224_vertices";
        title="disorder",
        save_plot = true,
        clims = (0, 2 ),
        x_y_lims = nothing,
        clims_from_mean = false,)

NA.plot_structure_factor_heatmap(
        structure_factor_crystal_dict,
        plot_path*"ctn_crystal_224_vertices";
        title="crystal",
        save_plot = true,
        clims = (0, 0.1 ),
        x_y_lims = nothing,
        clims_from_mean = false,)


spatial_network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\stern_vir\\"

filename = "stern_vir_blue_structure_factor_bonds_array.h5"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\stern_vir\\"

dict = GU.load_h5_dict(spatial_network_path*filename)

NA.plot_structure_factor_heatmap(
        dict,
        plot_path*"stern_vir_blue_y_z";
        save_plot = true,
        clims = (0, 10 ),
        x_y_lims = nothing,
        wavevector_component_to_fix = 1,
        clims_from_mean = false,)

filename = "stern_vir_green_structure_factor_bonds_array.h5"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\stern_vir\\"

dict = GU.load_h5_dict(spatial_network_path*filename)

NA.plot_structure_factor_heatmap(
        dict,
        plot_path*"stern_vir_green_y_z";
        save_plot = true,
        clims = (0, 10 ),
        x_y_lims = nothing,
        wavevector_component_to_fix = 1,
        clims_from_mean = false,)

        

load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

filename_blue = "pachy_blue_structure_factor_bonds_angle_averaged.h5"

load_path_generated = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\run_5\\"
#raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\Documents\\"
#

filename_generated = "ctn_beta_6.3359_t_max_8.6070_t_gradient_11.7170_structure_factor_bonds_angle_averaged.h5"

structure_factor_blue = GU.load_h5_dict(load_path*filename_blue)

structure_factor_generated = GU.load_h5_dict(load_path_generated*filename_generated)

first_index = 4

Plots.plot(structure_factor_generated["wavenumber_vec"], Measurements.value.(structure_factor_generated["structure_factor_vec"]),
    #ribbon =  Measurements.uncertainty.(structure_factor_generated["structure_factor_vec"]),
    label = "Generated",
    color = :gray
)
Plots.plot!(structure_factor_blue["wavenumber_vec"][first_index:end], Measurements.value.(structure_factor_blue["structure_factor_vec"][first_index:end]),
#ribbon =  Measurements.uncertainty.(structure_factor_blue["structure_factor_vec"]),
    label = "Biological",
    color = :blue
)

Plots.plot!(ylims=(0, 2), xlims=(0, 15))


# Custom ticks: multiples of π
xticks_positions = [π, 2π, 3π, 4π]  # adjust range as needed
xticks_labels = [Latex.L"\pi", Latex.L"2\pi", Latex.L"3\pi", Latex.L"4\pi"]

Plots.xticks!(xticks_positions, xticks_labels)

Plots.xlabel!("Wavenumber / (1/d)")
Plots.ylabel!("Structure factor")
Plots.savefig(plot_path*"pachy_generated_angle_averaged_structure_factor_comparison.png")



load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

filename_blue = "pachy_blue_structure_factor_bonds_array.h5"

load_path_generated = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\run_5\\"

filename_generated = "ctn_beta_6.3359_t_max_8.6070_t_gradient_11.7170_structure_factor_bonds_array.h5"

structure_factor_blue = GU.load_h5_dict(load_path*filename_blue)

structure_factor_generated = GU.load_h5_dict(load_path_generated*filename_generated)

t_range = (5e-2, 5)
min_wavenumber_to_consider = 0.4578267741033627 # pi/4
consider_spectral_density = false

time_vec = exp.(LinRange(log(t_range[1]), log(t_range[2]), 20))

function comparison_func(x, a) 
    return a * x ^(-3/2)
end

# calculate the excess spreadability for each t value
excess_spreadability_vec_blue = [NA.get_excess_spreadability(
        structure_factor_blue, time_vec[i]; 
        min_wavenumber_to_consider=min_wavenumber_to_consider,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

excess_spreadability_vec_generated = [NA.get_excess_spreadability(
        structure_factor_generated, time_vec[i]; 
        min_wavenumber_to_consider=0.0,
        consider_spectral_density = consider_spectral_density) for i in 
        eachindex(time_vec)]

Plots.plot(
    time_vec,
    excess_spreadability_vec_generated .* 3;
    xscale = :log10,
    yscale = :log10,
    xlabel = Latex.L"Time $t$",
    ylabel = "Excess spreadability",
    label = "Generated",
    ylim = (1, 1e4),
    bottom_margin = 2Plots.mm,
    color = :gray
)
Plots.plot!(
    time_vec,
    excess_spreadability_vec_blue .* 1e-3;
    label = "Biological",
    color = :blue
)

Plots.plot!(time_vec, comparison_func.(time_vec, 10), label = Latex.L"t^{-3/2}", color = :black, linestyle = :dash)

# place legend in the bottom left
#Plots.plot!(legend = :bottomleft)
Plots.savefig(plot_path*"pachy_generated_excess_spreadability_comparison.png")



"""
    scatter_order_metrics(
        order_metrics_bio_dict::Dict{String,<:AbstractVector},
        order_metrics_generated_dict::Dict{String,<:AbstractVector};
        metric_x::String,
        metric_y::String,
        bio_index::Int = 1,
        generated_label::AbstractString = "Generated",
        bio_label::AbstractString = "Biological (target)",
        generated_color = :steelblue,
        bio_color = :red,
        generated_marker = :circle,
        bio_marker = :star5,
        generated_ms = 5,
        bio_ms = 9,
        alpha_gen = 0.8,
        legend = :topright,
        title::AbstractString = "",
        savepath::Union{Nothing,String} = nothing
    ) -> Plots.Plot

Create a scatter plot for two order metrics:
- All generated networks (from order_metrics_generated_dict) are plotted as points.
- The target bio network (bio_index in order_metrics_bio_dict) is overlaid as a highlighted point.

Assumes keys exist and vectors are aligned in length.
"""
function scatter_order_metrics(
    order_metrics_bio_dict,
    order_metrics_generated_dict;
    metric_x::String,
    metric_y::String,
    x_label::AbstractString = "",
    y_label::AbstractString = "",
    xlim::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    ylim::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    bio_index::Int = 1,
    highlight_generated_index::Union{Nothing,Int} = nothing,
    generated_label::AbstractString = "Generated",
    bio_label::AbstractString = "Biological",
    generated_color = :gray,
    bio_color = :blue,
    generated_marker = :circle,
    bio_marker = :star5,
    generated_ms = 5,
    bio_ms = 12,
    alpha_gen = 0.8,
    legend = :topright,
    title::AbstractString = "",
    savepath::Union{Nothing,String} = nothing
)
    # Extract data (we rely on your guarantee that keys exist and sizes align)
    x_gen = order_metrics_generated_dict[metric_x]
    y_gen = order_metrics_generated_dict[metric_y]

    x_bio = order_metrics_bio_dict[metric_x][bio_index]
    y_bio = order_metrics_bio_dict[metric_y][bio_index]

    # Base scatter for generated networks
    plt = Plots.scatter(
        x_gen, y_gen;
        label = generated_label,
        color = generated_color,
        marker = (generated_marker, generated_ms),
        alpha = alpha_gen,
        legend = legend,
        title = title,
        xlabel = x_label == "" ? metric_x : x_label,
        ylabel = y_label == "" ? metric_y : y_label,
        xlim = xlim,
        ylim = ylim,
        grid = :on
    )

    # Overlay the target bio point
    Plots.scatter!(
        plt,
        [x_bio], [y_bio];
        label = bio_label,
        color = bio_color,
        marker = (bio_marker, bio_ms),
        alpha = 1.0
    )

    if highlight_generated_index !== nothing
        x_highlight = order_metrics_generated_dict[metric_x][highlight_generated_index]
        y_highlight = order_metrics_generated_dict[metric_y][highlight_generated_index]

        # Overlay the highlighted generated point
        Plots.scatter!(
            plt,
            [x_highlight], [y_highlight];
            label = "Best match",
            color = :black,
            marker = (generated_marker, generated_ms),
            alpha = 1.0
        )
    end

    # Optionally annotate the bio point (uncomment if desired)
    # Plots.annotate!(plt, x_bio, y_bio, text("bio i=$(bio_index)", :left, 10, bio_color))

    # Save if a path is provided
    if savepath !== nothing
        Plots.savefig(plt, savepath)
    end

    return plt
end


bio_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"

generated_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

# "pachy/pachy_blue"
# "pachy/pachy_red"
# "stern_ama/stern_ama_orange"
# "stern_vir/stern_vir_blue"
# "stern_vir/stern_vir_green"

filename = "all_order_metrics.h5"

order_metrics_bio_dict = GU.load_h5_dict(bio_path*filename)

order_metrics_generated_dict = GU.load_h5_dict(generated_path*filename)


metric_x = "hyperuniformity_alpha_vec"
metric_y = "bond_length_std_vec"

plt = scatter_order_metrics(
    order_metrics_bio_dict,
    order_metrics_generated_dict;
    metric_x = metric_x,
    metric_y = metric_y,
    x_label = "Hyperuniformity α",
    y_label = "Bond length "*Latex.L"\sigma",
    xlim = (-2.5, 4.0),
    bio_index = 1,
    highlight_generated_index = 2859,
    savepath = plot_path * "pachy_generated_hyperuniformity_alpha_vs_bond_length_std.png"
)

# save the plot
Plots.savefig(plt, plot_path * "pachy_generated_hyperuniformity_alpha_vs_bond_length_std.png")


bcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\bcu_cn_5_6_7_8\\"

ctn_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\\"

dia_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\dia\\"

lcs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\lcs\\"

pcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\pcu_cn_4_5_6\\"

srs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\srs\\"

bio_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

all_paths = [
    bcu_order_metrics_path,
    ctn_order_metrics_path,
    dia_order_metrics_path,
    lcs_order_metrics_path,
    pcu_order_metrics_path,
    srs_order_metrics_path,
    bio_order_metrics_path
]

all_dicts = [GU.load_h5_dict(path * "all_order_metrics.h5") for path in all_paths]


analysis_data_paths = [bio_order_metrics_path, ctn_order_metrics_path]

filenames = ["pachy\\pachy_blue_order_metrics.h5", "run_5\\ctn_beta_6.3359_t_max_8.6070_t_gradient_11.7170_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    push!(dicts, order_metrics_dict)
end

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "bond_orientation_entropy",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

all_order_metrics_dict = Dict{String, Vector{Float64}}()
for (i, k) in enumerate(key_list)
    current_key = k*"_vec"

    for (j, dict) in enumerate(all_dicts)
        if j == 1
            all_order_metrics_dict[current_key] = dict[current_key]
        else
            append!(all_order_metrics_dict[current_key], dict[current_key])
        end
    end
end

labels = [
    
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    Latex.L"Hyperuniformity $\alpha$",
    "Vertex homogeneity",
    "Critical pore radius",
    "Isotropy in bonds",
    "Mean coordination nr.",
    Latex.L"Ring radii $\sigma$",
    "Mean ring radius",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:end]
key_lists3 = key_list[n2+1:end]
yvals3 = 1:length(labels3)

# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

xlims = (-0.1, 1.1)
xticks = [0, 0.5, 1.0]

# Initialize the three subplots
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p3 = Plots.plot(
    legend = false,
    yticks = (yvals3, labels3),
    xticks = xticks,
    xlabel = "Order metric value",
    grid = true
)


# Plot the scatter series for each subplot
markersizes = [10, 8]
colors = [:blue, :gray]
for (i, dict) in enumerate(dicts)
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=markersizes[i], ylims=(0.5, n1+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=markersizes[i], ylims=(0.5, length(labels2)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p3, [Measurements.value.(dict[k]) for k in key_lists3], yvals3; label=labels3, markersize=markersizes[i], ylims=(0.5, length(labels3)+0.5), xlims=xlims, color=colors[i])
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, layout = (3, 1), size = (500, 600), link = :x,
bottom_margin = 0Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 10Plots.mm,
    right_margin = 0Plots.mm)

Plots.savefig(plot_path * "pachy_generated_order_metrics_comparison.png")


bcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\bcu_cn_5_6_7_8\\"

ctn_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\\"

dia_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\dia\\"

lcs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\lcs\\"

pcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\pcu_cn_4_5_6\\"

srs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\srs\\"

bio_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

all_paths = [
    bcu_order_metrics_path,
    ctn_order_metrics_path,
    dia_order_metrics_path,
    lcs_order_metrics_path,
    pcu_order_metrics_path,
    srs_order_metrics_path,
    bio_order_metrics_path
]

all_dicts = [GU.load_h5_dict(path * "all_order_metrics.h5") for path in all_paths]


analysis_data_paths = [bio_order_metrics_path, ctn_order_metrics_path]

filenames = ["pachy\\pachy_blue_order_metrics.h5", "run_5\\ctn_beta_6.3359_t_max_8.6070_t_gradient_11.7170_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    push!(dicts, order_metrics_dict)
end

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "bond_orientation_entropy",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

all_order_metrics_dict = Dict{String, Vector{Float64}}()
for (i, k) in enumerate(key_list)
    current_key = k*"_vec"

    for (j, dict) in enumerate(all_dicts)
        if j == 1
            all_order_metrics_dict[current_key] = dict[current_key]
        else
            append!(all_order_metrics_dict[current_key], dict[current_key])
        end
    end
end

labels = [
    
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    Latex.L"Hyperuniformity $\alpha$",
    "Vertex homogeneity",
    "Critical pore radius",
    "Isotropy in bonds",
    "Mean coordination nr.",
    Latex.L"Ring radii $\sigma$",
    "Mean ring radius",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:end]
key_lists3 = key_list[n2+1:end]
yvals3 = 1:length(labels3)

# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

# --- Build (min, max) per metric from all_order_metrics_dict (uses *_vec keys) ---
# key_list contains base names (e.g., "bond_length_std"), while all_order_metrics_dict uses "<key>_vec"
mins_maxs = Dict{String,Tuple{Float64,Float64}}()
for k in key_list
    k_vec = k * "_vec"
    vals_all = all_order_metrics_dict[k_vec]                 # Vector{Float64} from all networks
    mn = minimum(vals_all)
    mx = maximum(vals_all)
    mins_maxs[k] = (mn, mx)
end

# Helper: return a vector of relative values for a given dict and a list of base metric keys
rel_values = function (dict::Dict{String,Any}, keys::Vector{String})
    [begin
        x = Measurements.value.(dict[k])                     # scalar (possibly Measurement), take the value
        mn, mx = mins_maxs[k]
        range = mx - mn
        range > 0 ? (x - mn)/range : 0.0                     # guard against zero range
     end for k in keys]
end

# --- Axis cosmetics for relative scale ---
xlims = (-0.05, 1.05)
xticks = [0.0, 0.5, 1.0]

# Initialize the three subplots (unchanged except xticks labels on the last one)
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = (xticks, ["", "", ""]),
    grid = true
)
p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = (xticks, ["", "", ""]),
    grid = true
)
p3 = Plots.plot(
    legend = false,
    yticks = (yvals3, labels3),
    xticks = xticks,
    xlabel = "Relative metric",
    grid = true
)

# --- Plot the scatter series (now using relative values) ---
markersizes = [10, 8]
colors = [:blue, :gray]

for (i, dict) in enumerate(dicts)
    Plots.scatter!(
        p1,
        rel_values(dict, key_lists1), yvals1;
        label = labels1,
        markersize = markersizes[i],
        ylims = (0.5, length(labels1) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
    Plots.scatter!(
        p2,
        rel_values(dict, key_lists2), yvals2;
        label = labels2,
        markersize = markersizes[i],
        ylims = (0.5, length(labels2) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
    Plots.scatter!(
        p3,
        rel_values(dict, key_lists3), yvals3;
        label = labels3,
        markersize = markersizes[i],
        ylims = (0.5, length(labels3) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
end

# Combine and save (unchanged)
Plots.plot(p1, p2, p3,
           layout = (3, 1),
           size = (500, 600),
           link = :x,
           bottom_margin = 0Plots.mm,
           top_margin = 3Plots.mm,
           left_margin = 10Plots.mm,
           right_margin = 0Plots.mm)
Plots.savefig(plot_path * "pachy_generated_order_metrics_comparison_rel.png")


bcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\bcu_cn_5_6_7_8\\"

ctn_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\\"

dia_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\dia\\"

lcs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\lcs\\"

pcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\pcu_cn_4_5_6\\"

srs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\srs\\"

bio_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\stern_vir\\"

all_paths = [
    bcu_order_metrics_path,
    ctn_order_metrics_path,
    dia_order_metrics_path,
    lcs_order_metrics_path,
    pcu_order_metrics_path,
    srs_order_metrics_path,
    bio_order_metrics_path
]

all_dicts = [GU.load_h5_dict(path * "all_order_metrics.h5") for path in all_paths]


analysis_data_paths = [bio_order_metrics_path, bio_order_metrics_path, bcu_order_metrics_path]

filenames = ["stern_vir\\stern_vir_green_order_metrics.h5", "stern_vir\\stern_vir_blue_order_metrics.h5", "run_1\\bcu_cn_5_6_7_8_beta_0.0010_t_max_0.1983_t_gradient_0.1559_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    push!(dicts, order_metrics_dict)
end

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "bond_orientation_entropy",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

all_order_metrics_dict = Dict{String, Vector{Float64}}()
for (i, k) in enumerate(key_list)
    current_key = k*"_vec"

    for (j, dict) in enumerate(all_dicts)
        if j == 1
            all_order_metrics_dict[current_key] = dict[current_key]
        else
            append!(all_order_metrics_dict[current_key], dict[current_key])
        end
    end
end

labels = [
    
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    Latex.L"Hyperuniformity $\alpha$",
    "Vertex homogeneity",
    "Critical pore radius",
    "Isotropy in bonds",
    "Mean coordination nr.",
    Latex.L"Ring radii $\sigma$",
    "Mean ring radius",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:end]
key_lists3 = key_list[n2+1:end]
yvals3 = 1:length(labels3)

# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

# --- Build (min, max) per metric from all_order_metrics_dict (uses *_vec keys) ---
# key_list contains base names (e.g., "bond_length_std"), while all_order_metrics_dict uses "<key>_vec"
mins_maxs = Dict{String,Tuple{Float64,Float64}}()
for k in key_list
    k_vec = k * "_vec"
    vals_all = all_order_metrics_dict[k_vec]                 # Vector{Float64} from all networks
    mn = minimum(vals_all)
    mx = maximum(vals_all)
    mins_maxs[k] = (mn, mx)
end

# Helper: return a vector of relative values for a given dict and a list of base metric keys
rel_values = function (dict::Dict{String,Any}, keys::Vector{String})
    [begin
        x = Measurements.value.(dict[k])                     # scalar (possibly Measurement), take the value
        mn, mx = mins_maxs[k]
        range = mx - mn
        range > 0 ? (x - mn)/range : 0.0                     # guard against zero range
     end for k in keys]
end

# --- Axis cosmetics for relative scale ---
xlims = (-0.05, 1.05)
xticks = [0.0, 0.5, 1.0]

# Initialize the three subplots (unchanged except xticks labels on the last one)
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = (xticks, ["", "", ""]),
    grid = true
)
p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = (xticks, ["", "", ""]),
    grid = true
)
p3 = Plots.plot(
    legend = false,
    yticks = (yvals3, labels3),
    xticks = xticks,
    xlabel = "Relative metric",
    grid = true
)

# --- Plot the scatter series (now using relative values) ---
markersizes = [10, 10, 8]
colors = [:green, :blue, :gray]
alphas = [0.5, 0.5, 1]

for (i, dict) in enumerate(dicts)
    Plots.scatter!(
        p1,
        rel_values(dict, key_lists1), yvals1;
        label = labels1,
        markersize = markersizes[i],
        ylims = (0.5, length(labels1) + 0.5),
        xlims = xlims,
        color = colors[i],
        alpha = alphas[i]
    )
    Plots.scatter!(
        p2,
        rel_values(dict, key_lists2), yvals2;
        label = labels2,
        markersize = markersizes[i],
        ylims = (0.5, length(labels2) + 0.5),
        xlims = xlims,
        color = colors[i],
        alpha = alphas[i]
    )
    Plots.scatter!(
        p3,
        rel_values(dict, key_lists3), yvals3;
        label = labels3,
        markersize = markersizes[i],
        ylims = (0.5, length(labels3) + 0.5),
        xlims = xlims,
        color = colors[i],
        alpha = alphas[i]
    )
end

# Combine and save (unchanged)
Plots.plot(p1, p2, p3,
           layout = (3, 1),
           size = (500, 600),
           link = :x,
           bottom_margin = 0Plots.mm,
           top_margin = 3Plots.mm,
           left_margin = 10Plots.mm,
           right_margin = 0Plots.mm)
Plots.savefig(plot_path * "stern_vir_generated_bcu_order_metrics_comparison_rel.png")


bcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\bcu_cn_5_6_7_8\\"

ctn_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\\"

dia_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\dia\\"

lcs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\lcs\\"

pcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\pcu_cn_4_5_6\\"

srs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\srs\\"

bio_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\stern_ama\\"

all_paths = [
    bcu_order_metrics_path,
    ctn_order_metrics_path,
    dia_order_metrics_path,
    lcs_order_metrics_path,
    pcu_order_metrics_path,
    srs_order_metrics_path,
    bio_order_metrics_path
]

all_dicts = [GU.load_h5_dict(path * "all_order_metrics.h5") for path in all_paths]


analysis_data_paths = [bio_order_metrics_path, pcu_order_metrics_path, ctn_order_metrics_path]

filenames = ["stern_ama\\stern_ama_orange_order_metrics.h5", "run_3\\pcu_cn_4_5_6_beta_0.2195_t_max_0.4752_t_gradient_0.1001_order_metrics.h5", "run_5\\ctn_beta_6.3359_t_max_8.6070_t_gradient_11.7170_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    push!(dicts, order_metrics_dict)
end

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "bond_orientation_entropy",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

all_order_metrics_dict = Dict{String, Vector{Float64}}()
for (i, k) in enumerate(key_list)
    current_key = k*"_vec"

    for (j, dict) in enumerate(all_dicts)
        if j == 1
            all_order_metrics_dict[current_key] = dict[current_key]
        else
            append!(all_order_metrics_dict[current_key], dict[current_key])
        end
    end
end

labels = [
    
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    Latex.L"Hyperuniformity $\alpha$",
    "Vertex homogeneity",
    "Critical pore radius",
    "Isotropy in bonds",
    "Mean coordination nr.",
    Latex.L"Ring radii $\sigma$",
    "Mean ring radius",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:end]
key_lists3 = key_list[n2+1:end]
yvals3 = 1:length(labels3)

# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

# --- Build (min, max) per metric from all_order_metrics_dict (uses *_vec keys) ---
# key_list contains base names (e.g., "bond_length_std"), while all_order_metrics_dict uses "<key>_vec"
mins_maxs = Dict{String,Tuple{Float64,Float64}}()
for k in key_list
    k_vec = k * "_vec"
    vals_all = all_order_metrics_dict[k_vec]                 # Vector{Float64} from all networks
    mn = minimum(vals_all)
    mx = maximum(vals_all)
    mins_maxs[k] = (mn, mx)
end

# Helper: return a vector of relative values for a given dict and a list of base metric keys
rel_values = function (dict::Dict{String,Any}, keys::Vector{String})
    [begin
        x = Measurements.value.(dict[k])                     # scalar (possibly Measurement), take the value
        mn, mx = mins_maxs[k]
        range = mx - mn
        range > 0 ? (x - mn)/range : 0.0                     # guard against zero range
     end for k in keys]
end

# --- Axis cosmetics for relative scale ---
xlims = (-0.05, 1.05)
xticks = [0.0, 0.5, 1.0]

# Initialize the three subplots (unchanged except xticks labels on the last one)
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = (xticks, ["", "", ""]),
    grid = true
)
p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = (xticks, ["", "", ""]),
    grid = true
)
p3 = Plots.plot(
    legend = false,
    yticks = (yvals3, labels3),
    xticks = xticks,
    xlabel = "Relative metric",
    grid = true
)

# --- Plot the scatter series (now using relative values) ---
markersizes = [10, 8, 8]
colors = [:orange, :gray, :black]

for (i, dict) in enumerate(dicts)
    Plots.scatter!(
        p1,
        rel_values(dict, key_lists1), yvals1;
        label = labels1,
        markersize = markersizes[i],
        ylims = (0.5, length(labels1) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
    Plots.scatter!(
        p2,
        rel_values(dict, key_lists2), yvals2;
        label = labels2,
        markersize = markersizes[i],
        ylims = (0.5, length(labels2) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
    Plots.scatter!(
        p3,
        rel_values(dict, key_lists3), yvals3;
        label = labels3,
        markersize = markersizes[i],
        ylims = (0.5, length(labels3) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
end

# Combine and save (unchanged)
Plots.plot(p1, p2, p3,
           layout = (3, 1),
           size = (500, 600),
           link = :x,
           bottom_margin = 0Plots.mm,
           top_margin = 3Plots.mm,
           left_margin = 10Plots.mm,
           right_margin = 0Plots.mm)
Plots.savefig(plot_path * "stern_ama_generated_ctn_order_metrics_comparison_rel.png")


bcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\bcu_cn_5_6_7_8\\"

ctn_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\\"

dia_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\dia\\"

lcs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\lcs\\"

pcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\pcu_cn_4_5_6\\"

srs_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\srs\\"

bio_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\pachy\\"

all_paths = [
    bcu_order_metrics_path,
    ctn_order_metrics_path,
    dia_order_metrics_path,
    lcs_order_metrics_path,
    pcu_order_metrics_path,
    srs_order_metrics_path,
    bio_order_metrics_path
]

all_dicts = [GU.load_h5_dict(path * "all_order_metrics.h5") for path in all_paths]


analysis_data_paths = [bio_order_metrics_path, ctn_order_metrics_path]

filenames = ["pachy\\pachy_blue_order_metrics.h5", "run_5\\ctn_beta_6.3359_t_max_8.6070_t_gradient_11.7170_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    push!(dicts, order_metrics_dict)
end

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "bond_orientation_entropy",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

all_order_metrics_dict = Dict{String, Vector{Float64}}()
for (i, k) in enumerate(key_list)
    current_key = k*"_vec"

    for (j, dict) in enumerate(all_dicts)
        if j == 1
            all_order_metrics_dict[current_key] = dict[current_key]
        else
            append!(all_order_metrics_dict[current_key], dict[current_key])
        end
    end
end

labels = [
    
    "Dihedral angle entropy",
    Latex.L"Bond angles $\sigma$",
    Latex.L"Bond lengths $\sigma$",
    Latex.L"Hyperuniformity $\alpha$",
    "Vertex homogeneity",
    "Critical pore radius",
    "Isotropy in bonds",
    "Mean coordination nr.",
    Latex.L"Ring radii $\sigma$",
    "Mean ring radius",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:end]
key_lists3 = key_list[n2+1:end]
yvals3 = 1:length(labels3)

# Calculate relative heights based on number of labels
total_labels = length(labels)
heights = [length(labels1)/total_labels, length(labels2)/total_labels]


# Create the layout for the three subplots
l = @Plots.layout([a b; c d])

# --- Build (min, max) per metric from all_order_metrics_dict (uses *_vec keys) ---
# key_list contains base names (e.g., "bond_length_std"), while all_order_metrics_dict uses "<key>_vec"
mins_maxs = Dict{String,Tuple{Float64,Float64}}()
for k in key_list
    k_vec = k * "_vec"
    vals_all = all_order_metrics_dict[k_vec]                 # Vector{Float64} from all networks
    mn = minimum(vals_all)
    mx = maximum(vals_all)
    mins_maxs[k] = (mn, mx)
end

# Helper: return a vector of relative values for a given dict and a list of base metric keys
rel_values = function (dict::Dict{String,Any}, keys::Vector{String})
    [begin
        x = Measurements.value.(dict[k])                     # scalar (possibly Measurement), take the value
        mn, mx = mins_maxs[k]
        range = mx - mn
        range > 0 ? (x - mn)/range : 0.0                     # guard against zero range
     end for k in keys]
end

# --- Axis cosmetics for relative scale ---
xlims = (-0.05, 1.05)
xticks = [0.0, 0.5, 1.0]

# Initialize the three subplots (unchanged except xticks labels on the last one)
p1 = Plots.plot(
    legend = false,
    yticks = (yvals1, labels1),
    xticks = (xticks, ["", "", ""]),
    grid = true
)
p2 = Plots.plot(
    legend = false,
    yticks = (yvals2, labels2),
    xticks = (xticks, ["", "", ""]),
    grid = true
)
p3 = Plots.plot(
    legend = false,
    yticks = (yvals3, labels3),
    xticks = xticks,
    xlabel = "Relative metric",
    grid = true
)

# --- Plot the scatter series (now using relative values) ---
markersizes = [10, 8]
colors = [:blue, :gray]

for (i, dict) in enumerate(dicts)
    Plots.scatter!(
        p1,
        rel_values(dict, key_lists1), yvals1;
        label = labels1,
        markersize = markersizes[i],
        ylims = (0.5, length(labels1) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
    Plots.scatter!(
        p2,
        rel_values(dict, key_lists2), yvals2;
        label = labels2,
        markersize = markersizes[i],
        ylims = (0.5, length(labels2) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
    Plots.scatter!(
        p3,
        rel_values(dict, key_lists3), yvals3;
        label = labels3,
        markersize = markersizes[i],
        ylims = (0.5, length(labels3) + 0.5),
        xlims = xlims,
        color = colors[i]
    )
end

# Combine and save (unchanged)
Plots.plot(p1, p2, p3,
           layout = (3, 1),
           size = (500, 600),
           link = :x,
           bottom_margin = 0Plots.mm,
           top_margin = 3Plots.mm,
           left_margin = 10Plots.mm,
           right_margin = 0Plots.mm)
Plots.savefig(plot_path * "pachy_generated_order_metrics_comparison_rel.png")


spatial_network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\networks\stern_ama\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\biological\networks\stern_ama\\"

filename = "stern_ama_orange"

#spatial_network = NG.load_spatial_network_from_gml(
#    spatial_network_path*filename*".gml")


exclude_layer_thickness = 1.5
periodic_boundary_conditions = false

# get minimal and maximal vertex coordinates along all three axes in case
# of non-periodic boundary conditions
if !periodic_boundary_conditions
    min_vertex_coords, max_vertex_coords = NA.get_min_max_vertex_coords(
    spatial_network)
end
# get all considered vertices
considered_vertices = Vector{Int64}()
for vertex in MetaGraphsNext.labels(spatial_network)
    if periodic_boundary_conditions
        push!(considered_vertices, vertex)
    else
        # get the positions of the vertex
        vertex_pos = spatial_network[vertex]["position"]
        if (all(vertex_pos 
                .> (min_vertex_coords .+ exclude_layer_thickness))
            && all(vertex_pos 
                .< (max_vertex_coords .- exclude_layer_thickness)))
            push!(considered_vertices, vertex)
        end
    end
end
coordination_nr_vec = Vector{Int64}(
    undef, length(considered_vertices))
for (i, vertex) in enumerate(considered_vertices)
    coordination_nr_vec[i] = length(
        MetaGraphsNext.neighbor_labels(spatial_network, vertex))
end
# remove all elements of the coordination_nr_vec that are 1
coordination_nr_vec = filter(x -> x != 1, coordination_nr_vec)

# plot histogram of coordination numbers

# plot histogram of coordination numbers as relative frequencies
Plots.histogram(
    coordination_nr_vec;
    bins = minimum(coordination_nr_vec):maximum(coordination_nr_vec),
    xlabel = "Coordination number",
    ylabel = "Relative frequency",
    legend = false,
    xticks = minimum(coordination_nr_vec):maximum(coordination_nr_vec),
    xlim = (minimum(coordination_nr_vec)-0.5,
            maximum(coordination_nr_vec)+0.5),
    size = (800, 600),
    normalize = true,
    )


Plots.savefig(plot_path*filename*"_coord_nr_excluded_layer_"*string(exclude_layer_thickness)*".png")


plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\t_melt_network_types\\"

network_types = ["ctn", "dia", "lcs", "srs", "bcu_cn_5_6_7_8", "pcu_cn_4_5_6"]

beta_vec = collect(0.0:0.25:10.0)

Plots.plot()

for network_type in network_types
    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in beta_vec]

    Plots.plot!(beta_vec, t_melt_vec; label=network_type)
end

Plots.xlabel!("Bond bending constant β")
Plots.ylabel!("Melting temperature")
Plots.savefig(plot_path*"melting_temperatures_absolute.png")

Plots.plot()

for network_type in network_types
    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in beta_vec]
    # normalize by maximal melting temperature
    t_melt_vec = t_melt_vec ./ maximum(t_melt_vec)

    Plots.plot!(beta_vec, t_melt_vec; label=network_type)
end

Plots.xlabel!("Bond bending constant β")
Plots.ylabel!("Melting temperature (normalized)")
Plots.savefig(plot_path*"melting_temperatures_relative.png")



plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\t_melt_network_types\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\\"

filename = "all_order_metrics.h5"

network_types = ["ctn", "dia", "lcs", "srs", "bcu_cn_5_6_7_8", "pcu_cn_4_5_6"]

beta_vec = collect(0.0:0.25:10.0)

for network_type in network_types
    

    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    beta_vec = order_metrics_generated_dict["bond_bending_const_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in beta_vec]

    bond_orientation_entropy_vec = order_metrics_generated_dict["bond_orientation_entropy_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the bond orientation entropy
    
    # first, convert the bond orientation entropy to a color scale between 0 and 1
    bond_orientation_entropy_min = minimum(bond_orientation_entropy_vec)
    bond_orientation_entropy_max = maximum(bond_orientation_entropy_vec)
    bond_orientation_entropy_normalized_vec = (bond_orientation_entropy_vec .- bond_orientation_entropy_min) ./ (bond_orientation_entropy_max - bond_orientation_entropy_min)
    #colors = Plots.cgrad(:roma, 100)[round.(Int, bond_orientation_entropy_normalized_vec .* 99) .+ 1]
    colors = Plots.cgrad(:bluesreds)[bond_orientation_entropy_normalized_vec]

    #filter all vector to values of beta between 5 and 10
    filter_indices = findall(x -> x >= 0.0 && x <= 10.0, beta_vec)
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    t_gradient_normalized_vec_filtered = t_gradient_normalized_vec[filter_indices]
    colors_filtered = colors[filter_indices]
    bond_orientation_entropy_vec_filtered = bond_orientation_entropy_vec[filter_indices]

    xlims=(0.5, 2.0)
    ylims=(0.25, 2.0)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        t_gradient_normalized_vec_filtered,
        zcolor = bond_orientation_entropy_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (bond_orientation_entropy_min, 1.0),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n Bond Orientation Entropy",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 5,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"T_\mathrm{gradient} / T_\mathrm{melt}",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    Plots.savefig(p, plot_path*network_type*"_tmax_tgradient_bond_orientation_entropy_beta_greater_0_less_10.png")

    #filter all vector to values of beta between 5 and 10
    filter_indices = findall(x -> x >= 1.0 && x <= 10.0, beta_vec)
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    t_gradient_normalized_vec_filtered = t_gradient_normalized_vec[filter_indices]
    bond_orientation_entropy_vec_filtered = bond_orientation_entropy_vec[filter_indices]

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        t_gradient_normalized_vec_filtered,
        zcolor = bond_orientation_entropy_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (bond_orientation_entropy_min, 1.0),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n Bond Orientation Entropy",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 5,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"T_\mathrm{gradient} / T_\mathrm{melt}",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    Plots.savefig(p, plot_path*network_type*"_tmax_tgradient_bond_orientation_entropy_beta_greater_1_less_10.png")

    #filter all vector to values of beta between 5 and 10
    filter_indices = findall(x -> x >= 0.0 && x <= 1.0, beta_vec)
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    t_gradient_normalized_vec_filtered = t_gradient_normalized_vec[filter_indices]
    bond_orientation_entropy_vec_filtered = bond_orientation_entropy_vec[filter_indices]

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        t_gradient_normalized_vec_filtered,
        zcolor = bond_orientation_entropy_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (bond_orientation_entropy_min, 1.0),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n Bond Orientation Entropy",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 5,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"T_\mathrm{gradient} / T_\mathrm{melt}",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    Plots.savefig(p, plot_path*network_type*"_tmax_tgradient_bond_orientation_entropy_beta_greater_0_less_1.png")
end


plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\t_max_bond_length_beta_network_types\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\\"

filename = "all_order_metrics.h5"

network_types = ["ctn", "dia", "lcs", "srs", "bcu_cn_5_6_7_8", "pcu_cn_4_5_6"]

t_gradient_vec = collect(0.0:0.25:10.0)

for network_type in network_types
    

    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    bond_length_std_vec = order_metrics_generated_dict["bond_length_std_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the bond length std

    # first, convert the bond length std to a color scale between 0 and 1
    bond_length_std_min = minimum(bond_length_std_vec)
    bond_length_std_max = maximum(bond_length_std_vec)
    bond_length_std_normalized_vec = (bond_length_std_vec .- bond_length_std_min) ./ (bond_length_std_max - bond_length_std_min)
    
    colors = Plots.cgrad(:bluesreds)[bond_length_std_normalized_vec]

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    colors_filtered = colors[filter_indices]
    bond_length_std_vec_filtered = bond_length_std_vec[filter_indices]

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = bond_length_std_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (bond_length_std_min, bond_length_std_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n Bond length std.",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 5,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    Plots.savefig(p, plot_path*network_type*"_tmax_beta_bond_length_std_beta_greater_0.25_less_2.png")

end



plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\t_max_bond_angle_beta_network_types\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\\"

filename = "all_order_metrics.h5"

network_types = ["ctn", "dia", "lcs", "srs", "bcu_cn_5_6_7_8", "pcu_cn_4_5_6"]

t_gradient_vec = collect(0.0:0.25:10.0)

for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    bond_angle_std_vec = order_metrics_generated_dict["bond_angle_std_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the Bond angle std

    # first, convert the Bond angle std to a color scale between 0 and 1
    

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    bond_angle_std_vec_filtered = bond_angle_std_vec[filter_indices]

    bond_angle_std_min = minimum(bond_angle_std_vec_filtered)
    bond_angle_std_max = maximum(bond_angle_std_vec_filtered)
    bond_angle_std_normalized_vec = (bond_angle_std_vec_filtered .- bond_angle_std_min) ./ (bond_angle_std_max - bond_angle_std_min)

    #colors = Plots.cgrad(:bluesreds)[bond_angle_std_normalized_vec]

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = bond_angle_std_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (bond_angle_std_min, bond_angle_std_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n Bond angle std.",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 5,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    Plots.savefig(p, plot_path*network_type*"_tmax_beta_bond_angle_std_beta_greater_0.25_less_2.png")

end



plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\t_max_nr_accepted_moves_beta_network_types\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\\"
network_path = replace(analysis_data_path, "analysis_data" => "structures")

filename = "all_order_metrics_with_nr_accepted_moves.h5"

network_types = ["ctn", "dia", "lcs", "srs", "bcu_cn_5_6_7_8", "pcu_cn_4_5_6"]

t_gradient_vec = collect(0.0:0.25:10.0)

for network_type in network_types
    

    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]
    nr_accepted_moves_vec = order_metrics_generated_dict["nr_accepted_moves_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the Bond angle std

    # first, convert the Bond angle std to a color scale between 0 and 1
    

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    nr_accepted_moves_vec_filtered = nr_accepted_moves_vec[filter_indices]

    nr_accepted_moves_min = minimum(nr_accepted_moves_vec_filtered)
    nr_accepted_moves_max = maximum(nr_accepted_moves_vec_filtered)
    nr_accepted_moves_normalized_vec = (nr_accepted_moves_vec_filtered .- nr_accepted_moves_min) ./ (nr_accepted_moves_max - nr_accepted_moves_min)


    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    z = Float64.(copy(nr_accepted_moves_vec_filtered))
    z[z .< 10.0] .= NaN                      # mark zeros as NaN

    # Positive-only clims for log scale
    zmin = minimum(skipmissing(z[.!isnan.(z)]))  # smallest positive
    zmax = maximum(skipmissing(z[.!isnan.(z)]))

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = z,
        c = :bluesreds,
        clims = (10.0, 1000.0),
        colorbar = true,
        colorbar_scale = :log10,         # logarithmic colorbar for positive values
        nan_color = :black,              # <-- zeros (set to NaN) shown in black
        colorbar_title = "\n Nr accepted moves",
        colorbar_title_location = :right,
        markersize = 5,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
        topmargin = 3Plots.mm,
    )

    Plots.savefig(p, plot_path*network_type*"_tmax_beta_nr_accepted_moves.png")

end


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\structure_factor_comparison\\"

network_type_vec = ["ctn", "dia",  "bcu_cn_5_6_7_8", "lcs", "pcu_cn_4_5_6", "srs",] #["ctn",  "bcu_cn_5_6_7_8"] #"dia",  "lcs", "srs", "pcu_cn_4_5_6",

filename_ordered_vec = ["run_6/ctn_beta_1.7313_t_max_2.0550_t_gradient_2.1466", "run_3/dia_beta_6.3550_t_max_24.2721_t_gradient_8.6274", "run_2/bcu_cn_5_6_7_8_beta_1.9736_t_max_294.9035_t_gradient_702.9766", "run_2/lcs_beta_0.0127_t_max_0.0057_t_gradient_0.0139", "run_3/pcu_cn_4_5_6_beta_0.2258_t_max_0.2197_t_gradient_0.4606", "run_2/srs_beta_0.0394_t_max_0.0051_t_gradient_0.0020"]
filename_disordered_vec = ["run_5/ctn_beta_6.5833_t_max_13.7932_t_gradient_13.5214", "run_1/dia_beta_0.7422_t_max_0.7826_t_gradient_0.1969", "run_1/bcu_cn_5_6_7_8_beta_0.5163_t_max_158.6786_t_gradient_20.1955", "run_3/lcs_beta_1.2760_t_max_1.8437_t_gradient_1.8422", "run_3/pcu_cn_4_5_6_beta_0.3229_t_max_5.2818_t_gradient_2.0898", "run_2/srs_beta_1.2625_t_max_0.6123_t_gradient_0.5012"]
i=1
structure_factor_angle_averaged_dict_1 = GU.load_h5_dict(analysis_data_path*network_type_vec[i]*raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_angle_averaged.h5")
structure_factor_dict_1 = GU.load_h5_dict(analysis_data_path*network_type_vec[i] *raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_array.h5")

i=2
structure_factor_angle_averaged_dict_2 = GU.load_h5_dict(analysis_data_path*network_type_vec[i]*raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_angle_averaged.h5")
structure_factor_dict_2 = GU.load_h5_dict(analysis_data_path*network_type_vec[i] *raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_array.h5")

i=3
structure_factor_angle_averaged_dict_3 = GU.load_h5_dict(analysis_data_path*network_type_vec[i]*raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_angle_averaged.h5")
structure_factor_dict_3 = GU.load_h5_dict(analysis_data_path*network_type_vec[i] *raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_array.h5")

i=4
structure_factor_angle_averaged_dict_4 = GU.load_h5_dict(analysis_data_path*network_type_vec[i]*raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_angle_averaged.h5")
structure_factor_dict_4 = GU.load_h5_dict(analysis_data_path*network_type_vec[i] *raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_array.h5")

i=5
structure_factor_angle_averaged_dict_5 = GU.load_h5_dict(analysis_data_path*network_type_vec[i]*raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_angle_averaged.h5")
structure_factor_dict_5 = GU.load_h5_dict(analysis_data_path*network_type_vec[i] *raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_array.h5")

i=6
structure_factor_angle_averaged_dict_6 = GU.load_h5_dict(analysis_data_path*network_type_vec[i]*raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_angle_averaged.h5")
structure_factor_dict_6 = GU.load_h5_dict(analysis_data_path*network_type_vec[i] *raw"\\"*filename_ordered_vec[i]*"_structure_factor_bonds_array.h5")


wavenumbers_to_check_vec = (2*pi) .* [4/8, 5/8, 6/8, 7/8, 1.0]


Plots.plot()
Plots.plot(structure_factor_angle_averaged_dict_1["unfiltered_wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict_1["unfiltered_structure_factor_vec"]), ribbon=Measurements.uncertainty.(structure_factor_angle_averaged_dict_1["unfiltered_structure_factor_vec"]), label="Ordered $(network_type_vec[1])", xlabel="Wavenumber", ylabel="Structure Factor")
Plots.plot!(structure_factor_angle_averaged_dict_2["unfiltered_wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict_2["unfiltered_structure_factor_vec"]), ribbon=Measurements.uncertainty.(structure_factor_angle_averaged_dict_2["unfiltered_structure_factor_vec"]), label="Ordered $(network_type_vec[2])")
Plots.plot!(structure_factor_angle_averaged_dict_3["unfiltered_wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict_3["unfiltered_structure_factor_vec"]), ribbon=Measurements.uncertainty.(structure_factor_angle_averaged_dict_3["unfiltered_structure_factor_vec"]), label="Ordered $(network_type_vec[3])")
Plots.plot!(structure_factor_angle_averaged_dict_4["unfiltered_wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict_4["unfiltered_structure_factor_vec"]), ribbon=Measurements.uncertainty.(structure_factor_angle_averaged_dict_4["unfiltered_structure_factor_vec"]), label="Ordered $(network_type_vec[4])")
Plots.plot!(structure_factor_angle_averaged_dict_5["unfiltered_wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict_5["unfiltered_structure_factor_vec"]), ribbon=Measurements.uncertainty.(structure_factor_angle_averaged_dict_5["unfiltered_structure_factor_vec"]), label="Ordered $(network_type_vec[5])")
Plots.plot!(structure_factor_angle_averaged_dict_6["unfiltered_wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict_6["unfiltered_structure_factor_vec"]), ribbon=Measurements.uncertainty.(structure_factor_angle_averaged_dict_6["unfiltered_structure_factor_vec"]), label="Ordered $(network_type_vec[6])")

# plot vertical lines at the wavenumbers to check
for k in eachindex(wavenumbers_to_check_vec)
    Plots.vline!([wavenumbers_to_check_vec[k]], linestyle=:dash, color=:black, alpha=0.3, label="")
end

Plots.xlims!(0, 10.0)

Plots.savefig(plot_path*"structure_factors_angle_averaged_comparison.png")