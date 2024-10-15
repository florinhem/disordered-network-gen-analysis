
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
    corr_fct_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_correlation_functions.h5")
    Plots.plot!(corr_fct_dict["vertex_distance_vec"], corr_fct_dict["pair_correlation_fct_vec"], label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "distance", ylabel = "pair correlation function")

Plots.savefig(path*"pair_correlation_function_low_t.png")

Plots.plot()

for temperature in temperature_vec
    corr_fct_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_correlation_functions.h5")
    Plots.plot!(corr_fct_dict["vertex_distance_vec"], corr_fct_dict["ripley_k_vec"], label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "distance", ylabel = "Cumulative coord. nr.", ylims=(0,6))

Plots.savefig(path*"cumulative_coord_nr_low_t_wrong.png")

Plots.plot()

for temperature in temperature_vec
    corr_fct_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_correlation_functions.h5")

    supercell_edge_length = 11.547005383792516
    nr_vertices = 1000
    vertex_density = nr_vertices/supercell_edge_length^3
    cumulative_coord_nr_vec = Vector{Float64}(undef, 
    length(corr_fct_dict["vertex_distance_vec"]))
    vertex_distance_step = (corr_fct_dict["vertex_distance_vec"][2] - corr_fct_dict["vertex_distance_vec"][1])
    
    cumulative_coord_nr_vec = (vertex_density*4*pi *vertex_distance_step * cumsum(corr_fct_dict["vertex_distance_vec"] .^2 .* corr_fct_dict["pair_correlation_fct_vec"] ))
    
    
    Plots.plot!(corr_fct_dict["vertex_distance_vec"], cumulative_coord_nr_vec, label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "distance", ylabel = "Cumulative coord. nr.", ylims=(0,6))

Plots.savefig(path*"cumulative_coord_nr_low_t.png")

temperature_vec = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))[7:end]

Plots.plot()

for temperature in temperature_vec
    corr_fct_dict = GU.load_h5_dict(load_path * "1000_vertices_T_" * string(temperature) * "_heat_cool_0.1_per_mc_quenched_correlation_functions.h5")
    Plots.plot!(corr_fct_dict["vertex_distance_vec"], corr_fct_dict["pair_correlation_fct_vec"], label = "T = " * string(temperature))
end

Plots.plot!(xlabel = "distance", ylabel = "pair correlation function")

Plots.savefig(path*"pair_correlation_function_high_t.png")