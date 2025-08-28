
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex
import Measurements
import GLMakie

import CSV
import DataFrames
import Statistics

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
bottom_margin = 1Plots.mm,
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

## Add group labels using annotations
#Plots.annotate!(p1, -0.5, Statistics.mean(yvals1), Plots.text("Local", :center, 12, :bold, :rotation => 90))
#Plots.annotate!(p2, -0.5, Statistics.mean(yvals2), Plots.text("Global", :center, 12, :bold, :rotation => 90))
#Plots.annotate!(p3, -0.5, Statistics.mean(yvals3), Plots.text("Topological", :center, 12, :bold, :rotation => 90))
#

Plots.savefig(save_path * "order_relative_peak_height_pearson_correlations_dia_grouped.png")

