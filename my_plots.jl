
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
#Plots.pyplot()
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

function load_table_as_dict(filename::String)
    df = CSV.File(filename; delim='\t') |> DataFrames.DataFrame
    return Dict(col => collect(df[!, col]) for col in names(df))
end


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