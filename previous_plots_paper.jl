
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
colorbar_titlefontsize=fontsize,
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


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

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

colors = Plots.palette(:tab10)
color_high_T = colors[2]
color_low_T = colors[1]

Plots.plot()

Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.4, 0.7),  alpha=1.0, color=color_high_T)
Plots.plot!(x_vec, heat_cool_temperature_vec.(x_vec, 0.2, 1.0),  alpha=1.0, color=color_low_T)
Plots.plot!(legend = false, xlabel="Attempts per bond chain", ylabel=Latex.L"kT", right_margin = 4Plots.mm, size = (450, 300), xlims=(0,2))

Plots.savefig(path*"heat_cool_temperature_profile_5.pdf")


function get_coord_nr_histogram(
    spatial_network::MetaGraphsNext.MetaGraph;
    periodic_boundary_conditions::Bool=true,
    exclude_layer_thickness::Float64=0.0,
    max_coord_nr::Int64=5)

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

    histogram = StatsBase.fit(StatsBase.Histogram, coordination_nr_vec, 
        1.5:1:max_coord_nr+0.5)

    histogram = LinearAlgebra.normalize(histogram, mode=:probability)

    return histogram
end

path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

ctn_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\crystals\\"

pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\networks\pachy\\"

ctn_spatial_network = NG.load_spatial_network_from_gml(ctn_path*"ctn_1728_vertices.gml")

pachy_spatial_network = NG.load_spatial_network_from_gml(pachy_path*"pachy_blue.gml")

ctn_histogram = get_coord_nr_histogram(
    ctn_spatial_network;
    periodic_boundary_conditions=true,
    exclude_layer_thickness=0.0)

pachy_histogram = get_coord_nr_histogram(
    pachy_spatial_network;
    periodic_boundary_conditions=false,
    exclude_layer_thickness=1.5)


colors = Plots.palette(:tab10)

edges = ctn_histogram.edges[1][2:end-1] .+ 1
ctn_weights = ctn_histogram.weights[2:end]
pachy_weights = pachy_histogram.weights[2:end]

# Define a small offset for side-by-side bars
offset = 0.15  # adjust based on bar width and spacing
total_offset = 0.5

# Plot Pachy bars shifted left
Plots.bar(
    edges .- offset .- total_offset,
    pachy_weights;
    bar_width = 0.3,        # narrower bars
    color = colors[1],
    label = "Pachy",
    linecolor = :transparent,
    legend=false
)

# Overlay CTN bars shifted right
Plots.bar!(
    edges .+ offset .- total_offset,
    ctn_weights;
    bar_width = 0.3,
    color = :grey,
    label = "CTN",
    linecolor = :transparent
)

Plots.xlabel!("Coordination number")
Plots.ylabel!("Frequency")
Plots.xticks!(1:5)

# set the figsize
Plots.plot!(size = (450, 300))

Plots.savefig(path*"pachy_ctn_coordination_number_histogram.pdf")



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

pcu_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\crystals\\"

stern_ama_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\networks\stern_ama\\"

pcu_spatial_network = NG.load_spatial_network_from_gml(pcu_path*"pcu_cn_4_5_6_216_vertices.gml")

stern_ama_spatial_network = NG.load_spatial_network_from_gml(stern_ama_path*"stern_ama_orange.gml")

pcu_histogram = get_coord_nr_histogram(
    pcu_spatial_network;
    periodic_boundary_conditions=true,
    exclude_layer_thickness=0.0)

stern_ama_histogram = get_coord_nr_histogram(
    stern_ama_spatial_network;
    periodic_boundary_conditions=false,
    exclude_layer_thickness=1.5)


colors = Plots.palette(:tab10)

edges = pcu_histogram.edges[1][2:end-1] .+ 1
pcu_weights = pcu_histogram.weights[2:end]
stern_ama_weights = stern_ama_histogram.weights[2:end]

# Define a small offset for side-by-side bars
offset = 0.15  # adjust based on bar width and spacing
total_offset = 0.5

# Plot stern_ama bars shifted left
Plots.bar(
    edges .- offset .- total_offset,
    stern_ama_weights;
    bar_width = 0.3,        # narrower bars
    color = colors[2],
    label = "stern_ama",
    linecolor = :transparent,
    legend=false
)

# Overlay pcu bars shifted right
Plots.bar!(
    edges .+ offset .- total_offset,
    pcu_weights;
    bar_width = 0.3,
    color = :grey,
    label = "pcu",
    linecolor = :transparent
)

Plots.xlabel!("Coordination number")
Plots.ylabel!("Frequency")
Plots.xticks!(1:9)

# set the figsize
Plots.plot!(size = (450, 300))

Plots.savefig(path*"stern_ama_pcu_coordination_number_histogram.pdf")


bcu_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\crystals\\"

stern_vir_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\networks\stern_vir\\"

bcu_spatial_network = NG.load_spatial_network_from_gml(bcu_path*"bcu_cn_5_6_7_8_432_vertices.gml")

stern_vir_spatial_network = NG.load_spatial_network_from_gml(stern_vir_path*"stern_vir_blue.gml")

stern_vir_green_spatial_network = NG.load_spatial_network_from_gml(stern_vir_path*"stern_vir_green.gml")

max_coord_nr = 9

bcu_histogram = get_coord_nr_histogram(
    bcu_spatial_network;
    periodic_boundary_conditions=true,
    exclude_layer_thickness=0.0,
    max_coord_nr=9)

stern_vir_histogram = get_coord_nr_histogram(
    stern_vir_spatial_network;
    periodic_boundary_conditions=false,
    exclude_layer_thickness=1.5,
    max_coord_nr=max_coord_nr)

stern_vir_green_histogram = get_coord_nr_histogram(
    stern_vir_green_spatial_network;
    periodic_boundary_conditions=false,
    exclude_layer_thickness=1.5,
    max_coord_nr=max_coord_nr)


colors = Plots.palette(:tab10)

edges = bcu_histogram.edges[1][2:end-1] .+ 1
bcu_weights = bcu_histogram.weights[2:end]
stern_vir_weights = stern_vir_histogram.weights[2:end]
stern_vir_green_weights = stern_vir_green_histogram.weights[2:end]

# Define a small offset for side-by-side bars
offset = 0.2  # adjust based on bar width and spacing
total_offset = 0.5

# Plot stern_vir bars shifted left
Plots.bar(
    edges .- offset .- total_offset,
    stern_vir_green_weights;
    bar_width = 0.2,        # narrower bars
    color = colors[3],
    linecolor = :transparent,
    legend=false
)

# Overlay bcu bars shifted right
Plots.bar!(
    edges .- total_offset,
    stern_vir_weights;
    bar_width = 0.2,
    color = colors[1],
    linecolor = :transparent
)

# Overlay bcu bars shifted right
Plots.bar!(
    edges .+ offset .- total_offset,
    bcu_weights;
    bar_width = 0.2,
    color = :grey,
    linecolor = :transparent
)

Plots.xlabel!("Coordination number")
Plots.ylabel!("Frequency")
Plots.xticks!(1:max_coord_nr)

# set the figsize
Plots.plot!(size = (450, 300))

Plots.savefig(path*"stern_vir_bcu_coordination_number_histogram.pdf")


plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\\"

filename = "all_order_metrics.h5"

network_types = ["ctn", "dia", "lcs", "srs", "bcu_cn_5_6_7_8", "pcu_cn_4_5_6"]

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
        colorbar_title = "\n "*Latex.L"h_\mathrm{bond\ orientation}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\Delta T / T_\mathrm{melt}",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_tgradient_bond_orientation_entropy_beta_greater_0_less_10.pdf")

end


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
        colorbar_title = "\n "*Latex.L"\sigma_\mathrm{length}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_bond_length_std.pdf")

end


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

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = bond_angle_std_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (bond_angle_std_min, bond_angle_std_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"\sigma_\mathrm{angle}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_bond_angle_std.pdf")

end


filename = "all_order_metrics_with_nr_accepted_moves.h5"

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
        colorbar_title = "\n Nr. of accepted moves",
        colorbar_title_location = :right,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
        topmargin = 3Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_nr_accepted_moves.pdf")

end


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

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    t_gradient_normalized_vec_filtered = t_gradient_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    nr_accepted_moves_vec_filtered = nr_accepted_moves_vec[filter_indices]

    nr_accepted_moves_min = minimum(nr_accepted_moves_vec_filtered)
    nr_accepted_moves_max = maximum(nr_accepted_moves_vec_filtered)
    nr_accepted_moves_normalized_vec = (nr_accepted_moves_vec_filtered .- nr_accepted_moves_min) ./ (nr_accepted_moves_max - nr_accepted_moves_min)


    xlims=(0.5, 2.0)
    ylims=(0.25, 2.0)

    z = Float64.(copy(nr_accepted_moves_vec_filtered))
    z[z .< 10.0] .= NaN                      # mark zeros as NaN

    # Positive-only clims for log scale
    zmin = minimum(skipmissing(z[.!isnan.(z)]))  # smallest positive
    zmax = maximum(skipmissing(z[.!isnan.(z)]))

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        t_gradient_normalized_vec_filtered,
        zcolor = z,
        c = :bluesreds,
        clims = (10.0, 1000.0),
        colorbar = true,
        colorbar_scale = :log10,         # logarithmic colorbar for positive values
        nan_color = :black,              # <-- zeros (set to NaN) shown in black
        colorbar_title = "\n Nr. of accepted moves",
        colorbar_title_location = :right,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\Delta T / T_\mathrm{melt}",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
        topmargin = 3Plots.mm
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_tgradient_nr_accepted_moves.pdf")

end


for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    hyperuniformity_alpha_vec = order_metrics_generated_dict["hyperuniformity_alpha_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    hyperuniformity_alpha_vec_filtered = hyperuniformity_alpha_vec[filter_indices]

    hyperuniformity_alpha_min = minimum(hyperuniformity_alpha_vec_filtered)
    hyperuniformity_alpha_max = maximum(hyperuniformity_alpha_vec_filtered)
    hyperuniformity_alpha_normalized_vec = (hyperuniformity_alpha_vec_filtered .- hyperuniformity_alpha_min) ./ (hyperuniformity_alpha_max - hyperuniformity_alpha_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = hyperuniformity_alpha_vec_filtered,  # <-- use original (non-normalized) values
        c = :redsblues,                         # <-- same gradient you used
        clims = (-0.5, 0.5),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"\alpha",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_hyperuniformity_alpha.pdf")

end


for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    critical_pore_radius_vec = order_metrics_generated_dict["critical_pore_radius_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    critical_pore_radius_vec_filtered = critical_pore_radius_vec[filter_indices]

    critical_pore_radius_min = minimum(critical_pore_radius_vec_filtered)
    critical_pore_radius_max = maximum(critical_pore_radius_vec_filtered)
    critical_pore_radius_normalized_vec = (critical_pore_radius_vec_filtered .- critical_pore_radius_min) ./ (critical_pore_radius_max - critical_pore_radius_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = critical_pore_radius_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (critical_pore_radius_min, 0.62*critical_pore_radius_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n \n "*Latex.L"\delta_\mathrm{c}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 17Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_critical_pore_radius.pdf")

end


network_types = ["bcu_cn_5_6_7_8", "pcu_cn_4_5_6", "lcs" ]

for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    critical_pore_radius_vec = order_metrics_generated_dict["critical_pore_radius_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    critical_pore_radius_vec_filtered = critical_pore_radius_vec[filter_indices]

    critical_pore_radius_min = minimum(critical_pore_radius_vec_filtered)
    critical_pore_radius_max = maximum(critical_pore_radius_vec_filtered)
    critical_pore_radius_normalized_vec = (critical_pore_radius_vec_filtered .- critical_pore_radius_min) ./ (critical_pore_radius_max - critical_pore_radius_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = critical_pore_radius_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (critical_pore_radius_min, 0.62*critical_pore_radius_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"\delta_\mathrm{c}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_critical_pore_radius.pdf")

end


for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    ring_radius_std_vec = order_metrics_generated_dict["ring_radius_std_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    ring_radius_std_vec_filtered = ring_radius_std_vec[filter_indices]

    ring_radius_std_min = minimum(ring_radius_std_vec_filtered)
    ring_radius_std_max = maximum(ring_radius_std_vec_filtered)
    ring_radius_std_normalized_vec = (ring_radius_std_vec_filtered .- ring_radius_std_min) ./ (ring_radius_std_max - ring_radius_std_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = ring_radius_std_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (ring_radius_std_min, 0.75*ring_radius_std_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"\sigma_{r_\mathrm{ring}}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_ring_radius_std.pdf")
end

network_types = ["bcu_cn_5_6_7_8"]

for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    ring_radius_std_vec = order_metrics_generated_dict["ring_radius_std_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    ring_radius_std_vec_filtered = ring_radius_std_vec[filter_indices]

    ring_radius_std_min = minimum(ring_radius_std_vec_filtered)
    ring_radius_std_max = maximum(ring_radius_std_vec_filtered)
    ring_radius_std_normalized_vec = (ring_radius_std_vec_filtered .- ring_radius_std_min) ./ (ring_radius_std_max - ring_radius_std_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = ring_radius_std_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (ring_radius_std_min, 0.75*ring_radius_std_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n\n "*Latex.L"\sigma_{r_\mathrm{ring}}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 17Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_ring_radius_std.pdf")
end


for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    uncoordinated_neighbor_distance_vec = order_metrics_generated_dict["uncoordinated_neighbor_distance_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    uncoordinated_neighbor_distance_vec_filtered = uncoordinated_neighbor_distance_vec[filter_indices]

    uncoordinated_neighbor_distance_min = minimum(uncoordinated_neighbor_distance_vec_filtered)
    uncoordinated_neighbor_distance_max = maximum(uncoordinated_neighbor_distance_vec_filtered)
    uncoordinated_neighbor_distance_normalized_vec = (uncoordinated_neighbor_distance_vec_filtered .- uncoordinated_neighbor_distance_min) ./ (uncoordinated_neighbor_distance_max - uncoordinated_neighbor_distance_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = uncoordinated_neighbor_distance_vec_filtered,  # <-- use original (non-normalized) values
        c = :redsblues,                         # <-- same gradient you used
        clims = (uncoordinated_neighbor_distance_min, uncoordinated_neighbor_distance_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"R_\mathrm{uncoord}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_uncoordinated_neighbor_distance.pdf")

end


for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    vertex_homogeneity_metric_vec = order_metrics_generated_dict["vertex_homogeneity_metric_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    vertex_homogeneity_metric_vec_filtered = vertex_homogeneity_metric_vec[filter_indices]

    vertex_homogeneity_metric_min = minimum(vertex_homogeneity_metric_vec_filtered)
    vertex_homogeneity_metric_max = maximum(vertex_homogeneity_metric_vec_filtered)
    vertex_homogeneity_metric_normalized_vec = (vertex_homogeneity_metric_vec_filtered .- vertex_homogeneity_metric_min) ./ (vertex_homogeneity_metric_max - vertex_homogeneity_metric_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = vertex_homogeneity_metric_vec_filtered,  # <-- use original (non-normalized) values
        c = :redsblues,                         # <-- same gradient you used
        clims = (vertex_homogeneity_metric_min, vertex_homogeneity_metric_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"R_\mathrm{homogeneity}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_vertex_homogeneity_metric.pdf")

end

network_types = ["srs"]

for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    vertex_homogeneity_metric_vec = order_metrics_generated_dict["vertex_homogeneity_metric_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    # get the filter indices where t_gradient_normalized_vec is between 0.25 and 2.0
    # and t_max_normalized_vec is between 0.5 and 2.0
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))

    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    vertex_homogeneity_metric_vec_filtered = vertex_homogeneity_metric_vec[filter_indices]

    vertex_homogeneity_metric_min = minimum(vertex_homogeneity_metric_vec_filtered)
    vertex_homogeneity_metric_max = maximum(vertex_homogeneity_metric_vec_filtered)
    vertex_homogeneity_metric_normalized_vec = (vertex_homogeneity_metric_vec_filtered .- vertex_homogeneity_metric_min) ./ (vertex_homogeneity_metric_max - vertex_homogeneity_metric_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = vertex_homogeneity_metric_vec_filtered,  # <-- use original (non-normalized) values
        c = :redsblues,                         # <-- same gradient you used
        clims = (vertex_homogeneity_metric_min, vertex_homogeneity_metric_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n\n "*Latex.L"R_\mathrm{homogeneity}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 17Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_vertex_homogeneity_metric.pdf")

end


for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    anisotropy_metric_from_structure_factor_bonds_vec = order_metrics_generated_dict["anisotropy_metric_from_structure_factor_bonds_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the bond length std


    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    anisotropy_metric_from_structure_factor_bonds_vec_filtered = anisotropy_metric_from_structure_factor_bonds_vec[filter_indices]

    anisotropy_metric_from_structure_factor_bonds_min = minimum(anisotropy_metric_from_structure_factor_bonds_vec_filtered)
    anisotropy_metric_from_structure_factor_bonds_max = maximum(anisotropy_metric_from_structure_factor_bonds_vec_filtered)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)


    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = anisotropy_metric_from_structure_factor_bonds_vec_filtered,  # <-- use original (non-normalized) values
        c = :redsblues,                         # <-- same gradient you used
        clims = (anisotropy_metric_from_structure_factor_bonds_min, anisotropy_metric_from_structure_factor_bonds_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n \n "*Latex.L"A_\mathrm{bonds}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 17Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_anisotropy_metric_from_structure_factor_bonds.pdf")

end

network_types = ["bcu_cn_5_6_7_8", "pcu_cn_4_5_6", "dia", "srs"]

for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    anisotropy_metric_from_structure_factor_bonds_vec = order_metrics_generated_dict["anisotropy_metric_from_structure_factor_bonds_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the bond length std


    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    anisotropy_metric_from_structure_factor_bonds_vec_filtered = anisotropy_metric_from_structure_factor_bonds_vec[filter_indices]

    anisotropy_metric_from_structure_factor_bonds_min = minimum(anisotropy_metric_from_structure_factor_bonds_vec_filtered)
    anisotropy_metric_from_structure_factor_bonds_max = maximum(anisotropy_metric_from_structure_factor_bonds_vec_filtered)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)


    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = anisotropy_metric_from_structure_factor_bonds_vec_filtered,  # <-- use original (non-normalized) values
        c = :redsblues,                         # <-- same gradient you used
        clims = (anisotropy_metric_from_structure_factor_bonds_min, anisotropy_metric_from_structure_factor_bonds_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"A_\mathrm{bonds}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_anisotropy_metric_from_structure_factor_bonds.pdf")

end


for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    bond_orientation_entropy_vec = order_metrics_generated_dict["bond_orientation_entropy_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the bond length std


    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    bond_orientation_entropy_vec_filtered = bond_orientation_entropy_vec[filter_indices]

    bond_orientation_entropy_min = minimum(bond_orientation_entropy_vec_filtered)
    bond_orientation_entropy_max = maximum(bond_orientation_entropy_vec_filtered)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)


    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = bond_orientation_entropy_vec_filtered,  # <-- use original (non-normalized) values
        c = :bluesreds,                         # <-- same gradient you used
        clims = (bond_orientation_entropy_min, bond_orientation_entropy_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"h_\mathrm{bond\ orientation}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 12Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_bond_orientation_entropy.pdf")

end


for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    ring_radius_mean_vec = order_metrics_generated_dict["ring_radius_mean_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the bond length std

    # filter all indices where the t_gradient is between 0.25 and 2.0
    # and where t_max is between 0.5 and 2.0
    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    ring_radius_mean_vec_filtered = ring_radius_mean_vec[filter_indices]

    
    # first, convert the bond length std to a color scale between 0 and 1
    ring_radius_mean_min = minimum(ring_radius_mean_vec_filtered)
    ring_radius_mean_max = maximum(ring_radius_mean_vec_filtered)
    ring_radius_mean_normalized_vec = (ring_radius_mean_vec_filtered .- ring_radius_mean_min) ./ (ring_radius_mean_max - ring_radius_mean_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = ring_radius_mean_vec_filtered,  # <-- use original (non-normalized) values
        c = :redsblues,                         # <-- same gradient you used
        clims = (ring_radius_mean_min, ring_radius_mean_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n \n "*Latex.L"\overline{r}_\mathrm{ring}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 17Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_ring_radius_mean.pdf")

end

network_types = ["bcu_cn_5_6_7_8", "pcu_cn_4_5_6", "dia", "srs"]

for network_type in network_types
    
    order_metrics_generated_dict = GU.load_h5_dict(analysis_data_path*network_type*"\\"*filename)

    t_gradient_vec = order_metrics_generated_dict["t_gradient_vec"]
    t_max_vec = order_metrics_generated_dict["t_max_vec"]
    bond_bending_const_vec = order_metrics_generated_dict["bond_bending_const_vec"]

    t_melt_vec = [NA.get_melting_temperature(network_type, beta) for beta in bond_bending_const_vec]

    ring_radius_mean_vec = order_metrics_generated_dict["ring_radius_mean_vec"]

    #normalize temperatures with the melting temperature
    t_max_normalized_vec = t_max_vec ./ t_melt_vec
    t_gradient_normalized_vec = t_gradient_vec ./ t_melt_vec

    # create a scatter plot t_gradient against t_max where the color is given by the bond length std

    # filter all indices where the t_gradient is between 0.25 and 2.0
    # and where t_max is between 0.5 and 2.0
    filter_indices = findall(x -> x >= 0.25 && x <= 2.0, t_gradient_normalized_vec)
    filter_indices = intersect(filter_indices, findall(x -> x >= 0.5 && x <= 2.0, t_max_normalized_vec))
    t_max_normalized_vec_filtered = t_max_normalized_vec[filter_indices]
    bond_bending_const_vec_filtered = bond_bending_const_vec[filter_indices]
    ring_radius_mean_vec_filtered = ring_radius_mean_vec[filter_indices]

    
    # first, convert the bond length std to a color scale between 0 and 1
    ring_radius_mean_min = minimum(ring_radius_mean_vec_filtered)
    ring_radius_mean_max = maximum(ring_radius_mean_vec_filtered)
    ring_radius_mean_normalized_vec = (ring_radius_mean_vec_filtered .- ring_radius_mean_min) ./ (ring_radius_mean_max - ring_radius_mean_min)

    xlims=(0.5, 2.0)
    ylims=(minimum(bond_bending_const_vec)-0.1, maximum(bond_bending_const_vec)+0.1)

    p = Plots.scatter(
        t_max_normalized_vec_filtered,
        bond_bending_const_vec_filtered,
        zcolor = ring_radius_mean_vec_filtered,  # <-- use original (non-normalized) values
        c = :redsblues,                         # <-- same gradient you used
        clims = (ring_radius_mean_min, ring_radius_mean_max),  # <-- map correctly
        colorbar = true,                        # <-- shows on the right by default
        colorbar_title = "\n "*Latex.L"\overline{r}_\mathrm{ring}",
        colorbar_title_location = :right,  
        #colorbar_width = 0.01,
        markersize = 4,
        alpha = 0.6,
        xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
        ylabel = Latex.L"\beta",
        label = false,
        xlims = xlims,
        ylims = ylims,
        rightmargin = 10Plots.mm,
    )

    # set the figsize
    Plots.plot!(size = (500, 350))
    Plots.savefig(p, plot_path*network_type*"_tmax_beta_ring_radius_mean.pdf")

end



bio_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"
bcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\bcu_cn_5_6_7_8\\"
pcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\pcu_cn_4_5_6\\"
ctn_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\\"


analysis_data_paths = [bio_order_metrics_path, bcu_order_metrics_path]

filenames = ["stern_vir\\stern_vir_green_order_metrics.h5", "run_4\\bcu_cn_5_6_7_8_beta_0.0018_t_max_0.2401_t_gradient_0.3271_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    order_metrics_dict["coordination_nr_mean"] /= 10.0

    push!(dicts, order_metrics_dict)
end

all_order_metrics_dict = GU.load_h5_dict(bcu_order_metrics_path * "all_order_metrics_with_nr_accepted_moves.h5")
all_order_metrics_dict["coordination_nr_mean_vec"] ./= 10.0

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "uncoordinated_neighbor_distance",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor_bonds",
    "bond_orientation_entropy",
    "coordination_nr_std",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

# get all indices where the nr_accepted_moves_vec is greater than 10 
indices_accepted = findall(x -> x > 9, all_order_metrics_dict["nr_accepted_moves_vec"])

min_values = [minimum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]
max_values = [maximum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]

labels = [
    Latex.L"h_\mathrm{dihedral}",
    Latex.L"\sigma_\mathrm{angle}",
    Latex.L"\sigma_\mathrm{length}",
    Latex.L"\alpha",
    Latex.L"R_\mathrm{uncoord}",
    Latex.L"R_\mathrm{homogeneity}",
    Latex.L"\delta_\mathrm{c}",
    Latex.L"A_\mathrm{bonds}",
    Latex.L"h_\mathrm{bond\ orientation}",
    Latex.L"\sigma_{Z}",
    Latex.L"\overline{Z} / 10",
    Latex.L"\sigma_{r_\mathrm{ring}}",
    Latex.L"\overline{r}_\mathrm{ring}",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7
n3 = 9


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
min_values1 = min_values[1:n1]
max_values1 = max_values[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
min_values2 = min_values[n1+1:n2]
max_values2 = max_values[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:n3]
key_lists3 = key_list[n2+1:n3]
min_values3 = min_values[n2+1:n3]
max_values3 = max_values[n2+1:n3]
yvals3 = 1:length(labels3)

labels4 = labels[n3+1:end]
key_lists4 = key_list[n3+1:end]
min_values4 = min_values[n3+1:end]
max_values4 = max_values[n3+1:end]
yvals4 = 1:length(labels4)

# Calculate relative heights based on number of labels
total_labels = length(labels)



# Helper to add horizontal range bars [min, max] at each y in yvals
function add_range_bars!(plt, mins, maxs, yvals; color=:black, lw=2, offset=0.0)
    @assert length(mins) == length(maxs) == length(yvals)
    for i in eachindex(yvals)
        y = yvals[i] + offset  # offset can shift the bar slightly below the marker
        Plots.plot!(plt, [mins[i], maxs[i]], [y, y]; color=color, lw=lw, label=false)
    end
end


xlims = (-0.25, 1.25)
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
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p4 = Plots.plot(
    legend = false,
    yticks = (yvals4, labels4),
    xticks = xticks,
    xlabel = "Order metric value",
    grid = true
)


# Add horizontal bars first (so they appear underneath markers)
add_range_bars!(p1, min_values1, max_values1, yvals1; color=:black, lw=4, offset=0.0)
add_range_bars!(p2, min_values2, max_values2, yvals2; color=:black, lw=4, offset=0.0)
add_range_bars!(p3, min_values3, max_values3, yvals3; color=:black, lw=4, offset=0.0)
add_range_bars!(p4, min_values4, max_values4, yvals4; color=:black, lw=4, offset=0.0)


# Plot the scatter series for each subplot
markersizes = [12, 8]
color_palette = Plots.palette(:tab10)
colors = [color_palette[3], color_palette[9]]
for (i, dict) in enumerate(dicts)
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=markersizes[i], ylims=(0.5, n1+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=markersizes[i], ylims=(0.5, length(labels2)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p3, [Measurements.value.(dict[k]) for k in key_lists3], yvals3; label=labels3, markersize=markersizes[i], ylims=(0.5, length(labels3)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p4, [Measurements.value.(dict[k]) for k in key_lists4], yvals4; label=labels4, markersize=markersizes[i], ylims=(0.5, length(labels4)+0.5), xlims=xlims, color=colors[i])
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, p4, layout = (4, 1), size = (500, 700), link = :x,
bottom_margin = 0Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 10Plots.mm,
    right_margin = 0Plots.mm)

Plots.savefig(plot_path * "stern_vir_green_generated_order_metrics_comparison.pdf")



analysis_data_paths = [bio_order_metrics_path, bcu_order_metrics_path]

filenames = ["stern_vir\\stern_vir_blue_order_metrics.h5", "run_4\\bcu_cn_5_6_7_8_beta_0.0012_t_max_0.2314_t_gradient_0.4296_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    order_metrics_dict["coordination_nr_mean"] /= 10.0

    push!(dicts, order_metrics_dict)
end

all_order_metrics_dict = GU.load_h5_dict(bcu_order_metrics_path * "all_order_metrics_with_nr_accepted_moves.h5")
all_order_metrics_dict["coordination_nr_mean_vec"] ./= 10.0

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "uncoordinated_neighbor_distance",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor_bonds",
    "bond_orientation_entropy",
    "coordination_nr_std",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

# get all indices where the nr_accepted_moves_vec is greater than 10 
indices_accepted = findall(x -> x > 9, all_order_metrics_dict["nr_accepted_moves_vec"])

min_values = [minimum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]
max_values = [maximum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]

labels = [
    Latex.L"h_\mathrm{dihedral}",
    Latex.L"\sigma_\mathrm{angle}",
    Latex.L"\sigma_\mathrm{length}",
    Latex.L"\alpha",
    Latex.L"R_\mathrm{uncoord}",
    Latex.L"R_\mathrm{homogeneity}",
    Latex.L"\delta_\mathrm{c}",
    Latex.L"A_\mathrm{bonds}",
    Latex.L"h_\mathrm{bond\ orientation}",
    Latex.L"\sigma_{Z}",
    Latex.L"\overline{Z} / 10",
    Latex.L"\sigma_{r_\mathrm{ring}}",
    Latex.L"\overline{r}_\mathrm{ring}",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7
n3 = 9


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
min_values1 = min_values[1:n1]
max_values1 = max_values[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
min_values2 = min_values[n1+1:n2]
max_values2 = max_values[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:n3]
key_lists3 = key_list[n2+1:n3]
min_values3 = min_values[n2+1:n3]
max_values3 = max_values[n2+1:n3]
yvals3 = 1:length(labels3)

labels4 = labels[n3+1:end]
key_lists4 = key_list[n3+1:end]
min_values4 = min_values[n3+1:end]
max_values4 = max_values[n3+1:end]
yvals4 = 1:length(labels4)

# Calculate relative heights based on number of labels
total_labels = length(labels)



# Helper to add horizontal range bars [min, max] at each y in yvals
function add_range_bars!(plt, mins, maxs, yvals; color=:black, lw=2, offset=0.0)
    @assert length(mins) == length(maxs) == length(yvals)
    for i in eachindex(yvals)
        y = yvals[i] + offset  # offset can shift the bar slightly below the marker
        Plots.plot!(plt, [mins[i], maxs[i]], [y, y]; color=color, lw=lw, label=false)
    end
end


xlims = (-0.25, 1.25)
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
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p4 = Plots.plot(
    legend = false,
    yticks = (yvals4, labels4),
    xticks = xticks,
    xlabel = "Order metric value",
    grid = true
)


# Add horizontal bars first (so they appear underneath markers)
add_range_bars!(p1, min_values1, max_values1, yvals1; color=:black, lw=4, offset=0.0)
add_range_bars!(p2, min_values2, max_values2, yvals2; color=:black, lw=4, offset=0.0)
add_range_bars!(p3, min_values3, max_values3, yvals3; color=:black, lw=4, offset=0.0)
add_range_bars!(p4, min_values4, max_values4, yvals4; color=:black, lw=4, offset=0.0)


# Plot the scatter series for each subplot
markersizes = [12, 8]
color_palette = Plots.palette(:tab10)
colors = [color_palette[1], color_palette[10]]
for (i, dict) in enumerate(dicts)
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=markersizes[i], ylims=(0.5, n1+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=markersizes[i], ylims=(0.5, length(labels2)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p3, [Measurements.value.(dict[k]) for k in key_lists3], yvals3; label=labels3, markersize=markersizes[i], ylims=(0.5, length(labels3)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p4, [Measurements.value.(dict[k]) for k in key_lists4], yvals4; label=labels4, markersize=markersizes[i], ylims=(0.5, length(labels4)+0.5), xlims=xlims, color=colors[i])
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, p4, layout = (4, 1), size = (500, 700), link = :x,
bottom_margin = 0Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 10Plots.mm,
    right_margin = 0Plots.mm)

Plots.savefig(plot_path * "stern_vir_blue_generated_order_metrics_comparison.pdf")




analysis_data_paths = [bio_order_metrics_path, bio_order_metrics_path, bcu_order_metrics_path, bcu_order_metrics_path]

filenames = ["stern_vir\\stern_vir_green_order_metrics.h5", "stern_vir\\stern_vir_blue_order_metrics.h5",  
"run_4\\bcu_cn_5_6_7_8_beta_0.0018_t_max_0.2401_t_gradient_0.3271_order_metrics.h5",
"run_4\\bcu_cn_5_6_7_8_beta_0.0012_t_max_0.2314_t_gradient_0.4296_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    order_metrics_dict["coordination_nr_mean"] /= 10.0

    push!(dicts, order_metrics_dict)
end

all_order_metrics_dict = GU.load_h5_dict(bcu_order_metrics_path * "all_order_metrics_with_nr_accepted_moves.h5")
all_order_metrics_dict["coordination_nr_mean_vec"] ./= 10.0

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "uncoordinated_neighbor_distance",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor_bonds",
    "bond_orientation_entropy",
    "coordination_nr_std",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

# get all indices where the nr_accepted_moves_vec is greater than 10 
indices_accepted = findall(x -> x > 9, all_order_metrics_dict["nr_accepted_moves_vec"])

min_values = [minimum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]
max_values = [maximum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]

labels = [
    Latex.L"h_\mathrm{dihedral}",
    Latex.L"\sigma_\mathrm{angle}",
    Latex.L"\sigma_\mathrm{length}",
    Latex.L"\alpha",
    Latex.L"R_\mathrm{uncoord}",
    Latex.L"R_\mathrm{homogeneity}",
    Latex.L"\delta_\mathrm{c}",
    Latex.L"A_\mathrm{bonds}",
    Latex.L"h_\mathrm{bond\ orientation}",
    Latex.L"\sigma_{Z}",
    Latex.L"\overline{Z} / 10",
    Latex.L"\sigma_{r_\mathrm{ring}}",
    Latex.L"\overline{r}_\mathrm{ring}",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7
n3 = 9


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
min_values1 = min_values[1:n1]
max_values1 = max_values[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
min_values2 = min_values[n1+1:n2]
max_values2 = max_values[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:n3]
key_lists3 = key_list[n2+1:n3]
min_values3 = min_values[n2+1:n3]
max_values3 = max_values[n2+1:n3]
yvals3 = 1:length(labels3)

labels4 = labels[n3+1:end]
key_lists4 = key_list[n3+1:end]
min_values4 = min_values[n3+1:end]
max_values4 = max_values[n3+1:end]
yvals4 = 1:length(labels4)

# Calculate relative heights based on number of labels
total_labels = length(labels)



# Helper to add horizontal range bars [min, max] at each y in yvals
function add_range_bars!(plt, mins, maxs, yvals; color=:black, lw=2, offset=0.0)
    @assert length(mins) == length(maxs) == length(yvals)
    for i in eachindex(yvals)
        y = yvals[i] + offset  # offset can shift the bar slightly below the marker
        Plots.plot!(plt, [mins[i], maxs[i]], [y, y]; color=color, lw=lw, label=false)
    end
end


xlims = (-0.25, 1.25)
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
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p4 = Plots.plot(
    legend = false,
    yticks = (yvals4, labels4),
    xticks = xticks,
    xlabel = "Order metric value",
    grid = true
)


# Add horizontal bars first (so they appear underneath markers)
add_range_bars!(p1, min_values1, max_values1, yvals1; color=:black, lw=4, offset=0.0)
add_range_bars!(p2, min_values2, max_values2, yvals2; color=:black, lw=4, offset=0.0)
add_range_bars!(p3, min_values3, max_values3, yvals3; color=:black, lw=4, offset=0.0)
add_range_bars!(p4, min_values4, max_values4, yvals4; color=:black, lw=4, offset=0.0)


# Plot the scatter series for each subplot
markersizes = [12, 12, 8, 8]
color_palette = Plots.palette(:tab10)
colors = [color_palette[3], color_palette[1], color_palette[9], color_palette[10]]
alphas = [1.0, 0.8, 0.8, 0.8]
for (i, dict) in enumerate(dicts)
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=markersizes[i], ylims=(0.5, n1+0.5), xlims=xlims, color=colors[i], alpha=alphas[i])

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=markersizes[i], ylims=(0.5, length(labels2)+0.5), xlims=xlims, color=colors[i], alpha=alphas[i])
    Plots.scatter!(p3, [Measurements.value.(dict[k]) for k in key_lists3], yvals3; label=labels3, markersize=markersizes[i], ylims=(0.5, length(labels3)+0.5), xlims=xlims, color=colors[i], alpha=alphas[i])

    Plots.scatter!(p4, [Measurements.value.(dict[k]) for k in key_lists4], yvals4; label=labels4, markersize=markersizes[i], ylims=(0.5, length(labels4)+0.5), xlims=xlims, color=colors[i], alpha=alphas[i])
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, p4, layout = (4, 1), size = (500, 700), link = :x,
bottom_margin = 0Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 10Plots.mm,
    right_margin = 0Plots.mm)

Plots.savefig(plot_path * "stern_vir_generated_order_metrics_comparison.pdf")


analysis_data_paths = [bio_order_metrics_path, ctn_order_metrics_path]

filenames = ["pachy\\pachy_blue_order_metrics.h5", "run_5\\ctn_beta_6.3359_t_max_8.6070_t_gradient_11.7170_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    order_metrics_dict["coordination_nr_mean"] /= 10.0

    push!(dicts, order_metrics_dict)
end

all_order_metrics_dict = GU.load_h5_dict(ctn_order_metrics_path * "all_order_metrics_with_nr_accepted_moves.h5")
all_order_metrics_dict["coordination_nr_mean_vec"] ./= 10.0

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "uncoordinated_neighbor_distance",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor_bonds",
    "bond_orientation_entropy",
    "coordination_nr_std",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

# get all indices where the nr_accepted_moves_vec is greater than 10 
indices_accepted = findall(x -> x > 9, all_order_metrics_dict["nr_accepted_moves_vec"])

min_values = [minimum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]
max_values = [maximum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]

labels = [
    Latex.L"h_\mathrm{dihedral}",
    Latex.L"\sigma_\mathrm{angle}",
    Latex.L"\sigma_\mathrm{length}",
    Latex.L"\alpha",
    Latex.L"R_\mathrm{uncoord}",
    Latex.L"R_\mathrm{homogeneity}",
    Latex.L"\delta_\mathrm{c}",
    Latex.L"A_\mathrm{bonds}",
    Latex.L"h_\mathrm{bond\ orientation}",
    Latex.L"\sigma_{Z}",
    Latex.L"\overline{Z} / 10",
    Latex.L"\sigma_{r_\mathrm{ring}}",
    Latex.L"\overline{r}_\mathrm{ring}",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7
n3 = 9


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
min_values1 = min_values[1:n1]
max_values1 = max_values[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
min_values2 = min_values[n1+1:n2]
max_values2 = max_values[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:n3]
key_lists3 = key_list[n2+1:n3]
min_values3 = min_values[n2+1:n3]
max_values3 = max_values[n2+1:n3]
yvals3 = 1:length(labels3)

labels4 = labels[n3+1:end]
key_lists4 = key_list[n3+1:end]
min_values4 = min_values[n3+1:end]
max_values4 = max_values[n3+1:end]
yvals4 = 1:length(labels4)

# Calculate relative heights based on number of labels
total_labels = length(labels)



# Helper to add horizontal range bars [min, max] at each y in yvals
function add_range_bars!(plt, mins, maxs, yvals; color=:black, lw=2, offset=0.0)
    @assert length(mins) == length(maxs) == length(yvals)
    for i in eachindex(yvals)
        y = yvals[i] + offset  # offset can shift the bar slightly below the marker
        Plots.plot!(plt, [mins[i], maxs[i]], [y, y]; color=color, lw=lw, label=false)
    end
end


xlims = (-0.25, 1.25)
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
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p4 = Plots.plot(
    legend = false,
    yticks = (yvals4, labels4),
    xticks = xticks,
    xlabel = "Order metric value",
    grid = true
)


# Add horizontal bars first (so they appear underneath markers)
add_range_bars!(p1, min_values1, max_values1, yvals1; color=:black, lw=4, offset=0.0)
add_range_bars!(p2, min_values2, max_values2, yvals2; color=:black, lw=4, offset=0.0)
add_range_bars!(p3, min_values3, max_values3, yvals3; color=:black, lw=4, offset=0.0)
add_range_bars!(p4, min_values4, max_values4, yvals4; color=:black, lw=4, offset=0.0)


# Plot the scatter series for each subplot
markersizes = [12, 8]
color_palette = Plots.palette(:tab10)
colors = [color_palette[1], color_palette[10]]
for (i, dict) in enumerate(dicts)
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=markersizes[i], ylims=(0.5, n1+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=markersizes[i], ylims=(0.5, length(labels2)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p3, [Measurements.value.(dict[k]) for k in key_lists3], yvals3; label=labels3, markersize=markersizes[i], ylims=(0.5, length(labels3)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p4, [Measurements.value.(dict[k]) for k in key_lists4], yvals4; label=labels4, markersize=markersizes[i], ylims=(0.5, length(labels4)+0.5), xlims=xlims, color=colors[i])
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, p4, layout = (4, 1), size = (500, 700), link = :x,
bottom_margin = 0Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 10Plots.mm,
    right_margin = 0Plots.mm)

Plots.savefig(plot_path * "pachy_generated_order_metrics_comparison.pdf")


bio_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"
pcu_order_metrics_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\pcu_cn_4_5_6\\"

analysis_data_paths = [bio_order_metrics_path, pcu_order_metrics_path]

filenames = ["stern_ama\\stern_ama_orange_order_metrics.h5", "run_3\\pcu_cn_4_5_6_beta_0.0978_t_max_0.1788_t_gradient_0.1858_order_metrics.h5"]

dicts = []

for (i, filename) in enumerate(filenames)

    order_metrics_dict = GU.load_h5_dict(analysis_data_paths[i] * filename)
    order_metrics_dict["coordination_nr_mean"] /= 10.0

    push!(dicts, order_metrics_dict)
end

all_order_metrics_dict = GU.load_h5_dict(pcu_order_metrics_path * "all_order_metrics_with_nr_accepted_moves.h5")
all_order_metrics_dict["coordination_nr_mean_vec"] ./= 10.0

key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "uncoordinated_neighbor_distance",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor_bonds",
    "bond_orientation_entropy",
    "coordination_nr_std",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]

# get all indices where the nr_accepted_moves_vec is greater than 10 
indices_accepted = findall(x -> x > 9, all_order_metrics_dict["nr_accepted_moves_vec"])

min_values = [minimum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]
max_values = [maximum(all_order_metrics_dict[k*"_vec"][indices_accepted]) for k in key_list]

labels = [
    Latex.L"h_\mathrm{dihedral}",
    Latex.L"\sigma_\mathrm{angle}",
    Latex.L"\sigma_\mathrm{length}",
    Latex.L"\alpha",
    Latex.L"R_\mathrm{uncoord}",
    Latex.L"R_\mathrm{homogeneity}",
    Latex.L"\delta_\mathrm{c}",
    Latex.L"A_\mathrm{bonds}",
    Latex.L"h_\mathrm{bond\ orientation}",
    Latex.L"\sigma_{Z}",
    Latex.L"\overline{Z} / 10",
    Latex.L"\sigma_{r_\mathrm{ring}}",
    Latex.L"\overline{r}_\mathrm{ring}",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7
n3 = 9


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
min_values1 = min_values[1:n1]
max_values1 = max_values[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
min_values2 = min_values[n1+1:n2]
max_values2 = max_values[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:n3]
key_lists3 = key_list[n2+1:n3]
min_values3 = min_values[n2+1:n3]
max_values3 = max_values[n2+1:n3]
yvals3 = 1:length(labels3)

labels4 = labels[n3+1:end]
key_lists4 = key_list[n3+1:end]
min_values4 = min_values[n3+1:end]
max_values4 = max_values[n3+1:end]
yvals4 = 1:length(labels4)

# Calculate relative heights based on number of labels
total_labels = length(labels)



# Helper to add horizontal range bars [min, max] at each y in yvals
function add_range_bars!(plt, mins, maxs, yvals; color=:black, lw=2, offset=0.0)
    @assert length(mins) == length(maxs) == length(yvals)
    for i in eachindex(yvals)
        y = yvals[i] + offset  # offset can shift the bar slightly below the marker
        Plots.plot!(plt, [mins[i], maxs[i]], [y, y]; color=color, lw=lw, label=false)
    end
end


xlims = (-0.25, 1.25)
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
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p4 = Plots.plot(
    legend = false,
    yticks = (yvals4, labels4),
    xticks = xticks,
    xlabel = "Order metric value",
    grid = true
)


# Add horizontal bars first (so they appear underneath markers)
add_range_bars!(p1, min_values1, max_values1, yvals1; color=:black, lw=4, offset=0.0)
add_range_bars!(p2, min_values2, max_values2, yvals2; color=:black, lw=4, offset=0.0)
add_range_bars!(p3, min_values3, max_values3, yvals3; color=:black, lw=4, offset=0.0)
add_range_bars!(p4, min_values4, max_values4, yvals4; color=:black, lw=4, offset=0.0)


# Plot the scatter series for each subplot
markersizes = [12, 8]
color_palette = Plots.palette(:tab10)
colors = [color_palette[2], color_palette[4]]
for (i, dict) in enumerate(dicts)
    Plots.scatter!(p1, [Measurements.value.(dict[k]) for k in key_lists1], yvals1; label=labels1, markersize=markersizes[i], ylims=(0.5, n1+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p2, [Measurements.value.(dict[k]) for k in key_lists2], yvals2; label=labels2, markersize=markersizes[i], ylims=(0.5, length(labels2)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p3, [Measurements.value.(dict[k]) for k in key_lists3], yvals3; label=labels3, markersize=markersizes[i], ylims=(0.5, length(labels3)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p4, [Measurements.value.(dict[k]) for k in key_lists4], yvals4; label=labels4, markersize=markersizes[i], ylims=(0.5, length(labels4)+0.5), xlims=xlims, color=colors[i])
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, p4, layout = (4, 1), size = (500, 700), link = :x,
bottom_margin = 0Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 10Plots.mm,
    right_margin = 0Plots.mm)

Plots.savefig(plot_path * "stern_ama_generated_order_metrics_comparison.pdf")


function scatter_order_metrics(
    order_metrics_bio_dict,
    order_metrics_generated_dict;
    metric_x::String,
    metric_y::String,
    metric_color::String,
    x_label::AbstractString = "",
    y_label::AbstractString = "",
    xlim::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    ylim::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    bio_index::Int = 1,
    highlight_generated_index::Union{Nothing,Int} = nothing,
    generated_label::AbstractString = "Generated",
    bio_label::AbstractString = "Biological",
    clim = (10.0, 1000.0),
    bio_color = :blue,
    highlight_generated_color = :red,
    highlight_ms = 10,
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

    color_gen = Float64.(copy(order_metrics_generated_dict[metric_color]))
    color_gen[ color_gen .< 10] .= NaN
    color_gen_min = minimum(skipmissing(color_gen))
    color_gen_max = maximum(skipmissing(color_gen))

    x_bio = order_metrics_bio_dict[metric_x][bio_index]
    y_bio = order_metrics_bio_dict[metric_y][bio_index]

    # Base scatter for generated networks
    plt = Plots.scatter(
        x_gen, y_gen;
        label = generated_label,
        zcolor = color_gen,
        c = :bluesreds,
        clim = clim,
        colorbar = true,
        colorbar_scale = :log10,         # logarithmic colorbar for positive values
        nan_color = :black,              # <-- zeros (set to NaN) shown in black
        colorbar_title = "\n Nr. of accepted moves",
        marker = (generated_marker, generated_ms),
        alpha = alpha_gen,
        legend = legend,
        title = title,
        xlabel = x_label == "" ? metric_x : x_label,
        ylabel = y_label == "" ? metric_y : y_label,
        xlim = xlim,
        ylim = ylim,
        grid = :on,
        rightmargin = 10Plots.mm,
        topmargin = 3Plots.mm,
    )

    # Overlay the target bio point
    Plots.scatter!(
        plt,
        [x_bio], [y_bio];
        label = bio_label,
        color = bio_color,
        marker = (bio_marker, bio_ms),
        alpha = 0.8,
    )

    if highlight_generated_index !== nothing
        x_highlight = order_metrics_generated_dict[metric_x][highlight_generated_index]
        y_highlight = order_metrics_generated_dict[metric_y][highlight_generated_index]

        # Overlay the highlighted generated point
        Plots.scatter!(
            plt,
            [x_highlight], [y_highlight];
            label = "Best match",
            color = highlight_generated_color,
            marker = (generated_marker, highlight_ms),
            alpha = 0.8,
            legend=false
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


plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

bio_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"

generated_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\ctn\\"

order_metrics_bio_dict = GU.load_h5_dict(bio_path*"all_order_metrics.h5")

order_metrics_generated_dict = GU.load_h5_dict(generated_path*"all_order_metrics_with_nr_accepted_moves.h5")


metric_x = "hyperuniformity_alpha_vec"
metric_y = "bond_length_std_vec"
metric_color = "nr_accepted_moves_vec"

colors = Plots.palette(:tab10)

plt = scatter_order_metrics(
    order_metrics_bio_dict,
    order_metrics_generated_dict;
    metric_x = metric_x,
    metric_y = metric_y,
    metric_color = metric_color,
    x_label = "α",
    y_label = Latex.L"\sigma_\mathrm{length}",
    xlim = (-2.5, 4.0),
    bio_index = 1,
    highlight_generated_color = colors[10],
    bio_color = colors[1],
    bio_marker=:circle,
    bio_ms=11,
    generated_ms = 4,
    highlight_ms = 8,
    alpha_gen = 0.6,
    highlight_generated_index = 3091,
    savepath = plot_path * "pachy_generated_hyperuniformity_alpha_vs_bond_length_std.pdf"
)


metric_y = "bond_angle_std_vec"


plt = scatter_order_metrics(
    order_metrics_bio_dict,
    order_metrics_generated_dict;
    metric_x = metric_x,
    metric_y = metric_y,
    metric_color = metric_color,
    x_label = "α",
    y_label = Latex.L"\sigma_\mathrm{angle}",
    xlim = (-2.5, 4.0),
    bio_index = 1,
    highlight_generated_color = colors[10],
    bio_color = colors[1],
    bio_marker=:circle,
    bio_ms=11,
    generated_ms = 4,
    highlight_ms = 8,
    highlight_generated_index = 3091,
    savepath = plot_path * "pachy_generated_hyperuniformity_alpha_vs_bond_angle_std.pdf"
)




function scatter_order_metrics(
    order_metrics_bio_dict,
    order_metrics_generated_dict;
    metric_x::String,
    metric_y::String,
    metric_color::String,
    x_label::AbstractString = "",
    y_label::AbstractString = "",
    xlim::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    ylim::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    bio_index::Int = 1,
    highlight_generated_index::Union{Nothing,Int} = nothing,
    generated_label::AbstractString = "Generated",
    bio_label::AbstractString = "Biological",
    clim = (10.0, 1000.0),
    bio_color = :blue,
    highlight_generated_color = :red,
    highlight_ms = 10,
    generated_marker = :circle,
    bio_marker = :star5,
    generated_ms = 5,
    bio_ms = 12,
    bio_index_2::Union{Nothing,Int} = nothing,
    bio_color_2 = :green,
    alpha_gen = 0.8,
    legend = :topright,
    title::AbstractString = "",
    savepath::Union{Nothing,String} = nothing
)
    # Extract data (we rely on your guarantee that keys exist and sizes align)
    x_gen = order_metrics_generated_dict[metric_x]
    y_gen = order_metrics_generated_dict[metric_y]

    color_gen = Float64.(copy(order_metrics_generated_dict[metric_color]))
    color_gen[ color_gen .< 10] .= NaN
    color_gen_min = minimum(skipmissing(color_gen))
    color_gen_max = maximum(skipmissing(color_gen))

    x_bio = order_metrics_bio_dict[metric_x][bio_index]
    y_bio = order_metrics_bio_dict[metric_y][bio_index]

    # Base scatter for generated networks
    plt = Plots.scatter(
        x_gen, y_gen;
        label = generated_label,
        zcolor = color_gen,
        c = :bluesreds,
        clim = clim,
        colorbar = true,
        colorbar_scale = :log10,         # logarithmic colorbar for positive values
        nan_color = :black,              # <-- zeros (set to NaN) shown in black
        colorbar_title = "\n Nr. of accepted moves",
        marker = (generated_marker, generated_ms),
        alpha = alpha_gen,
        legend = legend,
        title = title,
        xlabel = x_label == "" ? metric_x : x_label,
        ylabel = y_label == "" ? metric_y : y_label,
        xlim = xlim,
        ylim = ylim,
        grid = :on,
        rightmargin = 10Plots.mm,
        topmargin = 3Plots.mm,
    )

    # Overlay the target bio point
    Plots.scatter!(
        plt,
        [x_bio], [y_bio];
        label = bio_label,
        color = bio_color,
        marker = (bio_marker, bio_ms),
        alpha = 0.8,
    )

    if bio_index_2 !== nothing
        x_bio_2 = order_metrics_bio_dict[metric_x][bio_index_2]
        y_bio_2 = order_metrics_bio_dict[metric_y][bio_index_2]

        # Overlay the second target bio point
        Plots.scatter!(
            plt,
            [x_bio_2], [y_bio_2];
            label = bio_label * " 2",
            color = bio_color_2,
            marker = (bio_marker, bio_ms),
            alpha = 0.8,
        )
    end

    if highlight_generated_index !== nothing
        x_highlight = order_metrics_generated_dict[metric_x][highlight_generated_index]
        y_highlight = order_metrics_generated_dict[metric_y][highlight_generated_index]

        # Overlay the highlighted generated point
        Plots.scatter!(
            plt,
            [x_highlight], [y_highlight];
            label = "Best match",
            color = highlight_generated_color,
            marker = (generated_marker, highlight_ms),
            alpha = 0.8,
            legend=false
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


plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

bio_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\\"

generated_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\random\bcu_cn_5_6_7_8\\"

order_metrics_bio_dict = GU.load_h5_dict(bio_path*"all_order_metrics.h5")

order_metrics_generated_dict = GU.load_h5_dict(generated_path*"all_order_metrics_with_nr_accepted_moves.h5")

# "pachy/pachy_blue"
# "pachy/pachy_red"
# "stern_ama/stern_ama_orange"
# "stern_vir/stern_vir_blue"
# "stern_vir/stern_vir_green"


metric_x = "hyperuniformity_alpha_vec"
metric_y = "bond_length_std_vec"
metric_color = "nr_accepted_moves_vec"

colors = Plots.palette(:tab10)

plt = scatter_order_metrics(
    order_metrics_bio_dict,
    order_metrics_generated_dict;
    metric_x = metric_x,
    metric_y = metric_y,
    metric_color = metric_color,
    x_label = "α",
    y_label = Latex.L"\sigma_\mathrm{length}",
    xlim = (-2.5, 4.0),
    bio_index = 4,
    bio_color = colors[1],
    bio_marker=:circle,
    bio_ms=11,
    bio_index_2 = 5,
    bio_color_2 = colors[3],
    generated_ms = 4,
    highlight_generated_color = colors[9],
    highlight_ms = 8,
    alpha_gen = 0.6,
    highlight_generated_index = 2,
    savepath = plot_path * "stern_vir_generated_hyperuniformity_alpha_vs_bond_length_std.pdf"
)


metric_y = "bond_angle_std_vec"


plt = scatter_order_metrics(
    order_metrics_bio_dict,
    order_metrics_generated_dict;
    metric_x = metric_x,
    metric_y = metric_y,
    metric_color = metric_color,
    x_label = "α",
    y_label = Latex.L"\sigma_\mathrm{angle}",
    xlim = (-2.5, 4.0),
    bio_index = 4,
    bio_color = colors[1],
    bio_marker=:circle,
    bio_ms=11,
    bio_index_2 = 5,
    bio_color_2 = colors[3],
    generated_ms = 4,
    highlight_generated_color = colors[9],
    highlight_ms = 8,
    alpha_gen = 0.6,
    highlight_generated_index = 2,
    savepath = plot_path * "stern_vir_generated_hyperuniformity_alpha_vs_bond_angle_std.pdf"
)



plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

network_type = "ctn"

predictions_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\predictions\local_relaxation\ctn\\"

prediction_filename = "ctn_predictions_nr_layers_4_nr_neurons_72_full_pca_10_reordered.h5"

predictions_dict = GU.load_h5_dict(predictions_path*prediction_filename)

# print shape of all keys in the predictions_dict
for (k, v) in predictions_dict
    println("Key: $k, Shape: ", size(v))
end

# t_max_over_t_melt_vec, Shape: (75,)
# Key: predictions_array, Shape: (75, 50, 42)
# Key: bond_bending_const_vec, Shape: (50,)

t_max_over_t_melt_vec = predictions_dict["t_max_over_t_melt_vec"]
bond_bending_const_vec = predictions_dict["bond_bending_const_vec"]
predictions_array = predictions_dict["predictions_array"]


Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 1]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n  "*Latex.L"\sigma_\mathrm{length}",
    c = :bluesreds, 
    clims = (0.01, 0.63),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 10Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_bond_length_std.pdf")

Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 2]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n  "*Latex.L"\sigma_\mathrm{angle}",
    c = :bluesreds, 
    clims = (0.12, 0.54),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 10Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_bond_angle_std.pdf")



Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 4]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n  "*Latex.L"h_\mathrm{bond\ orientation}",
    c = :bluesreds, 
    clims = (0.74, 0.999),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 10Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_bond_orientation_entropy.pdf")


Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 33]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n  "*Latex.L"R_\mathrm{homogeneity}",
    c = :redsblues, 
    clims = (0.59, 1.01),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 10Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_vertex_homogeneity_metric.pdf")

Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 34]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n  "*Latex.L"R_\mathrm{uncoord}",
    c = :redsblues, 
    clims = (0.59, 1.55),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 10Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_uncoordinated_neighbor_distance.pdf")


Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 37]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n \n  "*Latex.L"\overline{r}_\mathrm{ring}",
    c = :redsblues, 
    clims = (0.62, 0.85),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 17Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_ring_radius_mean.pdf")

Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 38]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n  "*Latex.L"\sigma_\mathrm{ring}",
    c = :bluesreds, 
    clims = (0.001, 0.091),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 10Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_ring_radius_std.pdf")


Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 39]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n \n "*Latex.L"\delta_c",
    c = :bluesreds, 
    clims = (0.298, 0.415),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 17Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_critical_pore_radius.pdf")

Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 41]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n \n "*Latex.L"A_\mathrm{bonds}",
    c = :redsblues, 
    clims = (0.447, 0.651),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 17Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_anisotropy_metric_from_structure_factor_bonds.pdf")


Plots.heatmap(
    t_max_over_t_melt_vec,
    bond_bending_const_vec,
    predictions_array[:, :, 42]', 
    ylabel = Latex.L"\beta",
    xlabel = Latex.L"T_\mathrm{max} / T_\mathrm{melt}",
    colorbar_title = "\n  "*Latex.L"\alpha",
    c = :redsblues, 
    clims = (-0.5, 0.5),
    xlims = (minimum(t_max_over_t_melt_vec), maximum(t_max_over_t_melt_vec)),
    ylims = (minimum(bond_bending_const_vec), maximum(bond_bending_const_vec)),
    xticks = collect(0.5:0.5:2.0),
    yticks = collect(0.0:2.5:10.0),
    right_margin = 10Plots.mm,
)

Plots.savefig(plot_path*network_type*"_predict_tmax_beta_hyperuniformity_alpha.pdf")


plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

network_type = "ctn"

predictions_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\predictions\local_relaxation\ctn\\"

r2_filename = "ctn_nr_layers_4_nr_neurons_72_r2_scores_pca_10.h5"

r2_dict = GU.load_h5_dict(predictions_path*r2_filename)


key_list = [
    "dihedral_angle_entropy",
    "bond_angle_std",
    "bond_length_std",
    "hyperuniformity_alpha",
    "uncoordinated_neighbor_distance",
    "vertex_homogeneity_metric",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor_bonds",
    "bond_orientation_entropy",
    "coordination_nr_std",
    "coordination_nr_mean",
    "ring_radius_std",
    "ring_radius_mean",
]


labels = [
    Latex.L"h_\mathrm{dihedral}",
    Latex.L"\sigma_\mathrm{angle}",
    Latex.L"\sigma_\mathrm{length}",
    Latex.L"\alpha",
    Latex.L"R_\mathrm{uncoord}",
    Latex.L"R_\mathrm{homogeneity}",
    Latex.L"\delta_\mathrm{c}",
    Latex.L"A_\mathrm{bonds}",
    Latex.L"h_\mathrm{bond\ orientation}",
    Latex.L"\sigma_{Z}",
    Latex.L"\overline{Z}",
    Latex.L"\sigma_{r_\mathrm{ring}}",
    Latex.L"\overline{r}_\mathrm{ring}",
]

# Map keys to numeric y positions
yvals = 1:length(key_list)

# Example split indices
n1 = 3
n2 = 7
n3 = 9


# Split the labels and corresponding yvals into three groups
labels1 = labels[1:n1]
key_lists1 = key_list[1:n1]
yvals1 = 1:n1

labels2 = labels[n1+1:n2]
key_lists2 = key_list[n1+1:n2]
yvals2 = 1:length(labels2)

labels3 = labels[n2+1:n3]
key_lists3 = key_list[n2+1:n3]
yvals3 = 1:length(labels3)

labels4 = labels[n3+1:end]
key_lists4 = key_list[n3+1:end]
yvals4 = 1:length(labels4)

# Calculate relative heights based on number of labels
total_labels = length(labels)



xlims = (0.57, 1.03)
xticks = [ 0.6, 0.8, 1.0]

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
    xticks = (xticks, ["", "", ""]),
    grid = true
)

p4 = Plots.plot(
    legend = false,
    yticks = (yvals4, labels4),
    xticks = xticks,
    xlabel = Latex.L"R^{2} \mathrm{\ score}",
    grid = true
)

dicts = [r2_dict]

# Plot the scatter series for each subplot
markersizes = [12]
color_palette = Plots.palette(:tab10)
colors = [:black]
for (i, dict) in enumerate(dicts)
    Plots.scatter!(p1, [dict[k*"_vec"] for k in key_lists1], yvals1; label=labels1, markersize=markersizes[i], ylims=(0.5, n1+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p2, [dict[k*"_vec"] for k in key_lists2], yvals2; label=labels2, markersize=markersizes[i], ylims=(0.5, length(labels2)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p3, [dict[k*"_vec"] for k in key_lists3], yvals3; label=labels3, markersize=markersizes[i], ylims=(0.5, length(labels3)+0.5), xlims=xlims, color=colors[i])

    Plots.scatter!(p4, [dict[k*"_vec"] for k in key_lists4], yvals4; label=labels4, markersize=markersizes[i], ylims=(0.5, length(labels4)+0.5), xlims=xlims, color=colors[i])
end


# Combine the subplots into a single plot
Plots.plot(p1, p2, p3, p4, layout = (4, 1), size = (500, 700), link = :x,
bottom_margin = 0Plots.mm,
    top_margin = 3Plots.mm,
    left_margin = 10Plots.mm,
    right_margin = 0Plots.mm)

Plots.savefig(plot_path * "ctn_nr_layers_4_nr_neurons_72_r2_scores_pca_10.pdf")



plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\writing\paper_disordered_networks\figs\\"

# Path to your PCA file
pca_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\datasets\local_relaxation\ctn\\"
pca_file = "ctn_all_order_metrics_metric_prediction_dataset_pca_10_pca_information.h5"

pca_dict = GU.load_h5_dict(pca_path*pca_file)

explained_variance_ratio = pca_dict["explained_variance_ratio_"]
components = pca_dict["components_"]

labels = [
    Latex.L"\sigma_\mathrm{length}",
    Latex.L"\sigma_\mathrm{angle}",
    Latex.L"h_\mathrm{dihedral}",
    Latex.L"h_\mathrm{bond\ orientation}",
    Latex.L"\overline{Z}",
    Latex.L"\sigma_{Z}",
    Latex.L"\overline{q}_{0}",
    Latex.L"\overline{q}_{1}",
    Latex.L"\overline{q}_{2}",
    Latex.L"\overline{q}_{3}",
    Latex.L"\overline{q}_{4}",
    Latex.L"\overline{q}_{5}",
    Latex.L"\overline{q}_{6}",
    Latex.L"\overline{q}_{7}",
    Latex.L"\overline{q}_{8}",
    Latex.L"\overline{q}_{9}",
    Latex.L"\overline{q}_{10}",
    Latex.L"\overline{q}_{11}",
    Latex.L"\overline{q}_{12}",
    Latex.L"\sigma_{q_{0}}",
    Latex.L"\sigma_{q_{1}}",
    Latex.L"\sigma_{q_{2}}",
    Latex.L"\sigma_{q_{3}}",
    Latex.L"\sigma_{q_{4}}",
    Latex.L"\sigma_{q_{5}}",
    Latex.L"\sigma_{q_{6}}",
    Latex.L"\sigma_{q_{7}}",
    Latex.L"\sigma_{q_{8}}",
    Latex.L"\sigma_{q_{9}}",
    Latex.L"\sigma_{q_{10}}",
    Latex.L"\sigma_{q_{11}}",
    Latex.L"\sigma_{q_{12}}",
    Latex.L"R_\mathrm{homogeneity}",
    Latex.L"R_\mathrm{uncoord}",
    Latex.L"\overline{s}",
    Latex.L"\sigma_{s}",
    Latex.L"\overline{r}_\mathrm{ring}",
    Latex.L"\sigma_{r_\mathrm{ring}}",
    Latex.L"\delta_\mathrm{c}",
    Latex.L"A_\mathrm{vertices  }",
    Latex.L"A_\mathrm{bonds}",
    Latex.L"\alpha",
]



# --- Plot 1: Explained variance per PC ---
pc_indices = collect(1:length(explained_variance_ratio))
plt1 = Plots.bar(
    pc_indices,
    explained_variance_ratio,
    xlabel = "PC",
    ylabel = "Explained variance ratio",
    legend = false,
)
#display(plt1)

# --- Plot 2: Cumulative explained variance ---
cumulative_variance = cumsum(explained_variance_ratio)
plt2 = Plots.plot(
    pc_indices,
    cumulative_variance,
    xlabel = "Number of PCs",
    ylabel = "Cumulative variance",
    label = "Cumulative variance",
    lw = 2,
    marker = :circle,
    legend = false
)
#display(plt2)

# --- Plot 3: Heatmap of component loadings ---
plt3 = Plots.heatmap(
    components,
    xlabel = "PC",
    colorbar_title = "Loading",
    c = :bluesreds,
    clim = (-0.6, 0.6),
    right_margin = 3Plots.mm,
    left_margin = 23Plots.mm,
    size = (500, 1300)
)

# put the defined labels on y axis
Plots.yticks!(plt3, (1:length(labels), labels))

#display(plt3)

# consider only the following indices of the labels in an additional plot p4
considered_indices = [1,2,3,39,33,34,42,4,41,37,38,5,6]

# --- Plot 4: Heatmap of component loadings (considered indices only) ---
plt4 = Plots.heatmap(
    reverse(components[considered_indices, :], dims=1),
    xlabel = "PC",
    colorbar_title = "Loading",
    c = :bluesreds,
    clim = (-0.6, 0.6),
    right_margin = 7Plots.mm,
    left_margin = 3Plots.mm,
    size = (500, 550)
)

# put the defined labels on y axis
Plots.yticks!(plt4, (1:length(considered_indices), reverse(labels[considered_indices])))


# Save Plots
Plots.savefig(plt1, plot_path*"explained_variance_per_pc.pdf")
Plots.savefig(plt2, plot_path*"cumulative_explained_variance.pdf")
Plots.savefig(plt3, plot_path*"pca_component_loadings_heatmap.pdf")
Plots.savefig(plt4, plot_path*"pca_component_loadings_heatmap_considered_indices.pdf")

