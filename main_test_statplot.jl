using StatsPlots # no need for `using Plots` as that is reexported here
gr(size=(400,300))

#%%

using DataFrames
#using IndexedTables
df = DataFrame(a = 1:10, b = 10 .* rand(10), c = 10 .* rand(10))
@df df plot(:a, [:b :c], colour = [:red :blue])
@df df scatter(:a, :b, markersize = 4 .* log.(:c .+ 0.1))
t = table(1:10, rand(10), names = [:a, :b]) # IndexedTable
@df t scatter(2 .* :b)

#%%

@df df plot(:a, cols(2:3), colour = [:red :blue])

#%%

s = :b
@df df plot(:a, cols(s))

#%%

df[:red] = rand(10)
@df df plot(:a, [:b :c], colour = ^([:red :blue]))

#%%