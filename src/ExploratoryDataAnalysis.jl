module ExploratoryDataAnalysis

# using Base: Real
using StatsBase
using DataFrames
using CategoricalArrays
using FreqTables
using NamedArrays

export infovalue
export ϕ
export mutualinfo
export eda

include("eda.jl")
include("utils.jl")

end
