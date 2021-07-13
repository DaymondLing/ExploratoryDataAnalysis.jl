module ExploratoryDataAnalysis

# using Base: Real
using StatsBase
using DataFrames
using CategoricalArrays
using FreqTables
using NamedArrays

export infovalue
export cramerv
export mutualinfo
export eda

include("eda.jl")
include("utils.jl")

end
