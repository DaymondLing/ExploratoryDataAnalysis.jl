module ExploratoryDataAnalysis

using StatsBase
using DataFrames
using CategoricalArrays
using FreqTables

export infovalue
export ϕ
export mutualinfo
export eda

include("eda.jl")
include("utils.jl")

end
