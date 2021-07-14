module ExploratoryDataAnalysis

using StatsBase
using DataFrames
using CategoricalArrays
using FreqTables
using NamedArrays

export infovalue
export cramerv
export mutualinfo
export eda
export ranks

include("eda.jl")
include("ranks.jl")

end
