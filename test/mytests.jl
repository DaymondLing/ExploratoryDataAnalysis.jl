using ExploratoryDataAnalysis
using DataFrames
using CategoricalArrays
# using Plots

df = DataFrame(t = [fill(0, 500); fill(1, 500)], x1 = collect(1:1000), x2 = randn(1000))

df.tc = categorical(df.t)

eda(df, :t)

eda(df, :x1)

f = [9800 9925 11825; 190200 190075 188175]'
infovalue(f)
mutualinfo(f)
