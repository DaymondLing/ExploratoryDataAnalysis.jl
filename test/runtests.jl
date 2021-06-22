using ExploratoryDataAnalysis
using Test

using StatsBase
using FreqTables

a = [1, 2, 3, 1, 2, 3, 1, 1, 2, 2]
b = [11, 11, 12, 12, 11, 11, 12, 12, 11, 12]

norm(v::AbstractVector) = v ./ sum(v)

@testset "phi" begin
    @test round(ϕ(a, b), digits=7) == 0.4472136
    @test round(ϕ(counts(a, b)), digits=7) == 0.4472136
    @test round(ϕ(freqtable(a, b)), digits=7) == 0.4472136
    @test round(ϕ(freqtable(a, b).array), digits=7) == 0.4472136
end

@testset "infovalue" begin
    @test round(ϕ(a, b), digits=7) == 0.4472136
    @test round(ϕ(counts(a, b)), digits=7) == 0.4472136
    @test round(ϕ(freqtable(a, b)), digits=7) == 0.4472136
    @test round(ϕ(freqtable(a, b).array), digits=7) == 0.4472136
end
