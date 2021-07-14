using ExploratoryDataAnalysis
using Test

using StatsBase
using FreqTables

a = [1, 2, 3, 1, 2, 3, 1, 1, 2, 2]
b = [11, 11, 12, 12, 11, 11, 12, 12, 11, 12]

@testset "Cramer's V" begin
    @test round(cramerv(a, b), digits = 7) == 0.4472136
    @test round(cramerv(counts(a, b)), digits = 7) == 0.4472136
    @test round(cramerv(freqtable(a, b)), digits = 7) == 0.4472136
    @test round(cramerv(freqtable(a, b).array), digits = 7) == 0.4472136
end

@testset "infovalue" begin
    @test round(infovalue(counts(a, b)), digits = 7) == 0.8788898
    @test round(infovalue(freqtable(a, b)), digits = 7) == 0.8788898
    @test round(infovalue(freqtable(a, b).array), digits = 7) == 0.8788898
end

@testset "mutual info" begin
    @test round(mutualinfo(a, b), digits = 7) == 0.1046496
    @test round(mutualinfo(counts(a, b)), digits = 7) == 0.1046496
    @test round(mutualinfo(freqtable(a, b)), digits = 7) == 0.1046496
    @test round(mutualinfo(freqtable(a, b).array), digits = 7) == 0.1046496
end
