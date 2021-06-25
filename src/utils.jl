"""
    ranks(x; groups = 10, rank = tiedrank, rev = false)

Return a variable which bins `x` into `groups` number of bins.
The `rank` keyword allows different ranking method;
use `rev = true` to reverse sort so that small bin number is large value of `x`.
Missing values are assigned to group `missing`.

Default values of `rank = tiedrank` and `rev = false` results
in similar grouping as SAS PROC RANK groups=n tied=mean.
"""
function ranks(x::Vector; groups::Int=10, rank=tiedrank, rev=false)
    nx = length(x)
    nx >= 2 || throw(ArgumentError("x length should be >= 2"))
    groups >= 2 || throw(ArgumentError("groups should be >= 2"))

    return floor.(Int32, rank(x; rev=rev) .* (groups / (nx + 1)))
end

function ranks(x::Vector{T} where {T>:Missing}; groups::Int=10, rank=tiedrank, rev=false)
    length(x) >= 2 || throw(ArgumentError("x length should be >= 2"))
    groups >= 2 || throw(ArgumentError("groups should be >= 2"))
    nx = count(!ismissing, x)

    return floor.(Union{Missing,Int32}, rank(x; rev=rev) .* (groups / (nx + 1)))
end
