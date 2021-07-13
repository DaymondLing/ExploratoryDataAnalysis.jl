"""
    norm1(v)

Normalize `v` to sum 1.
"""
function norm1(v)
    vs = sum(v)
    if vs ≈ 1.0
        return v
    else
        return v ./ vs
    end
end

"""
    infovalue(p, q)

Compute the symmetric relative entropy or Kullback-Liebler Divergence between `p` and `q`.\\
Value is bounded within [0, ∞].
"""
function infovalue(p, q)
    pn = norm1(p)
    qn = norm1(q)

    kldivergence(pn, qn) + kldivergence(qn, pn)
end

"""
    ϕ(f::Matrix{T}) where T <: Real

Compute the phi coefficient of a contingency table `f`, sqrt(χ² / n).\\
Value is bounded within [0, √min(row-1,col-1)].
"""
function ϕ(f::Matrix{T}) where {T <: Real}
    minimum(size(f)) >= 2 || throw(ArgumentError("Matrix needed, not single row or column"))

    p = norm1(f)
    rp = sum(p; dims=2)                 # row marginal
    cp = sum(p; dims=1)                 # column marginal
    ep = rp .* cp                       # expected probabilities
    ϕ² = sum((p .- ep).^2 ./ ep)        # ϕ²

    sqrt(ϕ²)
end
ϕ(f::NamedArray) = ϕ(f.array)

"""
    ϕ(x::Vector{T} where T<:Integer, y::Vector{T} where T<:Integer)

Compute the phi coefficient of two discrete vectors `x` and `y`.\\
Value is bounded within [0, √min(row-1,col-1)].
"""
ϕ(x::Vector{T} where {T <: Integer}, y::Vector{T} where {T <: Integer}) = ϕ(counts(x, y))
ϕ(x::AbstractVector, y::AbstractVector) = ϕ(freqtable(x, y).array)

"""
    mutualinfo(f::Matrix{T}) where T <: Real

Compute the mutual information of frequency matrix `f`.\\
Value is bounded within [0, ∞].
"""
function mutualinfo(f::Matrix{T} where {T <: Real})
    minimum(size(f)) >= 2 || throw(ArgumentError("Matrix needed, not single row or column"))

    ps = norm1(f)
    px = sum(ps; dims=2)
    py = sum(ps; dims=1)

    entropy(px) + entropy(py) - entropy(ps)
end
mutualinfo(f::NamedArray) = mutualinfo(f.array)

"""
    mutualinfo(x::Vector, y::Vector)

Compute the mutual information of two discrete vectors `x` and `y`.\\
Value is bounded within [0, ∞].
"""
function mutualinfo(x::Vector{T} where {T <: Integer}, y::Vector{S} where {S <: Integer})
    mutualinfo(counts(x, y))
end
mutualinfo(x::AbstractVector, y::AbstractVector) = mutualinfo(freqtable(x, y).array)

"""
    eda(df::AbstractDataFrame, target::Symbol)

Return a dataframe of Mutual Information and ϕ coefficient between `target` and other
variables in `df`. If `target` is binary, Information Value is also returned.

The dataframe is sorted by descending Mutual Information.
"""
function eda(df::AbstractDataFrame, target::Symbol; groups=20)::AbstractDataFrame
    t = df[!, target]
    tnlvl = length(unique(t))
    tnlvl <= 1 && throw(ArgumentError("Target is single valued"))

    if tnlvl > groups
        t = ranks(t, groups = groups)
        tnlvl = length(unique(t))
    end

    println("Target: $target   type: $(typeof(t))   Levels: $tnlvl")

    if tnlvl == 2
        out = DataFrame(Variable=Symbol[], Vartype=DataType[], Varlvls=Int[],
                MutualInfo=Float64[], Phi=Float64[], InfoValue=Float64[])
        for v in propertynames(df)
            v == target && continue

            vtype = eltype(df[!, v])
            if vtype <: Union{Missing, Real}
                vb =  ranks(df[!, v], groups = groups)
            else
                vb = df[!, v]
            end

            frq = freqtable(vb, t).array
            vnlvl = size(frq, 1)
            if vnlvl == 1
                println("Warning: [$v] is singled valued, skipped.")
                continue
            elseif vnlvl >= 100
                println("Warning: [$v] has more than 100 levels, suspicious")
            end

            mutin = mutualinfo(frq)
            phi = ϕ(frq)
            iv = infovalue(frq[:, 1], frq[:, 2])

            push!(out, (v, eltype(df[!, v]), vnlvl, mutin, phi, iv))
        end
    else
        out = DataFrame(Variable=Symbol[], Vartype=DataType[], Varlvls=Int[],
                MutualInfo=Float64[], Phi=Float64[])
        for v in propertynames(df)
            v == target && continue

            vtype = eltype(df[!, v])
            if vtype <: Real
                vb =  ranks(df[!, v], groups = groups)
            else
                vb = df[!, v]
            end

            frq = freqtable(vb, t).array
            vnlvl = size(frq, 1)
            if vnlvl == 1
                println("Warning: [$v] is singled valued, skipped.")
                continue
            elseif vnlvl >= 100
                println("Warning: [$v] has more than 100 levels, suspicious")
            end

            mutin = mutualinfo(frq)
            phi = ϕ(frq)

            push!(out, (v, eltype(df[!, v]), vnlvl, mutin, phi))
        end
    end

    sort!(out, [:MutualInfo, :Phi], rev = true)
end
