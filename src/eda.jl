"""
    norm1(v)

normalize `v` to sum 1
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

Computes the symmetric Kullback-Liebler Divergence between probabilities `p` and `q`
"""
function infovalue(p, q)
    pn = norm1(p)
    qn = norm1(q)

    return kldivergence(pn, qn) + kldivergence(qn, pn)
end

"""
    ϕ(f::Matrix{T}) where T <: Real

Computes the phi coefficient of a contingency table `f`, sqrt(χ² / n)
"""
function ϕ(f::Matrix{T}) where {T <: Real}
    minimum(size(f)) >= 2 || throw(ArgumentError("Matrix needed, not single row or column"))

    p = norm1(f)
    rp = sum(p; dims=2)                 # row marginal
    cp = sum(p; dims=1)                 # column marginal
    ep = rp .* cp                       # expected probabilities
    ϕ² = sum((p .- ep).^2 ./ ep)        # ϕ²

    return sqrt(ϕ²)
end

"""
    ϕ(x::Vector{T} where T<:Integer, y::Vector{T} where T<:Integer)

Computes the phi coefficient between two discrete vectors `x` and `y`
"""
ϕ(x::Vector{T} where {T <: Integer}, y::Vector{T} where {T <: Integer}) = ϕ(counts(x, y))
ϕ(x::AbstractVector, y::AbstractVector) = ϕ(freqtable(x, y).array)

"""
    mutualinfo(f::Matrix{T}) where T <: Real

Computes mutual information of frequency matrix `f`
"""
function mutualinfo(f::Matrix{T} where {T <: Real})
    minimum(size(f)) >= 2 || throw(ArgumentError("Matrix needed, not single row or column"))

    ps = norm1(f)
    px = sum(ps; dims=2)
    py = sum(ps; dims=1)

    return entropy(px) + entropy(py) - entropy(ps)
end

"""
    mutualinfo(x, y)

Computes mutual information between two discrete vectors `x` and `y`
"""
function mutualinfo(x::Vector{T} where {T <: Integer}, y::Vector{S} where {S <: Integer})
    return mutualinfo(counts(x, y))
end
mutualinfo(x::AbstractVector, y::AbstractVector) = mutualinfo(freqtable(x, y).array)

"""
    eda(df, target::Symbol)

Returns dataframe of measures of association

only real & string

- categoricalarray is good
-
"""
function eda(df::AbstractDataFrame, target::Symbol; groups=20)::AbstractDataFrame
    t = df[!, target]
    tnlvl = length(unique(t))
    tnlvl <= 1 && throw(ArgumentError("Target is single valued"))

    if tnlvl > nbins
        t = ranks(t, groups = groups)
        tnlvl = length(unique(t))
    end

    if tnlvl == 2
        out = DataFrame(MutualInfo=Float64[], Phi=Float64[], InfoValue=Float64[])
        for v in propertynames(df)
            v == target && continue
            if eltype(df[!, v]) <: Real
                vb =  ranks(df[!, v], groups = groups)
            else
                vb = df[!, v]
            end

            frq = freqtable(vb, t).array
            mutin = mutualinfo(frq)
            phi = ϕ(frq)
            iv = infovalue(frq[:, 1], frq[:, 2])

            push!(out, (mutin, phi, iv))
        end
    else
        out = DataFrame(MutualInfo=Float64[], Phi=Float64[])
        for v in propertynames(df)
            v == target && continue
            if eltype(df[!, v]) <: Real
                vb =  ranks(df[!, v], groups = groups)
            else
                vb = df[!, v]
            end

            frq = freqtable(vb, t).array
            mutin = mutualinfo(frq)
            phi = ϕ(frq)

            push!(out, (mutin, phi))
        end
    end

    out
end
