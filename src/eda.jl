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

Compute the symmetric Relative Entropy (Kullback-Liebler Divergence) between `p` and `q`,
two probability vectors or two frequency count vectors.
Value is bounded within [0, ∞].
"""
function infovalue(p, q)
    pn = norm1(p)
    qn = norm1(q)

    kldivergence(pn, qn) + kldivergence(qn, pn)
end
infovalue(f::Matrix) = infovalue(f[:, 1], f[:, 2])
infovalue(f::NamedArray) = infovalue(f.array)


"""
    cramerv(f::Matrix{T}) where T <: Real

Compute Cramer's V of a contingency table `f`.\\
ϕ coefficient is √(χ²/n), bounded within [0, √min(row-1,col-1)],
Cramer's V is ϕ/√min(row-1,col-1), bounded within [0,1].
"""
function cramerv(f::Matrix{T}) where {T<:Real}
    minimum(size(f)) >= 2 || throw(ArgumentError("Matrix needed, not single row or column"))

    p = norm1(f)
    rp = sum(p; dims = 2)                 # row marginal
    cp = sum(p; dims = 1)                 # column marginal
    ep = rp .* cp                       # expected probabilities
    ϕ² = sum((p .- ep) .^ 2 ./ ep)        # ϕ²

    sqrt(ϕ² / (minimum(size(f)) - 1))   # Cramer's V
end
cramerv(f::NamedArray) = cramerv(f.array)


"""
    cramerv(x::Vector{T} where T<:Integer, y::Vector{T} where T<:Integer)

Compute Cramer's V of two presumed discrete vectors `x` and `y`.\\
"""
cramerv(x::Vector{T} where {T<:Integer}, y::Vector{T} where {T<:Integer}) =
    cramerv(freqtable(x, y).array)


"""
    cramerv(x::Vector, y::Vector; groups = 10)

Compute Cramer's V of two presumed continuous vectors `x` and `y`,
vectors are binned into `groups` first.
"""
function cramerv(x::AbstractVector, y::AbstractVector; groups = 10)
    xr = ranks(x; groups = groups)
    yr = ranks(y, groups = groups)

    cramerv(freqtable(xr, yr).array)
end


"""
    mutualinfo(f::Matrix{T}) where T <: Real

Compute the mutual information of frequency matrix `f`; value is bounded within [0, ∞].
"""
function mutualinfo(f::Matrix{T} where {T<:Real})
    minimum(size(f)) >= 2 || throw(ArgumentError("Matrix needed, not single row or column"))

    ps = norm1(f)
    px = sum(ps; dims = 2)
    py = sum(ps; dims = 1)

    entropy(px) + entropy(py) - entropy(ps)
end
mutualinfo(f::NamedArray) = mutualinfo(f.array)


"""
    mutualinfo(x::Vector{T} where T<:Integer, y::Vector{T} where T<:Integer)

Compute mutual information of two presumed discrete vectors `x` and `y`.\\
"""
mutualinfo(x::Vector{T} where {T<:Integer}, y::Vector{T} where {T<:Integer}) =
    mutualinfo(freqtable(x, y).array)


"""
    mutualinfo(x::Vector, y::Vector)

Compute the mutual information of two presumed continuous vectors `x` and `y`.
"""
function mutualinfo(x::AbstractVector, y::AbstractVector; groups = 10)
    xr = ranks(x; groups = groups)
    yr = ranks(y, groups = groups)

    mutualinfo(freqtable(xr, yr).array)
end


"""
    eda(df::AbstractDataFrame, target::Symbol)

Return a dataframe of Mutual Information and Cramer's V between `target` and other
variables in `df`. If `target` is binary, Information Value is also returned.

The dataframe is sorted by descending Mutual Information.
"""
function eda(df::AbstractDataFrame, target::Symbol; groups = 20)::AbstractDataFrame
    t = df[!, target]
    tnlvl = length(unique(t))
    tnlvl <= 1 && throw(ArgumentError("Target is single valued"))

    if tnlvl > groups
        t = ranks(t, groups = groups)
        tnlvl = length(unique(t))
    end

    println("Target: $target   type: $(typeof(t))   Levels: $tnlvl")

    lk = ReentrantLock()

    if tnlvl == 2
        out = DataFrame(
            Variable   = Symbol[],
            Vartype    = String[],
            Varlvls    = Int[],
            MutualInfo = Float64[],
            CramerV    = Float64[],
            InfoValue  = Float64[],
        )

        Threads.@threads for v in propertynames(df)
            v == target && continue

            vtype = eltype(df[!, v])
            if vtype <: Union{Missing,Real}
                vb = ranks(df[!, v], groups = groups)
            else
                vb = df[!, v]
            end

            frq = freqtable(vb, t).array
            vnlvl = size(frq, 1)
            if vnlvl == 1
                println("Warning: [$v] is singled valued, skipped.")
                continue
            elseif vnlvl > 50
                println("Warning: [$v] has more than 100 levels, consider binning")
            end

            mutin = mutualinfo(frq)
            cramv = cramerv(frq)
            iv = infovalue(frq[:, 1], frq[:, 2])

            lock(lk)
            push!(out, (v, string(vtype), vnlvl, mutin, cramv, iv))
            unlock(lk)
        end
    else
        out = DataFrame(
            Variable   = Symbol[],
            Vartype    = String[],
            Varlvls    = Int[],
            MutualInfo = Float64[],
            CramerV    = Float64[],
        )

        Threads.@threads for v in propertynames(df)
            v == target && continue

            vtype = eltype(df[!, v])
            if vtype <: Union{Missing,Real}
                vb = ranks(df[!, v], groups = groups)
            else
                vb = df[!, v]
            end

            frq = freqtable(vb, t).array
            vnlvl = size(frq, 1)
            if vnlvl == 1
                println("Warning: [$v] is singled valued, skipped.")
                continue
            elseif vnlvl > 50
                println("Warning: [$v] has more than 100 levels, consider binning")
            end

            mutin = mutualinfo(frq)
            cramv = cramerv(frq)

            lock(lk)
            push!(out, (v, string(vtype), vnlvl, mutin, cramv))
            unlock(lk)
        end
    end

    sort!(out, [:MutualInfo, :CramerV], rev = true)
end
