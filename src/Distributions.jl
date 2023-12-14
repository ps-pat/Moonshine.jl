##########################
# Multivariate Bernoulli #
##########################

using Distributions

import Base: length, eltype

import Distributions: insupport, pdf,
                      params, partype,
                      mean, var, cov, entropy,
                      _logpdf,
                      _rand!

using Combinatorics

using LinearAlgebra

export BernoulliMulti
"""
    BernoulliMulti

Distribution of a 0-1 vector.

# Fields

  - `p::Vector{T}`: probability distribution of the vector. The probability of
    configuration c = b1 b2 ... bn is the binary representation of the index of
    its probability in `p`.

# Limitations

Due to implementation, the dimensionality of the distribution cannot exceed the
word size of the machine, i.e. `Sys.WORD_SIZE`.
"""
struct BernoulliMulti{T<:Real} <: DiscreteMultivariateDistribution
    p::Vector{T}
end

@generated eltype(::Type{BernoulliMulti}) = eltype(Bernoulli)
@generated eltype(::BernoulliMulti) = eltype(BernoulliMulti)

length(d::BernoulliMulti) = round(Int, log2(length(d.p)))

insupport(d::BernoulliMulti, x::AbstractVector) = all(z -> iszero(z) || isone(z), x)

params(d::Bernoulli) = (d.p,)
@generated partype(::Bernoulli{T}) where T = T

function pdf(d::BernoulliMulti, x::BitVector)
    getindex(d.p, bitrotate(bitreverse(first(x.chunks)), length(x)) + 1)
end

function pdf(d::BernoulliMulti, x::AbstractVector{Bool})
    idx = zero(Int)
    n = length(x)

    for (k, b) ∈ enumerate(x)
        idx |= b << (n - k)
    end

    getindex(d.p, idx + 1)
end

pdf(d::BernoulliMulti, x::AbstractVector{<:Real}) = pdf(d, convert(BitVector, x))

export marginal
"""
    marginal(d::MultivariateDistribution, idx)

Compute the marginal distribution of random variables from their joint
distribution.

# Type stability

Type stability of methods is not guaranteed.
`marginal(::BernoulliMulti, idx)` will return either a `BernoulliMulti` or
a `Bernoulli` depending on the length of `idx` for instance.
"""
function marginal(d::BernoulliMulti{T}, idx) where T
    (isone ∘ length)(idx) && return Bernoulli{T}(prob1(d, idx))

    sort!(idx, rev = true)

    n = length(d)
    idx_integrate = combinations(setdiff(1:n, idx))
    newp = Vector{T}(undef, 2^length(idx))

    for (k_newp, is) ∈ enumerate(combinations(idx))
        k_newp += 1

        k_dist = mapreduce(x -> 2^(n - x), ⊻, is)
        newp[k_newp] = d.p[k_dist + 1]

        for js ∈ idx_integrate
            kk = mapreduce(x -> 2^(n - x), ⊻, js)
            newp[k_newp] += d.p[k_dist ⊻ kk + 1]
        end
    end

    newp[1] = 1 - sum(view(newp, range(2, length(newp))))

    BernoulliMulti{T}(newp)
end

_logpdf(d::BernoulliMulti, x) = (log ∘ pdf)(d, x)

function mean(d::BernoulliMulti{T}) where T
    n = length(d)
    μ = zeros(T, n)

    range_step = 2^n
    range_length = one(Int)
    for k ∈ eachindex(μ)
        r1 = range(0, length = range_length, step = range_step)

        range_step ÷= 2
        range_length *= 2

        r2 = range(range_step + 1, length = range_step)
        it = Iterators.map(splat(+), Iterators.product(r1, r2))

        for i ∈ it
            μ[k] += d.p[i]
        end
    end

    μ
end

function var(d::BernoulliMulti)
    μ = mean(d)
    @. μ * (1 - μ)
end

function prob1(d, idx)
    n = length(d)
    k = mapreduce(x -> 2^(n - x), ⊻, idx)

    acc = d.p[k + 1]
    for is ∈ (combinations ∘ setdiff)(1:n, idx)
        kk = mapreduce(x -> 2^(n - x), ⊻, is)
        acc += d.p[k ⊻ kk + 1]
    end

    acc
end

function cov(d::BernoulliMulti{T}) where T
    n = length(d)
    μ = mean(d)

    M = Matrix{T}(undef, n, n)
    for i ∈ 1:n
        M[i, i] = prob1(d, i)
        for j ∈ range(i + 1, n)
            M[j, i] = prob1(d, (j, i))
        end
    end

    Symmetric(M, :L) - μ * μ'
end

entropy(d::BernoulliMulti{T}) where T = sum(p -> -p * log(p), d.p, init = zero(T))

function _rand!(rng::AbstractRNG,
                d::BernoulliMulti,
                x::AbstractArray{<:Real})
    n = length(d)
    m = length(x) ÷ n

    _cumsum = cumsum(d.p)
    for i ∈ 1:m
        xdi = bitrotate(bitreverse(findfirst(>(rand(rng)), _cumsum) - 1), n)
        for j ∈ 1:n
            x[j, i] = xdi & 1
            xdi >>= 1
        end
    end

    x
end
