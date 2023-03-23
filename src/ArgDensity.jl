abstract type AbstractDensity <: Function end
abstract type AbstractGraphDensity <: AbstractDensity end

##############
# Coalescent #
##############

struct CoalDensity <: AbstractGraphDensity
    ## Number of leaves.
    n::Int

    ## Probability of a topology on n leaves.
    ptopo_log::BigFloat

    function CoalDensity(n)
        N = big(n)
        ptopo_log = (log ∘ factorial)(N - 1) - (log ∘ factorial)(2(N - 1))

        new(n, ptopo_log)
    end
end

function (D::CoalDensity)(arg::Arg; logscale = false)
    n = big(D.n)
    iter = enumerate((diff ∘ latitudes)(arg))
    init = latitude(arg, n + 1) * log(n) * log(n - one(n)) / BigFloat(2)

    plat_log = mapreduce(+, iter, init = init) do p
        k, Δ = big(first(p)), big(last(p))

        l = 2n - k
        λ = l * (l + one(l)) / BigFloat(2)

        log(λ) - λ * Δ
    end

    ret = D.ptopo_log + plat_log

    logscale ? ret : exp(ret)
end

############################
# Coalescent with mutation #
############################

struct CoalMutDensity <: AbstractGraphDensity
    dens_coal::CoalDensity

    ## Scaled mutation rate.
    μ::BigFloat

    ## Length of sequences.
    seq_length::BigFloat

    CoalMutDensity(n, μ, seq_length) = new(CoalDensity(n), μ, seq_length)
end

function (D::CoalMutDensity)(arg::Arg; logscale = false)
    dens_coal = D.dens_coal
    μ = D.μ
    l = D.seq_length

    m = (BigFloat ∘ nmutations)(arg)
    bl = (big ∘ branchlength_tree)(arg)

    pmut_log = m * log(μ) - big(0.5) * μ * l * bl - (log ∘ factorial)(m)
    ret = dens_coal(arg, logscale = true) + pmut_log

    logscale ? ret : exp(ret)
end

###################################
# Fréchet Phenotypes Distribution #
###################################

"""
    FrechetCoalDensity{T <: Real}

Density of a vector of phenotypes of type `T` conditional on an
Ancestral Recombination Graph.

# Fields
- `leaves_phenotypes::Vector{T}`: Entry `k` is the phenotype of leaf `k`.
- `α::Function`: Function mapping the positive real numbers to [0, 1].
- `pars::Dict{Symbol, Any}`: Parameters.
"""
struct FrechetCoalDensity{T} <: AbstractGraphDensity where T <: Real
    leaves_phenotypes::Vector{T}

    α::Function
    pars::Dict{Symbol, Any}

    FrechetCoalDensity(leaves_phenotypes::Vector{T};
                       α = t -> 1 - exp(-t),
                       pars = Dict{Symbol, Any}()) where T =
                       new{T}(leaves_phenotypes, α, pars)
end

q(ψ::Bool, x) = ψ + (ψ ? -1 : 1) * x

"""
    joint_frechet(arg, parent, phenotypes, α, p; logscale = false)
    marginal_frechet(phenotype, p)

Joint/marginal densities. The first one is the density of two
phenotypes associated with two coalescing sequences conditional on the
phenotype of their parent. The second one is the marginal density of a
phenotype.

`phenotypes` must be a `Vector{T}` of length 3 containing the phenotypes of the
parent, left child and right child in that order.
"""
function joint_frechet(arg, parent, children, phenotypes::AbstractVector{Bool},
                       α, p;
                       logscale = false)
    ϕp, Φc = first(phenotypes), last(phenotypes, 2)
    lat_parent = latitude(arg, parent)

    Δs = map(c -> lat_parent - latitude(arg, c), children)

    qϕ = Fix1(q, ϕp)
    probs = map(Δ -> qϕ(qϕ(p) * α(Δ)), Δs)

    ret_log = sum(x -> log(1.0 - q(first(x), last(x))),
                  zip(Φc, probs),
                  init = zero(Float64))

    logscale ? ret_log : exp(ret_log)
end

marginal_frechet(phenotype::Bool, p; logscale = false) =
    phenotype ? p : 1.0 - p

using Turing: @model, sample, MH, PG, SMC
using Distributions: logpdf
using Random: AbstractRNG, GLOBAL_RNG
function (D::FrechetCoalDensity)(rng::AbstractRNG, arg;
                                 logscale = false, M = 1000)
    p = D.pars[:p]::Float64
    α = D.α
    ni = nivertices(arg)

    ## Compute posterior density of the phenotypes of the coalescence vertices.
    @model function ivertices_pheno(arg, φ, p)
        n = nleaves(arg)

        ψ = Vector{Bool}(undef, ni)

        ψ[mrca(arg) - n] ~ Bernoulli(p)

        vstack = Int[]
        push!(vstack, mrca(arg))
        while !isempty(vstack)
            parent = pop!(vstack)
            q_parent = q(ψ[parent - n], p)
            t_parent = latitude(arg, parent)

            for child ∈ children(arg, parent)
                Δt = t_parent - latitude(arg, child)

                if isleaf(arg, child)
                    φ[child] ~ Bernoulli(q(ψ[parent - n], q_parent * α(Δt)))
                else
                    ψ[child - n] ~ Bernoulli(q(ψ[parent - n], q_parent * α(Δt)))
                    push!(vstack, child)
                end
            end
        end
    end

    m = ivertices_pheno(arg, D.leaves_phenotypes, p)
    sampler = PG(20)
    ps_sample = sample(rng, m, sampler, M, discard_initial = 100)
    ivaluations = convert(Matrix{Bool}, ps_sample.value[:, 1:end-1, 1].data)

    ret_log = zero(Float64)
    for ivaluation ∈ eachrow(ivaluations)
        valuation = vcat(D.leaves_phenotypes, ivaluation)
        ret_log += sum(ivertices(arg)) do ψ
            _children = children(arg, ψ)
            phenos = view(valuation, vcat(ψ, _children))
            joint_frechet(arg, ψ, _children, phenos, α, p, logscale = true)::Float64
        end
        ret_log += marginal_frechet(valuation[mrca(arg)], p,
                                    logscale = true)
    end

    ret_log /= M

    logscale ? ret_log : exp(ret_log)
end

(D::FrechetCoalDensity)(arg; logscale = false, M = 1000) =
    D(GLOBAL_RNG, arg, logscale = logscale, M = M)
