abstract type AbstractDensity <: Function end
abstract type AbstractGraphDensity <: AbstractDensity end

##############
# Coalescent #
##############

export CoalDensity
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

export CoalMutDensity
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
struct FrechetCoalDensity{T} <: AbstractGraphDensity
    leaves_phenotypes::Vector{Union{Missing, T}}

    α::Function
    pars::Dict{Symbol, Any}

    FrechetCoalDensity(leaves_phenotypes::Vector{Union{Missing, S}};
                       α = t -> 1 - exp(-t),
                       pars = Dict{Symbol, Any}()) where S =
                       new{S}(leaves_phenotypes, α, pars)
end
export FrechetCoalDensity

q(ψ::Bool, x) = ψ + (ψ ? -1 : 1) * x

"""
    joint_frechet(arg, parent, children, phenotypes, α, p; logscale = false)
    joint_frechet(arg, parent, child, phenotypes, α, p; logscale = false)
    marginal_frechet(phenotype, p)

Joint/marginal densities. The first one is the density of two
phenotypes associated with two coalescing sequences conditional on the
phenotype of their parent. The second one is the marginal density of a
phenotype.

`phenotypes` must be a `Vector{T}` of length 3 containing the phenotypes of the
parent, left child and right child in that order.
"""
function dens_frechet end

function dens_frechet(arg, parent, child::AbstractVector,
                      phenotypes::AbstractVector{Bool},
                      α, p;
                      logscale = false)
    ϕp, Φc = first(phenotypes), last(phenotypes, 2)
    lat_parent = latitude(arg, parent)

    Δs = map(c -> lat_parent - latitude(arg, c), child)

    qϕ = Fix1(q, ϕp)
    probs = map(Δ -> qϕ(qϕ(p) * α(Δ)), Δs)

    ret_log = sum(x -> (log ∘ q)(!first(x), last(x)),
                  zip(Φc, probs),
                  init = zero(Float64))

    logscale ? ret_log : exp(ret_log)
end

## For a parent-child pair (eq. 6 of the paper).
function dens_frechet(arg, parent, child, phenotypes::AbstractVector{Bool},
                      α, p;
                      logscale = false)
    φp, φc = first(phenotypes), last(phenotypes)
    Δt = latitude(arg, parent) - latitude(arg, child)

    q_parent = Fix1(q, φp)
    prob = q_parent(q_parent(p) * α(Δt))

    ret_log = BigFloat((log ∘ q)(!φc, prob))

    logscale ? ret_log : exp(ret_log)
end

function dens_frechet(phenotype::Bool, p; logscale = false)
    ret_log = BigFloat(phenotype ? p : 1 - p)
    logscale ? ret_log : exp(ret_log)
end

function compute_integrand(arg, children, valuation, α, p; logscale = false)
    ret_log = sum(children) do child
        parent = (first ∘ parents)(arg, child)
        φs = convert(Vector{Bool}, view(valuation, [parent, child]))
        dens_frechet(arg, parent, child, φs, α, p, logscale = true)
    end

    logscale ? ret_log : exp(ret_log)
end

using Turing: @model, sample, MH, PG, SMC
using Distributions: logpdf
using Random: AbstractRNG, GLOBAL_RNG
using Statistics: mean
function (D::FrechetCoalDensity{Bool})(rng::AbstractRNG, arg; M = 1000)
    p = D.pars[:p]::Float64
    α = D.α
    leaves_phenotypes = D.leaves_phenotypes
    missing_phenotypes = findall(ismissing,leaves_phenotypes)
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
    ivaluations = convert(Matrix{eltype(leaves_phenotypes)},
                          ps_sample.value[:, 1:end-2, 1].data)

    int_partials_log::Vector{BigFloat} = map(eachrow(ivaluations)) do ivaluation
        valuation = vcat(leaves_phenotypes, ivaluation)
        vs = setdiff(vertices(arg), vcat(mrca(arg), missing_phenotypes))

        compute_integrand(arg, vs, valuation, α, p, logscale = true) +
            dens_frechet(Bool(valuation[mrca(arg)]), p, logscale = true)
    end

    function(phenotypes; logscale = false)
        φ_leaves = deepcopy(leaves_phenotypes)
        φ_leaves[missing_phenotypes] .= phenotypes

        ret_log = mean(zip(int_partials_log, eachrow(ivaluations))) do x
            int_partial_log, ivaluation = x

            valuation = vcat(φ_leaves, ivaluation)
            int_partial_log +
                compute_integrand(arg, missing_phenotypes, valuation, α, p,
                                  logscale = true)
        end

        logscale ? ret_log : exp(ret_log)
    end
end

(D::FrechetCoalDensity)(arg; M = 1000) =
    D(GLOBAL_RNG, arg, M = M)
