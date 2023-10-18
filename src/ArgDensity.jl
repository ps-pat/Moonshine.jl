using StaticArrays: SA

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
        ptopo_log = log(N) + (N - 1) * log(2) - 2 * sum(log, 2:N)

        new(n, ptopo_log)
    end
end

function (D::CoalDensity)(arg::Arg; logscale = false)
    n = D.n

    ## Probability of the first coalescence.
    plat1_log = logccdf(Exponential(inv(n - 1)), first(latitudes(arg)))

    iter = enumerate((diff ∘ latitudes)(arg))
    plat_log = mapreduce(+, iter, init = zero(Float64)) do p
        k, Δ = first(p) + 1, last(p)

        logccdf(Exponential(inv(n - k)), Δ)
    end

    ret = D.ptopo_log + plat1_log + plat_log

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

    m = nmutations(arg)
    bl = branchlength_tree(arg)

    pmut_log = m * log(μ) - 0.5 * μ * l * bl - (log ∘ factorial)(m)
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
    leaves_phenotypes::Vector{Union{Missing,T}}

    α::Function
    pars::Dict{Symbol,Any}

    function FrechetCoalDensity(leaves_phenotypes::Vector{Union{Missing, S}};
                                α = (t, λ) -> 1 - exp(-t / λ),
                                pars = Dict{Symbol,Any}()) where {S}
        new{S}(leaves_phenotypes, α, pars)
    end
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
                      α, p)
    φp, φc = first(phenotypes), last(phenotypes)
    Δt = latitude(arg, parent) - latitude(arg, child)

    q_parent = Fix1(q, φp)
    prob = q_parent(q_parent(p) * α(Δt))

    BigFloat(q(!φc, prob))
end

dens_frechet(phenotype::Bool, p) = BigFloat(phenotype ? p : 1 - p)

## Root.
function cmatrix_frechet(arg, p)
    reshape([dens_frechet(φ, p) for φ in SA[false, true]], 2, 1)
end

## This method assumes that both vectors of phenotypes are sorted
## (false < true).
function cmatrix_frechet(arg, σ, φs_σ, δ, φs_δ, α, p)
    map(Iterators.product(φs_σ, φs_δ)) do (φ_σ, φ_δ)
        dens_frechet(arg, δ, σ, SA[φ_δ, φ_σ], α, p)
    end
end

function cmatrix_frechet(arg, phenotypes::AbstractVector{Union{Missing,Bool}},
                         σ, α, p)
    _parents = parents(arg, σ)
    δ = isempty(_parents) ? (zero ∘ eltype)(arg) : first(_parents)

    ## Root
    iszero(δ) && return cmatrix_frechet(arg, p)

    ## Assume that the phenotype of δ is unknown since it is an
    ## internal vertex.
    φs_δ = [false, true]

    ## The phenotype of σ might be known if it is a leaf.
    if isleaf(arg, σ)
        φ_σ = phenotypes[σ]
        φs_σ = ismissing(φ_σ) ? [false, true] : [φ_σ]
    else # non-root internal vertex
        φs_σ = [false, true]
    end

    cmatrix_frechet(arg, σ, φs_σ, δ, φs_δ, α, p)
end

function (D::FrechetCoalDensity{Bool})(arg, perm = 1:nleaves(arg))
    p = D.pars[:p]::Float64
    α = D.α
    leaves_phenotypes = D.leaves_phenotypes[perm]
    missing_phenotypes = findall(ismissing, leaves_phenotypes)
    ni = nivertices(arg)
    nmiss = length(missing_phenotypes)

    ## Only 1 vertex in ARG.
    iszero(ni) && return first(leaves_phenotype) ? p : 1 - p

    ## Belief propagation through postorder traversal.

    ## Stack used to compute the postorder traversal.
    vertices_stack = CheapStack(eltype(arg), nv(arg))

    ## TODO: update explanations.
    ## Stack used to store messages. The first `nmiss` dimensions are
    ## indexed by the missing phenotypes. The last dimension is
    ## indexed by the current phenotype.
    messages_stack = CheapStack(Matrix{BigFloat}, nv(arg))

    push!(vertices_stack, (minimum ∘ children)(arg, mrca(arg)))
    push!(vertices_stack, mrca(arg))
    v = (maximum ∘ children)(arg, mrca(arg))

    res = zero(BigFloat)
    while !isempty(vertices_stack)
        while !iszero(v)
            if !isleaf(arg, v)
                v_children = children(arg, v)
                push!(vertices_stack, minimum(v_children))
                push!(vertices_stack, v)

                v = maximum(v_children)
            else
                push!(vertices_stack, v)
                v = zero(v)
            end
        end

        v = pop!(vertices_stack)
        if !isleaf(arg, v) &&
           (first ∘ children)(arg, v) == first(vertices_stack)
            v2 = pop!(vertices_stack)
            push!(vertices_stack, v)
            v = v2
        else # Compute message!
            μ = cmatrix_frechet(arg, leaves_phenotypes, v, α, p)

            if !isleaf(arg, v)
                μ = (pop!(messages_stack) ⊙ pop!(messages_stack)) * μ
            end
            push!(messages_stack, μ)

            v = zero(v)
        end
    end

    pop!(messages_stack)
end
