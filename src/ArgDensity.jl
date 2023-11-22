using StaticArrays: SA

using Distributions: logpdf

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

function (D::CoalDensity)(tree::Tree; logscale = false)
    n = D.n

    ## Probability of the first coalescence.
    plat1_log = logpdf(Exponential(inv(n - 1)), first(latitudes(tree)))

    iter = enumerate((diff ∘ latitudes)(tree))
    plat_log = mapreduce(+, iter, init = zero(Float64)) do p
        k, Δ = first(p) + 1, last(p)

        logpdf(Exponential(inv(n - k)), Δ)
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

function (D::CoalMutDensity)(tree::Tree; logscale = false)
    dens_coal = D.dens_coal
    μ = D.μ
    l = D.seq_length

    m = nmutations(tree)
    bl = branchlength(tree)

    pmut_log = m * log(μ) - 0.5 * μ * l * bl - sum(log, 2:m)
    ret = dens_coal(tree, logscale = true) + pmut_log

    logscale ? ret : exp(ret)
end

############################
# Phenotypes Distributions #
############################

export PhenotypeDensity
"""
    PhenotypeDensity

Density of a vector of phenotypes of type `T` conditional on a genealogy.

# Fields

  - `Φ::Vector{Union{Missing, T}}`: entry `k` is the phenotype of leaf `k`
  - `copula::C`: `<:AbstractΦCopula`
"""
struct PhenotypeDensity{T,C<:AbstractΦCopula} <: AbstractGraphDensity
    Φ::Vector{Union{Missing,T}}
    copula::C
end

function (D::PhenotypeDensity)(tree::Tree)
    ## TODO: Manage this case more cleanly
    nleaves(tree) <= 1 && error("tree must have at least two leaves")

    ## Belief propagation through postorder traversal.
    copula = D.copula
    Φ = D.Φ

    ## Stack used to compute the postorder traversal.
    vertices_stack = CheapStack(eltype(tree), nv(tree))

    push!(vertices_stack, (minimum ∘ children)(tree, mrca(tree)))
    push!(vertices_stack, mrca(tree))
    v = (maximum ∘ children)(tree, mrca(tree))

    ## Stack used to store messages. Entry (i, j) contains the value of
    ## f(φ = g(i) | ψ = g(j)) where φ and ψ are the source and destination of
    ## the message respectively and g is some functions that compute phenotype
    ## from an integer.
    messages_stack = CheapStack(Matrix{Float64}, nv(tree))

    res = zero(Float64)
    while !isempty(vertices_stack)
        while !iszero(v)
            if !isleaf(tree, v)
                v_children = children(tree, v)
                push!(vertices_stack, minimum(v_children))
                push!(vertices_stack, v)

                v = maximum(v_children)
            else
                push!(vertices_stack, v)
                v = zero(v)
            end
        end

        v = pop!(vertices_stack)
        if !isleaf(tree, v) &&
           (first ∘ children)(tree, v) == first(vertices_stack)
            v2 = pop!(vertices_stack)
            push!(vertices_stack, v)
            v = v2
        else # Compute message!
            μ = cmatrix(tree, copula, v, Φ)

            if !isleaf(tree, v)
                μ = (pop!(messages_stack) ⊙ pop!(messages_stack)) * μ
            end
            push!(messages_stack, μ)

            v = zero(v)
        end
    end

    pop!(messages_stack)
end

"""
    cmatrix(tree, copula, σ, phenotypes)

Matrix containing relevant points of the conditional distribution of a vector
of phenotypes. Used by the belief propagation algorithm.
"""
function cmatrix(tree::Tree, copula::AbstractΦCopula{PhenotypeBinary}, σ, phenotypes)
    isroot(tree, σ) &&
        return reshape([pdf(copula, φ) for φ in (false, true)], 2, 1)

    δ = dad(tree, σ)

    ## Assume that the phenotype of δ is unknown since it is an
    ## internal vertex.
    φsδ = (false, true)

    ## The phenotype of σ might be known if it is a leaf.
    if isleaf(tree, σ) && !ismissing(phenotypes[σ])
        φsσ = [phenotypes[σ]]
    else # non-root internal vertex
        φsσ = (false, true)
    end

    cmatrix(tree, copula, σ, φsσ, δ, φsδ)
end

function cmatrix(tree::Tree, copula::AbstractΦCopula{PhenotypeBinary}, σ, φsσ, δ, φsδ)
    map(Iterators.product(φsσ, φsδ)) do (φσ, φδ)
        conditional_pdf(copula, φσ, φδ, distance(tree, σ, δ))
    end
end
