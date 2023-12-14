using StaticArrays: SA

using Distributions: logpdf

#############
# Interface #
#############

export genpars
"""
    genpars(D)

Returns the genetical parameters associated with a density in the form of a
named tuple.

# Names
Names of the parameters are standardized as follow:
- `seq_length`: length of the haplotypes
- `Ne`: effective population size
- `μ_loc`: per locus mutation rate
- `ρ_loc`: per locus recombination rate
- `positions`: positions of the markers normalized in [0, 1]
"""
function genpars end

##############
# Coalescent #
##############

export CoalDensity
struct CoalDensity
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

genpars(D::CoalDensity) = (Ne = zero(Float64),
                           μ_loc = zero(Float64),
                           ρ_loc = zero(Float64),
                           seq_length = zero(Float64),
                           positions = [0])

############################
# Coalescent with mutation #
############################

export CoalMutDensity
"""
    CoalMutDensity

Density of a genealogy with mutations according to the coalescent.
Mutations are assumed to be selection neutral i.e. independent of the
coalescent process.

See also [`CoalDensity`](@ref).

# Fields
- `fC::CoalDensity`: density of the genealogy (excluding mutations)
- `Ne::BigFloat`: effective population size
- `μ_loc::BigFloat`: per locus mutation rate
- `seq_length`: length of the sequences
"""
struct CoalMutDensity
    fC::CoalDensity

    Ne::BigFloat

    μ_loc::BigFloat

    seq_length::BigFloat

    CoalMutDensity(n, Ne, μ, seq_length) = new(CoalDensity(n), Ne, μ, seq_length)
end

genpars(D::CoalMutDensity) = (Ne = D.Ne,
                              μ_loc = D.μ_loc,
                              ρ_loc = zero(Float64),
                              seq_length = D.seq_length,
                              positions = [0])

function (D::CoalMutDensity)(tree::Tree; logscale = false)
    fC = D.fC
    μ = D.μ_loc * D.Ne
    l = D.seq_length

    m = nmutations(tree)
    bl = branchlength(tree)

    pmut_log = m * log(μ) - 0.5 * μ * l * bl - sum(log, 2:m, init = zero(BigFloat))
    ret = fC(tree, logscale = true) + pmut_log

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
struct PhenotypeDensity{T,C<:AbstractΦCopula}
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

    ## Stack used to store messages. Entry (i, j) contains the value
    ## of f(φ = g(i) | ψ = g(j), t) where φ and ψ are the source and
    ## destination of the message respectively and g is some functions
    ## that compute phenotype from an integer.
    messages_stack = CheapStack(Matrix{Float64}, nv(tree))

    ## Conditional and marginal pdfs for the phenotype.
    marginal_pdf = pdf_marginal(copula)
    conditional_pdf = pdf_conditional(copula)

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
            !isempty(vertices_stack) &&
            (first ∘ children)(tree, v) == first(vertices_stack)
            v2 = pop!(vertices_stack)
            push!(vertices_stack, v)
            v = v2
        else # Compute message!
            pdf = isempty(vertices_stack) ? marginal_pdf : conditional_pdf
            μ = cmatrix(tree, pdf, v, Φ)

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
    cmatrix(tree, pdf, σ, Φ)
    cmatrix(tree, pdf, σ, φsσ, δ, φsδ)

Matrix containing relevant points of the conditional distribution of a vector
of phenotypes. Used by the belief propagation algorithm.
"""
function cmatrix(tree::Tree, pdf, σ, Φ)
    isroot(tree, σ) &&
        return reshape([pdf(φ) for φ in (false, true)], 2, 1)

    δ = dad(tree, σ)

    ## Assume that the phenotype of δ is unknown since it is an
    ## internal vertex.
    φsδ = (false, true)

    ## The phenotype of σ might be known if it is a leaf.
    if isleaf(tree, σ)
        φσ = Φ[σ]
        φsσ = ismissing(φσ) ? [false, true] : [φσ]
    else # non-root internal vertex
        φsσ = [false, true]
    end

    cmatrix(tree, pdf, σ, φsσ, δ, φsδ)
end

function cmatrix(tree::Tree, pdf, σ, φsσ, δ, φsδ)
    ## Julia issue #46331.
    ret = Matrix{Float64}(undef, length(φsσ), length(φsδ))

    map!(ret, collect(Iterators.product(φsσ, φsδ))) do (φσ, φδ)
        Δlat = latitude(tree, δ) - latitude(tree, σ)
        pdf(φσ, φδ, Δlat)
    end

    ret
end
