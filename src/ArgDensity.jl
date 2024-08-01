using Distributions: logpdf

using DataStructures: Stack, Queue, enqueue!, dequeue!

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
- `μloc`: per locus mutation rate
- `ρloc`: per locus recombination rate
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

    ## Positions of the markers on [0,1].
    positions::Vector{Float64}

    function CoalDensity(n, positions = [0])
        N = big(n)
        ptopo_log = log(N) + (N - 1) * logtwo - 2 * sum(log, 2:N)

        new(n, ptopo_log, positions)
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
                           μloc = zero(Float64),
                           seq_length = zero(Float64),
                           positions = positions)

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

    μloc::BigFloat

    seq_length::BigFloat

    CoalMutDensity(n, Ne, μ, seq_length, positions = [0]) =
        new(CoalDensity(n, positions), Ne, μ, seq_length)
end

genpars(D::CoalMutDensity) = (Ne = D.Ne,
                              μloc = D.μloc,
                              seq_length = D.seq_length,
                              positions = D.fC.positions)

function (D::CoalMutDensity)(tree::Tree; logscale = false)
    fC = D.fC
    μ = D.μloc * D.Ne
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

function postorder_traversal(tree::Tree)
    ## Stack used for postorder traversal.
    vstack = Stack{eltype(tree)}(min(nleaves(tree) - 1, 1024))
    push!(vstack, mrca(tree))

    ## Queue storing vertices as visited following postorder traversal.
    vqueue = Queue{eltype(tree)}(min(nv(tree), 1024))

    ## Current and last visited vertices.
    v = (minimum ∘ children)(tree, mrca(tree))
    lastvisited = 0

    while !isempty(vstack)
        if iszero(v)
            peekv = first(vstack)
            isleaf(tree, peekv) && @goto enqueue

            peekv_maxchild = (maximum ∘ children)(tree, peekv)

            if peekv_maxchild ≠ lastvisited
                v = peekv_maxchild
            else
                @label enqueue
                enqueue!(vqueue, peekv)
                lastvisited = pop!(vstack)
            end
        else
            push!(vstack, v)
            v = isleaf(tree, v) ? 0 : (minimum ∘ children)(tree, v)
        end
    end

    vqueue
end

function (D::PhenotypeDensity)(tree::Tree)
    n = nleaves(tree)
    ## TODO: Manage this case more cleanly
    n <= 1 && error("tree must have at least two leaves")

    ## Queue of vertices to process.
    vqueue = postorder_traversal(tree)

    ## TODO: rewrite explanation.
    ## Stack used to store messages. Entry (i, j) contains the value
    ## of f(φ = g(i) | ψ = g(j), t) where φ and ψ are the source and
    ## destination of the message respectively and g is some functions
    ## that compute phenotype from an integer.
    mstack = Stack{Matrix{BigFloat}}(min(n, 1024))

    ## Necessary to keep track of the order.
    idstack = Stack{eltype(vqueue)}(min(n, 1024))

    copula = D.copula
    marginal_pdf = φ -> pdf_marginal(D.copula, φ)
    conditional_pdf = (φσ, φδ, Δlat) -> pdf_conditional(D.copula, φσ, φδ, Δlat)
    Φ = D.Φ

    for _ ∈ 1:ne(tree)
        v = dequeue!(vqueue)
        μ = cmatrix(tree, conditional_pdf, v, Φ)

        ## If v is an internal vertex, multiply μ by the K-R product
        ## of the messages of its children.
        if !isleaf(tree, v)
            idx, idy = pop!(idstack), pop!(idstack)
            x, y = pop!(mstack), pop!(mstack)

            μ = (x ⊙ y) * μ
        end

        push!(mstack, μ)
        push!(idstack, v)
    end

    ## Final message, computed with the marginal pdf.
    v = dequeue!(vqueue)
    (pop!(mstack) ⊙ pop!(mstack)) * cmatrix_root(marginal_pdf)
end

phenotypes(D::PhenotypeDensity) = D.Φ

"""
    cmatrix(tree, pdf, σ, Φ)
    cmatrix(tree, pdf, σ, φsσ, δ, φsδ)
    cmatrix_root(pdf)

Matrix containing relevant points of the conditional distribution of a vector
of phenotypes. Used by the belief propagation algorithm.
"""
function cmatrix end,
function cmatrix_root end

cmatrix_root(pdf) = reshape([pdf(φ) for φ in (false, true)], 2, 1)

function cmatrix(tree::Tree, pdf, σ, Φ)
    δ = dad(tree, σ)

    ## Assume that the phenotype of δ is unknown since it is an
    ## internal vertex.
    φsδ = (false, true)

    ## The phenotype of σ might be known if it is a leaf.
    if isleaf(tree, σ)
        φσ = Φ[σ]
        ismissing(φσ) && @goto unknownpheno2
        φsσ = [φσ]
    else # non-root internal vertex
        @label unknownpheno2
        φsσ = [false, true]
    end

    cmatrix(tree, pdf, σ, φsσ, δ, φsδ)
end

function cmatrix(tree::Tree, pdf, σ, φsσ, δ, φsδ)
    ## Julia issue #46331.
    ret = Matrix{Float64}(undef, length(φsσ), length(φsδ))
    Δlat = latitude(tree, δ) - latitude(tree, σ)

    @inbounds for (k, (φσ, φδ)) ∈ enumerate(Iterators.product(φsσ, φsδ))
        ret[k] = pdf(φσ, φδ, Δlat)
    end

    ret
end
