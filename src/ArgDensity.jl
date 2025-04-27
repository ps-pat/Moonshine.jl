using DataStructures: enqueue!, dequeue!

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
