abstract type AbstractGraphDensity <: Function end

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

#####################
# End of coalescent #
#####################

############################
# Coalescent with mutation #
############################

struct CoalMutDensity <:AbstractGraphDensity
    dens_coal::CoalDensity

    ## Scaled mutation rate.
    μ::BigFloat

    ## Length of sequences.
    seq_length::BigFloat

    CoalMutDensity(n, μ, seq_length) = new(CoalDensity(n), μ, seq_length)
end

function (D::CoalMutDensity)(arg::Arg, logscale = false)
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
# End of Coalescent with mutation #
###################################
