export lik_coal_is
"""
    lik_coal_is

Compute the likelihood of a phenotype using importance sampling.
"""
function lik_coal_is end

using Random: AbstractRNG, GLOBAL_RNG
using Statistics: mean
using Base.Threads
using Random: TaskLocalRNG
function lik_coal_is(haplotypes, phenotypes;
                     n = 1000, μ = 0, seq_length = 1,
                     α = t -> 1 - exp(-t), pars = Dict())
    dens_tree = CoalMutDensity(length(haplotypes), μ, seq_length)
    dens_φ = FrechetCoalDensity(phenotypes, α = α, pars = pars)

    terms = Vector{Function}(undef, n)
    @threads for k ∈ 1:n
        rng = TaskLocalRNG()
        tree = buildtree!(rng, Arg(haplotypes))

        dens_φ_marg = dens_φ(tree)
        c = dens_φ_marg([true]) + dens_φ_marg([false])
        logweight = dens_tree(tree, logscale = true) - tree.logprob
        terms[k] = φ -> exp(logweight) * dens_φ_marg(φ) / c
    end

    function(phenotypes; logscale = false)
        ret = mean(f -> f(phenotypes), terms)

        logscale ? log(ret) : ret
    end
end
