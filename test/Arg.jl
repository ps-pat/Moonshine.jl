using Test

arg_qc = Quickcheck("Test ARGs Properties", serialize_fails = false)

const A2 = Arg2{UInt}

let ftree = CoalDensity(2),
    α = t -> -expm1(-t),
    pars = Dict(:p => 0.15),
    φs_false = [missing, false],
    φs_true = [missing, true],
    fφ_false = FrechetCoalDensity(φs_false, α = α, pars = pars),
    fφ_true = FrechetCoalDensity(φs_true, α = α, pars = pars)

    @add_predicate(arg_qc, "Weights equal for 2 leaves.",
                   a2::A2 -> ≈(ftree(a2.arg, logscale = true),
                               a2.arg.logprob,
                               atol = 1e-5))

    # @add_predicate(arg_qc, "Density for 2 leaves (φ = 0)",
    #                a2::A2 -> )
end

@testset "ARG" begin
    @quickcheck arg_qc
end
