for G ∈ (:Tree, :Arg)
    Gargname = Symbol(lowercase(string(G)))

    ## Random constructors.
    @eval function $G(rng::AbstractRNG, n, μ, ρ, Ne, sequence_length)
        sample = Sample(rng, n, μ, ρ, Ne, sequence_length)
        $G(sample)
    end

    ## Methods required by AbstractGenealogy interface that are only
    ## getters.
    for field ∈ (:(:graph), :(:latitudes), :(:sequences))
        fun_name = eval(field)
        @eval $fun_name($Gargname::$G) = getfield($Gargname, $field)
    end

    ## Positions of the markers.
    @eval positions($Gargname::$G) = positions(sam($Gargname))

    ## Sample of genetic sequences.
    export sam
    @eval sam($Gargname::$G) = getfield($Gargname, :sample)

    ## Latitude of a vertex
    @eval function latitude($Gargname::$G, v)
        n = nleaves($Gargname)
        v ≤ n ? zero(Float64) : latitudes($Gargname)[v - n]
    end

    ## Probability of a genealogy.
    @eval function prob($Gargname::$G; logscale = false)
        ret = getfield($Gargname, :logprob)[]

        logscale ? ret : exp(ret)
    end

    ## Other methods.
    @eval isempty($Gargname::$G) = isempty(getfield($Gargname, :sequences))

    for par ∈ ("mut", "rec")
        f = Symbol(par * "_rate")
        @eval export $f

        @eval $f($Gargname::$G, scaled = true) = $f(sam($Gargname), scaled)
    end

# -- MRCA --------------------------------------------------------------

    @eval function mrca($Gargname::$G)
    any(iszero, latitudes($Gargname)) && return zero(VertexType)
    isone(nv($Gargname)) && return one(VertexType)
    argmax(latitudes($Gargname)) + nleaves($Gargname)
    end

    @eval function mrca($Gargname::$G, vs)
        μ = mrca($Gargname)
        iszero(μ) && return μ

        while true
            @label beg
            for c ∈ children($Gargname, μ)
                vs ⊆ descendants($Gargname, c) || continue
                μ = c
                @goto beg
            end
            break
        end

        μ
    end
end
