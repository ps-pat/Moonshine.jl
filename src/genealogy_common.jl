for G ∈ (:Tree, :Arg)
    Gcore = Symbol(string(G) * "Core")
    Gargname = Symbol(lowercase(string(G)))
    Gargnamecore = Symbol(string(Gargname) * "core")

    ## Random constructors.
    @eval function $Gcore(rng::AbstractRNG,
                          nmin::Integer, minlength::Integer,
                          nmax::Integer = 0, maxlength::Integer = 0;
                          genpars...)
        n = iszero(nmax) ? nmin : rand(rng, nmin:nmax)
        nmarkers = iszero(maxlength) ? minlength : rand(rng, minlength:maxlength)

        $Gcore([Sequence(rng, nmarkers) for _ ∈ 1:n]; genpars...)
    end

    @eval function $Gcore(nmin::Integer, minlength::Integer,
                          nmax::Integer = 0, maxlength::Integer = 0;
                          genpars...)
        $Gcore(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
    end

    @eval function $G(rng::AbstractRNG,
                      nmin::Integer, minlength::Integer,
                      nmax::Integer = 0, maxlength::Integer = 0;
                      genpars...)
        $G($Gcore(rng, nmin, minlength, nmax, maxlength; genpars...),
           zero(BigFloat))
    end

    @eval function $G(nmin::Integer, minlength::Integer,
                      nmax::Integer = 0, maxlength::Integer = 0;
                      genpars...)
        $G(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
    end

    ## Methods required by AbstractGenealogy interface that are only
    ## getters.
    for field ∈ (:(:graph), :(:latitudes), :(:sequences), :(:positions))
        fun_name = eval(field)
        @eval $fun_name($Gargname::$G) = getfield($(Gargname).core, $field)
    end

    ## Latitude of a vertex
    @eval function latitude($Gargname::$G, v)
        n = nleaves($Gargname)
        v ≤ n ? zero(Float64) : latitudes($Gargname)[v - n]
    end

    ## Probability of a genealogy.
    @eval function prob($Gargname::$G; logscale = false)
        ret = $(Gargname).logprob

        logscale ? ret : exp(ret)
    end

    ## Other methods.
    @eval isempty($Gargnamecore::$Gcore) = isempty($Gargnamecore.sequences)

    @eval isempty($Gargname::$G) = isempty($Gargname.core)

    for field ∈ [:(:seq_length), :(:Ne)]
        fun_name = eval(field)
        @eval begin
            export $fun_name

            function $fun_name($Gargname::$G)
                getfield($(Gargname).core, $field)
            end
        end
    end

    export mut_rate
    @eval mut_rate($Gargname::$G, scaled = true) =
        $Gargname.core.μloc * (scaled ? 4 * Ne($Gargname) : 1.)
end
