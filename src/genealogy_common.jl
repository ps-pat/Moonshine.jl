for G ∈ (:Tree, :Arg)
    Gargname = Symbol(lowercase(string(G)))

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

    @eval function mrca($Gargname::$G, vs::AbstractVector, x;
                        buffer = default_buffer())
        μ = mrca($Gargname)
        iszero(μ) && return μ

        @no_escape buffer begin
            descendants_ptr = unsafe_convert(Ptr{VertexType},
                                             @alloc_ptr(nv(arg) * sizeof(VertexType)))
            while true
                flag = false
                for c ∈ children($Gargname, μ, x)
                    if vs ⊆ descendants!(descendants_ptr, $Gargname, c, x)
                        μ = c
                        flag = true
                        break
                    end
                end
                flag || break
            end
        end

        μ
    end

    @eval mrca($Gargname::$G, vs::AbstractVector) = mrca($Gargname, vs, -∞..∞)

    @eval tmrca($Gargname::$G, vs::AbstractVector, ωs) =
        latitude($Gargname, mrca($Gargname, vs, ωs))
end

