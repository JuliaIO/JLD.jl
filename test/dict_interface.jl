using JLD, Test

# Regression tests for the AbstractDict `iterate` contract that JldFile/JldGroup now
# inherit transitively via HDF5.H5DataStore <: AbstractDict{String,Any} (HDF5.jl >= 0.18).
# `iterate` must yield `key => value` pairs, not bare values, or generic AbstractDict
# operations built on it (pairs, values, Dict(x), filter, ...) silently misbehave.
mktempdir() do d
    fn = joinpath(d, "dict_interface.jld")

    jldopen(fn, "w") do fid
        fid["a"] = 1
        fid["b"] = 2
        g = create_group(fid, "g")
        g["x"] = 10
        g["y"] = 20
    end

    jldopen(fn, "r") do fid
        @testset "iterate yields pairs" begin
            # getindex on a JldFile/JldGroup is lazy: it returns a JldDataset/JldGroup
            # handle, not the materialized value, so `v` is one of those handle types.
            for (k, v) in fid
                @test k isa String
                @test v isa Union{JLD.JldDataset,JLD.JldGroup}
            end
        end

        @testset "Dict(jldfile) / Dict(jldgroup)" begin
            d = Dict(k => (v isa JLD.JldDataset ? read(v) : v) for (k, v) in fid)
            @test d["a"] == 1
            @test d["b"] == 2
            @test d["g"] isa JLD.JldGroup

            g = fid["g"]
            dg = Dict(k => read(v) for (k, v) in g)
            @test dg == Dict("x" => 10, "y" => 20)
        end

        @testset "collect(...) over a group" begin
            # Deliberately tests plain `collect`/iteration directly, not `pairs(g)`:
            # `pairs(::AbstractDict) === d` only kicks in once JldGroup is actually
            # recognized as an AbstractDict (HDF5.jl >= 0.18); prior to that, generic
            # `pairs` falls back to a key/value zip that would double-wrap this method's
            # own pair-yielding `iterate`. What's being guarded here is `iterate` itself,
            # which must behave the same regardless of the HDF5.jl version in use.
            g = fid["g"]
            ps = collect(g)
            @test Set(first.(ps)) == Set(["x", "y"])
            @test Dict(k => read(v) for (k, v) in ps) == Dict("x" => 10, "y" => 20)
        end
    end

    @testset "delete! on a group still removes datasets" begin
        jldopen(fn, "r+") do fid
            g = fid["g"]
            @test haskey(g, "x")
            delete!(g)
        end
        jldopen(fn, "r") do fid
            @test !haskey(fid, "g")
        end
    end
end
