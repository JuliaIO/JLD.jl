using JLD, Random, Test

@testset "Float16 array conversion" begin
    mktempdir() do directory
        filename = joinpath(directory, "float16.jld")
        rng = MersenneTwister(337)
        arrays = [rand(rng, Float16, dims...) for dims in ((0,), (0, 2), (1,), (17, 19), (3, 5, 7))]
        large = rand(rng, Float16, 256, 256)
        jldopen(filename, "w") do file
            for (i, array) in enumerate(arrays)
                write(file, "array$i", array)
            end
            write(file, "large", large)
        end
        jldopen(filename, "r") do file
            for (i, array) in enumerate(arrays)
                restored = read(file, "array$i")
                @test restored == array
                @test size(restored) == size(array)
                @test eltype(restored) === Float16
            end
            @test read(file, "large") == large
            # Allow room for HDF5 metadata while ruling out per-element boxing.
            @test (@allocated read(file, "large")) < 8 * sizeof(large)
        end
    end
end
