using JLD, Random, Test

primitive type ConvertedWritePrimitive 16 end

@testset "Converted array writes" begin
    mktempdir() do directory
        filename = joinpath(directory, "converted.jld")
        bits = rand(MersenneTwister(339), UInt16, 256, 256)
        values = copy(reinterpret(ConvertedWritePrimitive, bits))
        jldopen(filename, "w") do file
            write(file, "warm", values)
            @test read(file, "warm") == values
            # Warm the conversion, then exclude per-element boxing from the budget.
            @test (@allocated write(file, "measured", values)) < 8 * sizeof(values)
            restored = read(file, "measured")
            @test restored == values
            @test size(restored) == size(values)
            @test eltype(restored) === ConvertedWritePrimitive
        end
    end
end
