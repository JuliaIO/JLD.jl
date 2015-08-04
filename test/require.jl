using HDF5, JLD

module JLDTemp
using HDF5, JLD
include("JLDTest.jl")  # Private module inside JLDTemp

function create()
    x = JLDTest.Object(convert(Int16, 5))  # int16 makes this work on 0.2
    jldopen("require.jld", "w") do file
        addrequire(file, :JLDTest)
        truncate_module_path(file, JLDTemp)
        write(file, "x", x)
    end
end
end

JLDTemp.create()

push!(LOAD_PATH, splitdir(@__FILE__)[1])
x = jldopen("require.jld") do file
    read(file, "x")
end
@assert typeof(x) == JLDTest.Object
@assert x.data == 5
pop!(LOAD_PATH)
rm("require.jld")
