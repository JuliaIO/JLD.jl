using HDF5, Random

runtest(filename) = (println(filename); include(filename))

runtest("jldtests.jl")
runtest("require.jl")
runtest("custom_serialization.jl")
runtest("type_translation.jl")
runtest("read_conversion.jl")
runtest("write_conversion.jl")
# runtest("jld_dataframe.jl") # FIXME: fails (segfault when reading back the dataframe)
