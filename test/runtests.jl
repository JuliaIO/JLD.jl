using HDF5

runtest(filename) = (println(filename); include(filename))

runtest("jldtests.jl")
runtest("require.jl")
runtest("custom_serialization.jl")
runtest("type_translation.jl")
if Pkg.installed("DataFrames") != nothing
    runtest("jld_dataframe.jl")
end
