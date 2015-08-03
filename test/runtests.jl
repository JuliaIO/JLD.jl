using HDF5

runtest(filename) = (println(filename); include(filename))

runtest("jldtests.jl")
runtest("require.jl")
if Pkg.installed("DataFrames") != nothing
    runtest("jld_dataframe.jl")
end
