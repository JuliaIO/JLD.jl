using HDF5

runtest(filename) = (println(filename); include(filename))

runtest("jld.jl")
runtest("require.jl")
if Pkg.installed("DataFrames") != nothing
    runtest("jld_dataframe.jl")
end
