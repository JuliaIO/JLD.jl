module MyTypes

import Base: ==
export MyType, MyContainer

## Objects we want to save
# data in MyType is always of length 5, and that is the basis for a more efficient serialization
struct MyType{T}
    data::Vector{T}
    id::Int

    function MyType{T}(v::Vector{T}, id::Integer) where T
        length(v) == 5 || error("All vectors must be of length 5")
        new{T}(v, id)
    end
end
MyType(v::Vector{T}, id::Integer) where {T} = MyType{T}(v, id)
Base.eltype(::Type{MyType{T}}) where {T} = T
==(a::MyType, b::MyType) = a.data == b.data && a.id == b.id

struct MyContainer{T}
    objs::Vector{MyType{T}}
end
Base.eltype(::Type{MyContainer{T}}) where {T} = T
==(a::MyContainer, b::MyContainer) = length(a.objs) == length(b.objs) && all(i->a.objs[i]==b.objs[i], 1:length(a.objs))

end # MyTypes


### Here are the definitions needed to implement the custom serialization
# If you prefer, you could include these definitions in the MyTypes module
module MySerializer

using HDF5, JLD, ..MyTypes
using Compat

## Defining the serialization format
mutable struct MyContainerSerializer{T}
    data::Matrix{T}
    ids::Vector{Int}
end
MyContainerSerializer(data::Matrix{T},ids) where {T} = MyContainerSerializer{T}(data, ids)
Base.eltype(::Type{MyContainerSerializer{T}}) where {T} = T
Base.eltype(::MyContainerSerializer{T}) where {T} = T

JLD.readas(serdata::MyContainerSerializer) =
    MyContainer([MyType(serdata.data[:,i], serdata.ids[i]) for i = 1:length(serdata.ids)])
function JLD.writeas(data::MyContainer{T}) where T
    ids = [obj.id for obj in data.objs]
    n = length(data.objs)
    vectors = Matrix{T}(undef, 5, n)
    for i = 1:n
        vectors[:,i] = data.objs[i].data
    end
    MyContainerSerializer(vectors, ids)
end

end # MySerializer

using ..MyTypes, JLD, Compat.Test

obj1 = MyType(rand(5), 2)
obj2 = MyType(rand(5), 17)
container = MyContainer([obj1,obj2])
filename = joinpath(tempdir(), "customserializer-$(randstring()).jld")
jldopen(filename, "w") do file
    write(file, "mydata", container)
end

container_r = jldopen(filename) do file
    obj = file["mydata"]
    dtype = JLD.datatype(obj.plain)
    @test JLD.jldatatype(JLD.file(obj), dtype) === MySerializer.MyContainerSerializer{Float64}
    read(file, "mydata")
end

@test container_r == container
