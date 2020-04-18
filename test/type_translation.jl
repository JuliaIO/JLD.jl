module Translation

using Random

const filename = joinpath(tempdir(), "translation-$(randstring()).jld")

module Writing

using JLD
import ..Translation: filename

mutable struct MyType
    a::Int
end

jldopen(filename, "w") do file
    truncate_module_path(file, Writing)
    write(file, "x", MyType(3))
end

end # Writing


module Reading

using JLD, Test
import ..Translation: filename

mutable struct MyType
    a::Int
    b::Float32
end

mutable struct MyOldType
    a::Int
end

translate("MyType", "Translation.Reading.MyOldType")

t = jldopen(filename) do file
    read(file, "x")
end

@test isa(t, MyOldType)

JLD.readas(x::MyOldType) = MyType(x.a, 0f0)

t = jldopen(filename) do file
    read(file, "x")
end

@test isa(t, MyType)

end # Reading

end
