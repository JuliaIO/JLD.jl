module Translation

const filename = joinpath(tempdir(), "translation.jld")

module Writing

using JLD
import ..Translation: filename

type MyType
    a::Int
end

jldopen(filename, "w") do file
    truncate_module_path(file, Writing)
    write(file, "x", MyType(3))
end

end


module Reading

using JLD, Base.Test
import ..Translation: filename

type MyType
    a::Int
    b::Float32
end

type MyOldType
    a::Int
end

translate("MyType", "MyOldType")

t = jldopen(filename) do file
    truncate_module_path(file, Reading)
    read(file, "x")
end

@test isa(t, MyOldType)

JLD.readas(x::MyOldType) = MyType(x.a, 0f0)

t = jldopen(filename) do file
    read(file, "x")
end

@test isa(t, MyType)

end

end
