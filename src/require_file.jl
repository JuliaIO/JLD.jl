function require_file(fn::AbstractString)
    dir, bn = splitdir(fn)
    mod = symbol(endswith(bn, ".jl") ? bn : splitext(bn)[1])
    if !in(dir, LOAD_PATH)
        push!(LOAD_PATH, dir)
        @eval import $mod
        pop!(LOAD_PATH)
    else
        @eval import $mod
    end
end
