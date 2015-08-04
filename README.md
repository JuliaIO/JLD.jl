# Saving and loading variables in Julia Data format (JLD)

[![Build Status](https://travis-ci.org/JuliaLang/JLD.jl.svg?branch=master)](https://travis-ci.org/JuliaLang/JLD.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/p2wkp87jfr0q7g4q?svg=true)](https://ci.appveyor.com/project/timholy/jld-jl)

JLD, for which files conventionally have the extension `.jld`, is a
widely-used format for data storage with the Julia programming
language.  JLD is a specific "dialect" of [HDF5][HDF5], a
cross-platform, multi-language data storage format most frequently
used for scientific data.  By comparison with "plain" HDF5, JLD files
automatically add attributes and naming conventions to preserve type
information for each object.

For lossless storage of arbitrary Julia objects, the only other
complete solution appears to be Julia's serializer, which can be
accessed via the `serialize` and `deserialize` commands.  However,
because the serializer is also used for inter-process communication,
long-term backwards compatibility is currently uncertain.  (The
[JLDArchives](https://github.com/timholy/JLDArchives.jl) package exists to test compatibility of older JLD file
formats.) If you choose to save data using the serializer, please use
the file extension `.jls` to distinguish the files from `.jld` files.


## Installation

Within Julia, use the package manager:
```julia
Pkg.add("JLD")
```

Currently this also requires the [HDF5 package](https://github.com/timholy/HDF5.jl).

## Quickstart

To use the JLD module, begin your code with

```julia
using JLD
```

If you just want to save a few variables and don't care to use the more
advanced features, then a simple syntax is:

```
t = 15
z = [1,3]
save("/tmp/myfile.jld", "t", t, "arr", z)
```
Here we're explicitly saving `t` and `z` as `"t"` and `"arr"` within
`myfile.jld`. You can alternatively pass `save` a dictionary; the keys must be
strings and are saved as the variable names of their values within the JLD
file. You can read these variables back in with
```
d = load("/tmp/myfile.jld")
```
which reads the entire file into a returned dictionary `d`. Or you can be more
specific and just request particular variables of interest. For example, `z =
load("/tmp/myfile.jld", "arr")` will return the value of `arr` from the file
and assign it back to z.

There are also convenience macros `@save` and `@load` that work on the
variables themselves. `@save "/tmp/myfile.jld" t z` will create a file with
just `t` and `z`; if you don't mention any variables, then it saves all the
variables in the current module. Conversely, `@load` will pop the saved
variables directly into the global workspace of the current module.
However, keep in mind that these macros have significant limitations: for example,
you can't use `@load` inside a function, there are constraints on using string
interpolation inside filenames, etc. These limitations stem
from the fact that Julia compiles functions to machine code before evaluation,
so you cannot introduce new variables at runtime or evaluate expressions
in other workspaces.
The `save` and `load` functions do not have these limitations, and are therefore
recommended as being considerably more robust, at the cost of some slight
reduction in convenience.

More fine-grained control can be obtained using functional syntax:

```julia
jldopen("mydata.jld", "w") do file
    write(file, "A", A)  # alternatively, say "@write file A"
end

c = jldopen("mydata.jld", "r") do file
    read(file, "A")
end
```
This allows you to add variables as they are generated to an open JLD file.
You don't have to use the `do` syntax (`file = jldopen("mydata.jld", "w")` works
just fine), but an advantage is that it will automatically close the file (`close(file)`)
for you, even in cases of error.

Julia's high-level wrapper, providing a dictionary-like interface, may
also be of interest:

```julia
using JLD, HDF5

jldopen("test.jld", "w") do file
    g = g_create(file, "mygroup") # create a group
    g["dset1"] = 3.2              # create a scalar dataset inside the group
    g["dest2"] = rand(2,2)
end
```

Note that the features of HDF5 generally can also be used on JLD files.

## Types and their definitions

You can save objects that have user-defined type; in a fresh Julia session, before loading those objects these types need to be defined. If no definition is avaiable, the JLD module will automatically create the types for you.  However, it's important to note that `MyType`, defined automatically by JLD, is not the same `MyType` as defined in an external module---in particular, module functions will not work for types defined by JLD.  To ensure that loaded types have the full suite of behaviors provided by their definition in external modules, you should ensure that such modules are available before reading such variables from a `.jld` file.

To ensure automatic loading of modules, use `addrequire` to specify any dependencies. For example, suppose you have a file `"MyTypes.jl"` somewhere on your default `LOAD_PATH`, defined this way:
```
module MyTypes

export MyType

type MyType
    value::Int
end

end
```
and you have an object `x` of type `MyType`. Then save `x` in the following way:

```
jldopen("somedata.jld", "w") do file
    addrequire(file, MyTypes)
    write(file, "x", x)
end
```
This will cause `"MyTypes.jl"` to be loaded automatically whenever `"somedata.jld"` is opened.

## Complete documentation

More extensive documentation, including information about the JLD
format conventions, can be found in the [`doc/`](doc/) directory.

The `test/` directory contains a number of test scripts that also
demonstrate usage.

## Credits

- Simon Kornblith and Tim Holy (co-maintainers and primary authors)

- Tom Short contributed to string->type conversion

- Thanks also to the users who have reported bugs and tested fixes


[Julia]: http://julialang.org "Julia"
[HDF5]: http://www.hdfgroup.org/HDF5/ "HDF5"
