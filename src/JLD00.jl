###############################################
## Reading and writing Julia data .jld files ##
###############################################

module JLD00
using Printf
using HDF5
# Add methods to...
import HDF5: file, create_group, open_group, delete_object, name, ismmappable, readmmap
import Base: close, dump, getindex, iterate, length, read, setindex!, size, show, delete!, write
import ..JLD: JLD, _joinpath

# See julia issue #8907
replacements = Any[]
push!(replacements, :(s = replace(s, r"Uint(?=\d{1,3})" => "UInt")))
push!(replacements, :(s = replace(s, r"ASCIIString|UTF8String|ByteString" => "String")))
ex = Expr(:block, replacements...)
@eval function julia_type(s::AbstractString)
    $ex
    _julia_type(s)
end

const magic_base = "Julia data file (HDF5), version "
const version_current = "0.0.2"
const pathrefs = "/_refs"
const pathtypes = "/_types"
const pathrequire = "/_require"
const name_type_attr = "julia type"

### Dummy types used for converting attribute strings to Julia types
mutable struct UnsupportedType; end
mutable struct UnconvertedType; end
mutable struct CompositeKind; end   # here this means "a type with fields"

# The Julia Data file type
# Purpose of the nrefs field:
# length(group) only returns the number of _completed_ items in a group. Since
# we'll write recursively, we need to keep track of the number of reference
# objects _started_.
mutable struct JldFile <: HDF5.H5DataStore
    plain::HDF5.File
    version::String
    toclose::Bool
    writeheader::Bool
    mmaparrays::Bool

    function JldFile(plain::HDF5.File, version::AbstractString=version_current, toclose::Bool=true,
                     writeheader::Bool=false, mmaparrays::Bool=false)
        f = new(plain, version, toclose, writeheader, mmaparrays)
        if toclose
            finalizer(close, f)
        end
        f
    end
end

struct JldGroup
    plain::HDF5.Group
    file::JldFile
end

struct JldDataset
    plain::HDF5.Dataset
    file::JldFile
end

file(x::JldFile) = x
file(x::Union{JldGroup, JldDataset}) = x.file

function close(f::JldFile)
    if f.toclose
        close(f.plain)
        if f.writeheader
            magic = zeros(UInt8, 512)
            tmp = string(magic_base, f.version)
            magic[1:length(tmp)] = Vector{UInt8}(codeunits(tmp))
            rawfid = open(f.plain.filename, "r+")
            write(rawfid, magic)
            close(rawfid)
        end
        f.toclose = false
    end
    nothing
end
close(g::Union{JldGroup, JldDataset}) = close(g.plain)
show(io::IO, fid::JldFile) = isvalid(fid.plain) ? print(io, "Julia data file version ", fid.version, ": ", fid.plain.filename) : print(io, "Julia data file (closed): ", fid.plain.filename)

function jldopen(filename::AbstractString, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool; mmaparrays::Bool=false)
    local fj
    if ff && !wr
        error("Cannot append to a write-only file")
    end
    if !cr && !isfile(filename)
        error("File ", filename, " cannot be found")
    end
    version = version_current
    pa = create_property(HDF5.H5P_FILE_ACCESS)
    try
        pa[:fclose_degree] = HDF5.H5F_CLOSE_STRONG
        if cr && (tr || !isfile(filename))
            # We're truncating, so we don't have to check the format of an existing file
            # Set the user block to 512 bytes, to save room for the header
            p = create_property(HDF5.H5P_FILE_CREATE)
            local f
            try
                p[:userblock] = 512
                f = HDF5.h5f_create(filename, HDF5.H5F_ACC_TRUNC, p.id, pa.id)
            finally
                close(p)
            end
            fj = JldFile(HDF5.File(f, filename), version, true, true, mmaparrays)
            # initialize empty require list
            write(fj, pathrequire, String[])
        else
            # Test whether this is a jld file
            sz = filesize(filename)
            if sz < 512
                error("File size indicates $filename cannot be a Julia data file")
            end
            magic = Vector{UInt8}(undef, 512)
            rawfid = open(filename, "r")
            try
                magic = read!(rawfid, magic)
            finally
                close(rawfid)
            end
            if length(magic) ≥ ncodeunits(magic_base) && view(magic, 1:ncodeunits(magic_base)) == Vector{UInt8}(codeunits(magic_base))
                f = HDF5.h5f_open(filename, wr ? HDF5.H5F_ACC_RDWR : HDF5.H5F_ACC_RDONLY, pa.id)
                version = unsafe_string(pointer(magic) + length(magic_base))
                fj = JldFile(HDF5.File(f, filename), version, true, true, mmaparrays)
                # Load any required files/packages
                if haskey(fj, pathrequire)
                    r = read(fj, pathrequire)
                    for fn in r
                        require(fn)
                    end
                end
            else
                if HDF5.ishdf5(filename)
                    println("$filename is an HDF5 file, but it is not a recognized Julia data file. Opening anyway.")
                    # inverse of logic in jldopen() below
                    mode =  rd & !wr & !cr & !tr & !ff ? "r" :
                            rd &  wr & !cr & !tr & !ff ? "r+" :
                           !rd &  wr &  cr &  tr & !ff ? "w" : error("invalid mode")
                    fj = JldFile(h5open(filename, mode), version_current, true, false, mmaparrays)
                else
                    error("$filename does not seem to be a Julia data or HDF5 file")
                end
            end
        end
    finally
        close(pa)
    end
    return fj
end

function jldopen(fname::AbstractString, mode::AbstractString="r"; mmaparrays::Bool=false)
    mode == "r"  ? jldopen(fname, true , false, false, false, false, mmaparrays=mmaparrays) :
    mode == "r+" ? jldopen(fname, true , true , false, false, false, mmaparrays=mmaparrays) :
    mode == "w"  ? jldopen(fname, false, true , true , true , false, mmaparrays=mmaparrays) :
#     mode == "w+" ? jldopen(fname, true , true , true , true , false) :
#     mode == "a"  ? jldopen(fname, false, true , true , false, true ) :
#     mode == "a+" ? jldopen(fname, true , true , true , false, true ) :
    error("invalid open mode: ", mode)
end

function jldopen(f::Function, args...)
    jld = jldopen(args...)
    try
        f(jld)
    finally
        close(jld)
    end
end

function jldobject(obj_id::HDF5.hid_t, parent)
    obj_type = HDF5.h5i_get_type(obj_id)
    obj_type == HDF5.H5I_GROUP ? JldGroup(HDF5.Group(obj_id, file(parent.plain)), file(parent)) :
    obj_type == HDF5.H5I_DATATYPE ? HDF5.Datatype(obj_id) :
    obj_type == HDF5.H5I_DATASET ? JldDataset(HDF5.Dataset(obj_id, file(parent.plain)), file(parent)) :
    error("Invalid object type for path ", path)
end

getindex(parent::Union{JldFile, JldGroup}, path::String) =
    jldobject(HDF5.h5o_open(parent.plain.id, path, HDF5.H5P_DEFAULT), parent)

function getindex(parent::Union{JldFile, JldGroup, JldDataset}, r::HDF5.Reference)
    r == HDF5.Reference() && error("Reference is null")
    obj_id = HDF5.h5r_dereference(parent.plain.id, HDF5.H5P_DEFAULT, HDF5.H5R_OBJECT, r)
    jldobject(obj_id, parent)
end

### "Inherited" behaviors
create_group(parent::Union{JldFile, JldGroup}, args...) = JldGroup(create_group(parent.plain, args...), file(parent))
function create_group(f::Function, parent::Union{JldFile, JldGroup}, args...)
    g = JldGroup(create_group(parent.plain, args...), file(parent))
    try
        f(g)
    finally
        close(g)
    end
end
open_group(parent::Union{JldFile, JldGroup}, args...) = JldGroup(open_group(parent.plain, args...), file(parent))
name(p::Union{JldFile, JldGroup, JldDataset}) = name(p.plain)
Base.haskey(p::Union{JldFile, JldGroup, JldDataset}, path::String) = haskey(p.plain, path)
root(p::Union{JldFile, JldGroup, JldDataset}) = open_group(file(p), "/")
delete_object(parent::Union{JldFile, JldGroup}, args...) = delete_object(parent.plain, args...)
function ensurepathsafe(path::String)
    if any([startswith(path, s) for s in (pathrefs,pathtypes,pathrequire)])
        error("$name is internal to the JLD format, use delete_object if you really want to delete it")
    end
end
function delete!(o::JldDataset)
    fullpath = name(o)
    ensurepathsafe(fullpath)
    delete_object(o.file, fullpath)
    refspath = _joinpath(pathrefs, fullpath[2:end])
    haskey(o.file, refspath) && delete_object(o.file, refspath)
end
function delete!(g::JldGroup)
    fullpath = name(g)
    ensurepathsafe(fullpath)
    for o in g typeof(o) == JldDataset && delete!(o) end
    delete_object(g.file,name(g))
end
function delete!(parent::Union{JldFile, JldGroup}, path::String)
    haskey(parent, path) || error("$path does not exist in $parent")
    delete!(parent[path])
end
delete!(parent::Union{JldFile, JldGroup}, args::Tuple{Vararg{String}}) = for a in args delete!(parent,a) end
ismmappable(obj::JldDataset) = ismmappable(obj.plain)
readmmap(obj::JldDataset, args...) = readmmap(obj.plain, args...)
setindex!(parent::Union{JldFile, JldGroup}, val, path::String) = write(parent, path, val)

Base.iterate(parent::Union{JldFile, JldGroup}, state=(keys(parent), 1)) = state[2] > length(state[1]) ? nothing :
                                                     (parent[state[1][state[2]]], (state[1], state[2]+1))


### Julia data file format implementation ###


### Read ###

function read(parent::Union{JldFile, JldGroup}, name::String)
    local val
    obj = parent[name]
    try
        val = read(obj)
    finally
        close(obj)
    end
    val
end
read(parent::Union{JldFile,JldGroup}, name::Symbol) = read(parent, string(name))

# read and readsafely differ only in how they handle CompositeKind
function read(obj::Union{JldFile, JldDataset})
    if !haskey(attributes(obj.plain), name_type_attr)
        # Fallback to plain read
        return read(obj.plain)
    end
    # Read the type
    typename = read_attribute(obj.plain, name_type_attr)
    if typename == "Tuple"
        return read_tuple(obj)
    end
    # Convert to Julia type
    T = julia_type(typename)
    if T == CompositeKind
        # Use type information in the file to ensure we find the right module
        typename = read_attribute(obj.plain, "CompositeKind")
        try
            gtypes = root(obj)[pathtypes]
            try
                objtype = gtypes[typename]
                try
                    modnames = read_attribute(objtype.plain, "Module")
                    mod = Main
                    for mname in modnames
                        mod = Core.eval(mod, Symbol(mname))
                    end
                    T = Core.eval(mod, Symbol(typename))
                finally
                    close(objtype)
                end
            finally
                close(gtypes)
            end
        catch
            error("Type ", typename, " is not recognized. As a fallback, you can load ", name(obj), " with readsafely().")
        end
    end
    if T == UnsupportedType
        error("Type expression ", typename, ", is not recognized.")
    end
    read(obj, T)
end
function readsafely(obj::Union{JldFile, JldDataset})
    if !haskey(attributes(obj.plain), name_type_attr)
        # Fallback to plain read
        return read(obj.plain)
    end
    # Read the type
    typename = read_attribute(obj.plain, name_type_attr)
    if typename == "Tuple"
        return read_tuple(obj)
    end
    # Convert to Julia type
    T = julia_type(typename)
    local ret
    if T == CompositeKind
        # Read as a dict
        typename = read_attribute(obj.plain, "CompositeKind")
        gtypes = root(obj)[pathtypes]
        try
            objtype = gtypes[typename]
            try
                n = read(objtype)
                v = getrefs(obj, Any)
                ret = Dict(zip(n[1,:], v))
            finally
                close(objtype)
            end
        finally
            close(gtypes)
        end
    else
        ret = read(obj, T)
    end
    ret
end
function readsafely(parent::Union{JldFile, JldGroup}, name::String)
    local ret
    obj = parent[name]
    try
        ret = readsafely(JldDataset(obj.plain, obj.file))
    finally
        close(obj)
    end
    return ret
end
readsafely(parent::Union{JldFile,JldGroup}, name::Symbol) = readsafely(parent, string(symbol))

# Basic types
const BitsKindOrString = Union{HDF5.BitsType, String}
read(obj::JldDataset, ::Type{T}) where {T<:BitsKindOrString} = read(obj.plain, T)
function read(obj::JldDataset, ::Type{Array{T}}) where T<:HDF5.BitsType
    A = obj.file.mmaparrays && HDF5.iscontiguous(obj.plain) ? readmmap(obj.plain, T) : read(obj.plain, T)
    if isempty(A) && haskey(obj, "dims")
        dims = read_attribute(obj.plain, "dims")
        A = size(A) == () ? Array{T}(undef, dims...) : reshape(A, dims...)
    end
    A
end
function read(obj::JldDataset, ::Type{Array{T}}) where {T<:String}
    A = read(obj.plain, T)
    if isempty(A) && haskey(obj, "dims")
        dims = read_attribute(obj.plain, "dims")
        A = size(A) == () ? Array{T}(undef, dims...) : reshape(A, dims...)
    end
    A
end
function read(obj::JldDataset, ::Type{Array{T,N}}) where {T<:BitsKindOrString,N}
    A = read(obj.plain, T)
    if isempty(A) && haskey(obj, "dims")
        dims = read_attribute(obj.plain, "dims")
        A = size(A) == () ? Array{T}(undef, dims...) : reshape(A, dims...)
    end
    A
end

# Arrays-of-arrays of basic types
function read(obj::JldDataset, ::Type{Array{Array{T,N},M}}) where {T<:HDF5.BitsType,M,N}
    # fallback for backwards compatibility with pre-v0.2.27 format
    HDF5.get_jl_type(datatype(obj.plain)) == HDF5.Reference &&
        return getrefs(obj, Array{T,N})
    A = read(obj.plain, HDF5.VariableArray{T})
    if isempty(A) && haskey(obj, "dims")
        dims = read_attribute(obj.plain, "dims")
        A = size(A) == () ? Array{T}(undef, dims...) : reshape(A, dims...)
    end
    convert(Array{Array{T,N},M}, A)
end

# Nothing
read(obj::JldDataset, ::Type{Nothing}) = nothing
read(obj::JldDataset, ::Type{Bool}) = read(obj, UInt8) != 0

# Types
read(obj::JldDataset, ::Type{Type{T}}) where {T} = T

# Bool
function read(obj::JldDataset, ::Type{Array{Bool,N}}) where N
    format = read_attribute(obj.plain, "julia_format")
    if format == "EachUint8"
        bool(read(obj.plain, UInt8))
    else
        error("bool format not recognized")
    end
end

# Complex
function read(obj::JldDataset, ::Type{Complex{T}}) where T
    a = read(obj.plain, T)
    a[1]+a[2]*im
end
function read(obj::JldDataset, ::Type{Array{T,N}}) where {T<:Complex,N}
    A = read(obj, realtype(T))
    reshape(reinterpret(T, vec(A)), ntuple(i->size(A, i+1), ndims(A)-1))
end

# Symbol
read(obj::JldDataset, ::Type{Symbol}) = Symbol(read(obj.plain, String))
read(obj::JldDataset, ::Type{Array{Symbol,N}}) where {N} = map(Symbol, read(obj.plain, String))

# Char
read(obj::JldDataset, ::Type{Char}) = Char(read(obj.plain, UInt32))

# General arrays
read(obj::JldDataset, t::Type{Array{T,N}}) where {T,N} = getrefs(obj, T)

# Tuple
function read_tuple(obj::JldDataset)
    t = read(obj, Array{Any, 1})
    return tuple(t...)
end
function read_tuple(obj::JldDataset, indices::AbstractVector)
    t = read(obj, Array{Any, 1})
    return tuple(t...)
end

# Dict
function read(obj::JldDataset, ::Type{T}) where T<:AbstractDict
    kv = getrefs(obj, Any)
    ret = T()
    for (cn, c) in zip(kv[1], kv[2])
        ret[cn] = c
    end
    ret
end

# Expressions
function read(obj::JldDataset, ::Type{Expr})
    a = getrefs(obj, Any)
    Expr(a[1], a[2]...)
end

# CompositeKind
read(obj::JldDataset, T::UnionAll) = read(obj, Base.unwrap_unionall(T))
function read(obj::JldDataset, T::DataType)
    if isempty(fieldnames(T)) && T.size > 0
        return read_bitstype(obj, T)
    end
    local x
    # Add the parameters
    if haskey(obj, "TypeParameters")
        params = read_attribute(obj.plain, "TypeParameters")
        if !isempty(params)
            p = Vector{Any}(undef, length(params))
            for i = 1:length(params)
                p[i] = Core.eval(@__MODULE__, Meta.parse(params[i]))
            end
            T = T.name.wrapper
            T = T{p...}
        end
    end
    v = getrefs(obj, Any)
    if length(v) == 0
        x = ccall(:jl_new_struct, Any, (Any,Any...), T)
    else
        n = fieldnames(T)
        if length(v) != length(n)
            error("Wrong number of fields")
        end
        if !T.mutable
            x = ccall(:jl_new_structv, Any, (Any,Ptr{Cvoid},UInt32), T, v, length(fieldnames(T)))
        else
            x = ccall(:jl_new_struct_uninit, Any, (Any,), T)
            for i = 1:length(v)
                if isassigned(v, i)
                    setfield!(x, n[i], v[i])
                end
            end
        end
    end
    x
end

function read_bitstype(obj::JldDataset, T::DataType)
    a = read(obj.plain)
    reinterpret(T, a[1])
end

# Read an array of references
function getrefs(obj::JldDataset, ::Type{T}) where T
    refs = read(obj.plain, HDF5.Reference)
    if isempty(refs) && size(refs) == () # HDF5.EmptyArray-like
        out = Array{T}(undef, 0)
    else
        out = Array{T}(undef, size(refs))
    end
    f = file(obj)
    for i = 1:length(refs)
        if refs[i] != HDF5.Reference()
            ref = f[refs[i]]
            try
                out[i] = read(ref)
            finally
                close(ref)
            end
        end
    end
    return out
end
function getrefs(obj::JldDataset, ::Type{T}, indices::Union{Integer, AbstractVector}...) where T
    refs = read(obj.plain, HDF5.Reference)
    refs = refs[indices...]
    f = file(obj)
    local out
    if isa(refs, HDF5.Reference)
        # This is a scalar, not an array
        ref = f[refs]
        try
            out = read(ref)
        finally
            close(ref)
        end
    else
        out = Array{T}(undef, size(refs))
        for i = 1:length(refs)
            ref = f[refs[i]]
            try
                out[i] = read(ref)
            finally
                close(ref)
            end
        end
    end
    return out
end

### Writing ###

# Write "basic" types
function write(parent::Union{JldFile, JldGroup}, name::String,
                   data::Union{T, StridedArray{T}}, astype::String) where T<:Union{HDF5.BitsType, String}
    # Create the dataset
    dset, dtype = create_dataset(parent.plain, name, data)
    try
        # Write the attribute
        write_attribute(dset, name_type_attr, astype)
        isa(data, StridedArray) && isempty(data) && write_attribute(dset, "dims", [size(data)...])
        # Write the data
        HDF5.writearray(dset, dtype.id, data)
    finally
        close(dset)
        close(dtype)
    end
end
write(parent::Union{JldFile, JldGroup}, name::String, data::Union{T, Array{T}}) where {T<:Union{HDF5.BitsType, String}} =
    write(parent, name, data, full_typename(typeof(data)))

# Arrays-of-arrays of basic types
write(parent::Union{JldFile, JldGroup}, name::String,
            data::Array{Array{T,1}}, astype::String) where {T<:Union{HDF5.BitsType, String}} =
    write(parent, name, HDF5.VLen(data), astype)
write(parent::Union{JldFile, JldGroup}, name::String,
            data::Array{Array{T,1}}) where {T<:Union{HDF5.BitsType, String}} =
    write(parent, name, data, full_typename(typeof(data)))
function write(parent::Union{JldFile, JldGroup}, name::String,
                  data::HDF5.VLen{T}, astype::String) where T
    # Create the dataset
    dset, dtype = create_dataset(parent.plain, name, data)
    try
        # Write the attribute
        write_attribute(dset, name_type_attr, astype)
        isa(data, Array) && isempty(data) && write_attribute(dset, "dims", [size(data)...])
        # Write the data
        HDF5.writearray(dset, dtype.id, data)
    finally
        close(dset)
        close(dtype)
    end
end


# Write nothing
function write(parent::Union{JldFile, JldGroup}, name::String, n::Nothing, astype::String)
    local dspace, dset
    try
        dspace = dataspace(nothing)
        dset = HDF5.Dataset(HDF5.h5d_create(HDF5.parents_create(HDF5.checkvalid(parent.plain), name, HDF5.H5T_NATIVE_UINT8, dspace.id,
                           HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT, HDF5.H5P_DEFAULT)...), file(parent.plain))
        write_attribute(dset, name_type_attr, astype)
    finally
        close(dspace)
        close(dset)
    end
end
write(parent::Union{JldFile, JldGroup}, name::String, n::Nothing) = write(parent, name, n, "Nothing")

# Types
# the first is needed to avoid an ambiguity warning
write(parent::Union{JldFile, JldGroup}, name::String, t::Tuple{Vararg{Type{T}}}) where {T} =
    write(parent, name, Any[t...], "Tuple")
write(parent::Union{JldFile, JldGroup}, name::String, t::Type{T}) where {T} =
    write(parent, name, nothing, string("Type{", full_typename(t), "}"))

# Bools
write(parent::Union{JldFile, JldGroup}, name::String, tf::Bool) = write(parent, name, uint8(tf), "Bool")
function write(parent::Union{JldFile, JldGroup}, name::String, tf::Array{Bool})
    write(parent, name, uint8(tf), full_typename(typeof(tf)))
    write_attribute(parent[name].plain, "julia_format", "EachUint8")
end

# Complex
realtype(::Type{Complex{T}}) where {T} = T
function write(parent::Union{JldFile, JldGroup}, name::String, c::Complex)
    reim = [real(c), imag(c)]
    write(parent, name, reim, full_typename(typeof(c)))
end
function write(parent::Union{JldFile, JldGroup}, name::String, C::Array{T}) where T<:Complex
    reim = reinterpret(realtype(T), reshape(C, ntuple(i->i==1 ? 1 : size(C,i-1), ndims(C)+1)))
    write(parent, name, reim, full_typename(typeof(C)))
end

# Int128/UInt128

# Symbols
write(parent::Union{JldFile, JldGroup}, name::String, sym::Symbol) = write(parent, name, string(sym), "Symbol")
write(parent::Union{JldFile, JldGroup}, name::String, syms::Array{Symbol}) = write(parent, name, map(string, syms), full_typename(typeof(syms)))

# Char
write(parent::Union{JldFile, JldGroup}, name::String, char::Char) = write(parent, name, UInt32(char), "Char")

# General array types (as arrays of references)
function write(parent::Union{JldFile, JldGroup}, path::String, data::Array{T}, astype::AbstractString) where T
    local gref  # a group, inside /_refs, for all the elements in data
    local refs
    # Determine whether parent already exists in /_refs, so we can avoid group/dataset conflict
    pname = name(parent)
    if startswith(pname, pathrefs)
        gref = create_group(parent, path*"g")
    else
        pathr = _joinpath(pathrefs, pname, path)
        if haskey(file(parent), pathr)
            gref = open_group(file(parent), pathr)
        else
            gref = create_group(file(parent), pathr)
        end
    end
    grefname = name(gref)
    try
        # Write the items to the reference group
        refs = Array{HDF5.Reference}(undef, size(data))
        # pad with zeros to keep in order
        nd = ndigits(length(data))
        z = "0"
        z = z[ones(Int, nd-1)]
        nd = 1
        for i = 1:length(data)
            if isassigned(data, i)
                if ndigits(i) > nd
                    nd = ndigits(i)
                    z = z[1:end-1]
                end
                itemname = z*string(i)
                write(gref, itemname, data[i])
                # Extract references
                tmp = gref[itemname]
                refs[i] = HDF5.Reference(tmp.plain, grefname*"/"*itemname)
                close(tmp)
            else
                refs[i] = HDF5.Reference()
            end
        end
    finally
        close(gref)
    end
    # Write the references as the chosen variable
    cset, ctype = create_dataset(parent.plain, path, refs)
    try
        HDF5.writearray(cset, ctype.id, refs)
        write_attribute(cset, name_type_attr, astype)
    finally
        close(ctype)
        close(cset)
    end
end
write(parent::Union{JldFile, JldGroup}, path::String, data::Array{T}) where {T} =
    write(parent, path, data, full_typename(typeof(data)))

# Tuple
write(parent::Union{JldFile, JldGroup}, name::String, t::Tuple) = write(parent, name, Any[t...], "Tuple")

# AbstractDict
function write(parent::Union{JldFile, JldGroup}, name::String, d::AbstractDict)
    tn = full_typename(typeof(d))
    if tn == "DataFrame"
        return write_composite(parent, name, d)
    end
    n = length(d)
    ks = Vector{keytype(d)}(undef, n)
    vs = Vector{valtype(d)}(undef, n)
    i = 0
    for (k,v) in d
        ks[i+=1] = k
        vs[i] = v
    end
    da = Any[ks, vs]
    write(parent, name, da, tn)
end

# Expressions
function write(parent::Union{JldFile, JldGroup}, name::String, ex::Expr)
    args = ex.args
    # Discard "line" expressions
    keep = trues(length(args))
    for i = 1:length(args)
        if isa(args[i], Expr) && args[i].head == :line
            keep[i] = false
        end
    end
    args = args[keep]
    a = Any[ex.head, args]
    write(parent, name, a, "Expr")
end

# CompositeKind
write(parent::Union{JldFile, JldGroup}, name::String, s; rootmodule="") = write_composite(parent, name, s; rootmodule=rootmodule)

function write_composite(parent::Union{JldFile, JldGroup}, name::String, s; rootmodule="")
    T = typeof(s)
    if isempty(fieldnames(T))
        if T.size > 0
            return write_bitstype(parent, name, s)
        end
        isdefined(T, :instance) || error("This is the write function for CompositeKind, but the input is of type ", T)
    end
    if has_pointer_field(s, name)
        return
    end
    Tname = string(T.name.name)
    n = fieldnames(T)
    local gtypes
    if !haskey(file(parent), pathtypes)
        gtypes = create_group(file(parent), pathtypes)
    else
        gtypes = parent[pathtypes]
    end
    try
        if !haskey(gtypes, Tname)
            # Write names to a dataset, so that other languages reading this file can
            # at least create a sensible dict
            nametype = Matrix{String}(2, length(n))
            t = T.types
            for i = 1:length(n)
                nametype[1, i] = string(n[i])
                nametype[2, i] = string(t[i])
            end
            write(gtypes.plain, Tname, nametype)
            obj = gtypes[Tname]
            # Write the module name as an attribute
            mod = Base.fullname(T.name.module)
            modnames = [map(string, mod)...]
            indx = findfirst(x->x==rootmodule, modnames)
            if indx > 0
                modnames = modnames[indx+1:end]
            end
            write_attribute(obj.plain, "Module", modnames)
            close(obj)
        end
    finally
        close(gtypes)
    end
    # Write the data
    v = Vector{Any}(undef, length(n))
    for i = 1:length(v)
        if isdefined(s, n[i])
            v[i] = getfield(s, n[i])
        end
    end
    write(parent, name, v, "CompositeKind")
    obj = parent[name]
    write_attribute(obj.plain, "CompositeKind", Tname)
    params = [map(full_typename, T.parameters)...]
    write_attribute(obj.plain, "TypeParameters", params)
    close(obj)
end

function write_bitstype(parent::Union{JldFile, JldGroup}, name::String, s)
    T = typeof(s)
    if T.size == 1
        ub = reinterpret(UInt8, s)
    elseif T.size == 2
        ub = reinterpret(UInt16, s)
    elseif T.size == 4
        ub = reinterpret(UInt32, s)
    elseif T.size == 8
        ub = reinterpret(UInt64, s)
    else
        error("Unsupported bitstype $T of size $(T.size)")
    end
    write(parent, name, [ub], "$(full_typename(T))")
end

function has_pointer_field(obj::Tuple, name)
    for o in obj
        if has_pointer_field(o, name)
            return true
        end
    end
    false
end

function has_pointer_field(obj, name)
    names = fieldnames(typeof(obj))
    for fieldname in names
        if isdefined(obj, fieldname)
            x = getfield(obj, fieldname)
            if isa(x, Ptr)
                @warn("Skipping $name because field \"$fieldname\" is a pointer")
                return true
            end
            if !isa(x, AbstractDict) && has_pointer_field(x, name)
                return true
            end
        end
    end
    false
end

### Size, length, etc ###
function size(dset::JldDataset)
    if !haskey(attributes(dset.plain), name_type_attr)
        return size(dset.plain)
    end
    # Read the type
    typename = read_attribute(dset.plain, name_type_attr)
    if typename == "Tuple"
        return size(dset.plain)
    end
    # Convert to Julia type
    T = julia_type(typename)
    if T == CompositeKind || T <: AbstractDict || T == Expr
        return ()
    elseif T <: Complex
        return ()
    elseif isarraycomplex(T)
        sz = size(dset.plain)
        return sz[2:end]
    end
    size(dset.plain)
end
length(dset::JldDataset) = prod(size(dset))
lastindex(dset::JldDataset) = length(dset)

isarraycomplex(::Type{Array{T, N}}) where {T<:Complex, N} = true
isarraycomplex(t) = false

### Read/write via getindex/setindex! ###
function getindex(dset::JldDataset, indices::Union{Integer, Base.RangeIndex}...)
    if !haskey(attributes(dset.plain), name_type_attr)
        # Fallback to plain read
        return getindex(dset.plain, indices...)
    end
    # Read the type
    typename = read_attribute(dset.plain, name_type_attr)
    if typename == "Tuple"
        return read_tuple(dset, indices...)
    end
    # Convert to Julia type
    T = julia_type(typename)
    _getindex(dset, T, indices...)
end

_getindex(dset::JldDataset, ::Type{Array{T,N}}, indices::Base.RangeIndex...) where {T<:HDF5.BitsType,N} =
    HDF5._getindex(dset.plain, T, indices...)
function _getindex(dset::JldDataset, ::Type{Array{T,N}}, indices::Base.RangeIndex...) where {T<:Complex,N}
    reinterpret(T, HDF5._getindex(dset.plain, realtype(T), 1:2, indices...), ntuple(i->length(indices[i]), length(indices)))
end
function _getindex(dset::JldDataset, ::Type{Array{Bool,N}}, indices::Base.RangeIndex...) where N
    tf = HDF5._getindex(dset.plain, UInt8, indices...)
    bool(tf)
end
_getindex(dset::JldDataset, ::Type{Array{T,N}}, indices::Union{Integer, Base.RangeIndex}...) where {T,N} =
    getrefs(dset, T, indices...)
function setindex!(dset::JldDataset, X::Array, indices::Base.RangeIndex...)
    if !haskey(attributes(dset.plain), name_type_attr)
        # Fallback to plain read
        return setindex!(dset.plain, X, indices...)
    end
    # Read the type
    typename = read_attribute(dset.plain, name_type_attr)
    if typename == "Tuple"
        return read_tuple(dset, indices...)
    end
    # Convert to Julia type
    T = julia_type(typename)
    HDF5._setindex!(dset, T, X, indices...)
end

length(x::Union{JldFile, JldGroup}) = length(keys(x))

### Dump ###
function dump(io::IO, parent::Union{JldFile, JldGroup}, n::Int, indent)
    nms = keys(parent)
    println(io, typeof(parent), " len ", length(nms))
    if n > 0
        i = 1
        for k in nms
            print(io, indent, "  ", k, ": ")
            v = parent[k]
            if isa(v, HDF5.Group)
                dump(io, v, n-1, string(indent, "  "))
            else
                if haskey(attributes(v.plain), name_type_attr)
                    typename = read_attribute(v.plain, name_type_attr)
                    if length(typename) >= 5 && (typename[1:5] == "Array" || typename[1:5] == "Tuple")
                        println(io, typename, " ", size(v))
                    else
                        println(io, typename)
                    end
                else
                    dump(io, v, 1, indent)
                end
            end
            close(v)
            if i > n
                println(io, indent, "  ...")
                break
            end
            i += 1
        end
    end
end



### Converting attribute strings to Julia types

const _typedict = Dict{String,Type}()
_typedict["CompositeKind"] = CompositeKind
function _julia_type(s::AbstractString)
    typ = get(_typedict, s, UnconvertedType)
    if typ == UnconvertedType
        e = Meta.parse(s)
        e = JLD.fixtypes(e)
        typ = UnsupportedType
        if JLD.is_valid_type_ex(e)
            try     # try needed to catch undefined symbols
                typ = eval(e)
                if !isa(typ, Type)
                    typ = UnsupportedType
                end
            catch
                try
                    typ = Core.eval(Main, e)
                catch
                    typ = UnsupportedType
                    if !isa(typ, Type)
                        typ = UnsupportedType
                    end
                end
            end
        else
            typ = UnsupportedType
        end
        if typ != UnsupportedType
            _typedict[s] = typ
        end
    end
    typ
end

### Converting Julia types to fully qualified names
full_typename(jltype::Union) = @sprintf "Union(%s)" join(map(full_typename, jltype.types), ",")
function full_typename(tv::TypeVar)
    if is(tv.lb, Union{}) && is(tv.ub, Any)
        "TypeVar(:$(tv.name))"
    elseif is(tv.lb, Union{})
        "TypeVar(:$(tv.name),$(full_typename(tv.ub)))"
    else
        "TypeVar(:$(tv.name),$(full_typename(tv.lb)),$(full_typename(tv.ub)))"
    end
end
full_typename(jltype::Tuple{Vararg{Type}}) =
    length(jltype) == 1 ? @sprintf("(%s,)", full_typename(jltype[1])) :
                          @sprintf("(%s)", join(map(full_typename, jltype), ","))
full_typename(x) = string(x)
function full_typename(jltype::DataType)
    #tname = "$(jltype.name.module).$(jltype.name)"
    tname = string(jltype.name.module, ".", jltype.name.name)  # NOTE: performance bottleneck
    if isempty(jltype.parameters)
        tname
    else
        @sprintf "%s{%s}" tname join([full_typename(x) for x in jltype.parameters], ",")
    end
end

### Version number utilities
versionnum(v::AbstractString) = map(int, split(v, '.'))
versionstring(v::Array{Int}) = join(v, '.')
function isversionless(l::Array{Int}, r::Array{Int})
    len = min(length(l), length(r))
    for i = 1:len
        if l[i] < r[i]
            return true
        end
    end
    if length(r) > len
        for i = len+1:length(r)
            if r[i] > 0
                return true
            end
        end
    end
    false
end

function Base.keys(parent::Union{JldFile, JldGroup})
    n = keys(parent.plain)
    keep = trues(length(n))
    reserved = [pathrefs[2:end], pathtypes[2:end], pathrequire[2:end]]
    for i = 1:length(n)
        if in(n[i], reserved)
            keep[i] = false
        end
    end
    n[keep]
end

macro save(filename, vars...)
    if isempty(vars)
        # Save all variables in the current module
        writeexprs = Vector{Expr}(undef, 0)
        m = @__MODULE__
        for vname in names(m)
            s = string(vname)
            if !ismatch(r"^_+[0-9]*$", s) # skip IJulia history vars
                v = Core.eval(m, vname)
                if !isa(v, Module)
                    push!(writeexprs, :(if !isa($(esc(vname)), Function) write(f, $s, $(esc(vname))) end))
                end
            end
        end
    else
        writeexprs = Vector{Expr}(undef, length(vars))
        for i = 1:length(vars)
            writeexprs[i] = :(write(f, $(string(vars[i])), $(esc(vars[i]))))
        end
    end
    Expr(:block,
         :(local f = jldopen($(esc(filename)), "w")),
         Expr(:try, Expr(:block, writeexprs...), false, false,
              :(close(f))))
end

macro load(filename, vars...)
    if isempty(vars)
        if isa(filename, Expr)
            filename = Core.eval(@__MODULE__, filename)
        end
        # Load all variables in the top level of the file
        readexprs = Vector{Expr}(undef, 0)
        vars = Vector{Expr}(undef, 0)
        f = jldopen(filename)
        nms = keys(f)
        for n in nms
            obj = f[n]
            if isa(obj, JldDataset)
                sym = esc(Symbol(n))
                push!(readexprs, :($sym = read($f, $n)))
                push!(vars, sym)
            end
        end
        return Expr(:block,
                    Expr(:global, vars...),
                    Expr(:try,  Expr(:block, readexprs...), false, false,
                         :(close($f))),
                    Symbol[v.args[1] for v in vars]) # "unescape" vars
    else
        readexprs = Vector{Expr}(undef, length(vars))
        for i = 1:length(vars)
            readexprs[i] = :($(esc(vars[i])) = read(f, $(string(vars[i]))))
        end
        return Expr(:block,
                    Expr(:global, map(esc, vars)...),
                    :(local f = jldopen($(esc(filename)))),
                    Expr(:try,  Expr(:block, readexprs...), false, false,
                         :(close(f))),
                    Symbol[v for v in vars]) # vars is a tuple
    end
end

# Save all the key-value pairs in the dict as top-level variables of the JLD
function save(filename::AbstractString, dict::AbstractDict)
    jldopen(filename, "w") do file
        for (k,v) in dict
            write(file, String(k), v)
        end
    end
end
# Or the names and values may be specified as alternating pairs
function save(filename::AbstractString, name::AbstractString, value, pairs...)
    if isodd(length(pairs)) || !isa(pairs[1:2:end], (AbstractString...))
        throw(ArgumentError("arguments must be in name-value pairs"))
    end
    jldopen(filename, "w") do file
        write(file, String(name), value)
        for i = 1:2:length(pairs)
            write(file, String(pairs[i]), pairs[i+1])
        end
    end
end

# load with just a filename returns a dictionary containing all the variables
function load(filename::AbstractString)
    jldopen(filename, "r") do file
        Dict{String,Any}([(var, read(file, var)) for var in keys(file)])
    end
end
# When called with explicitly requested variable names, return each one
function load(filename::AbstractString, varname::AbstractString)
    jldopen(filename, "r") do file
        read(file, varname)
    end
end
load(filename::AbstractString, varnames::AbstractString...) = load(filename, varnames)
function load(filename::AbstractString, varnames::Tuple{Vararg{AbstractString}})
    jldopen(filename, "r") do file
        map((var)->read(file, var), varnames)
    end
end

function addrequire(file::JldFile, filename::AbstractString)
    files = read(file, pathrequire)
    push!(files, filename)
    delete_object(file, pathrequire)
    write(file, pathrequire, files)
end

export
    addrequire,
    ismmappable,
    jldopen,
    delete_object,
    readmmap,
    readsafely,
    @load,
    @save,
    load,
    save

###
### v0.12.0 deprecations
###

@noinline function Base.names(parent::Union{JldFile, JldGroup})
    Base.depwarn("`names(parent)` is deprecated, use `keys(parent)` instead.", :names)
    return keys(parent)
end

end
