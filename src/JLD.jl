module JLD
using Printf
using HDF5, H5Zblosc, FileIO
using Compat

import HDF5: file, create_group, open_group, delete_object, name, ismmappable, readmmap
import Base: close, convert, datatype_pointerfree, delete!, dump, eltype, getindex, iterate,
             length, ndims, read, setindex!, show, size, sizeof, unsafe_convert, write

@noinline gcuse(x) = x # because of use of `pointer`, need to mark gc-use end explicitly

const magic_base = "Julia data file (HDF5), version "
const version_current = v"0.1.3"
const pathrefs = "/_refs"
const pathtypes = "/_types"
const pathrequire = "/_require"
const pathcreator = "/_creator"
const name_type_attr = "julia type"

# HDF5 paths are Unix-like, even on Windows
_joinpath(path::AbstractString)::String = path
function _joinpath(path::AbstractString, paths::AbstractString...)::String
    for p in paths
        if startswith(p, '/')
            path = p
        elseif isempty(path) || path[end] == '/'
            path *= p
        else
            path *= "/" * p
        end
    end
    return path
end

const BitsKindOrString = Union{HDF5.BitsType, String}

function julia_type(s::AbstractString)
    s = replace(s, r"ASCIIString|UTF8String|ByteString" => "String")
    if occursin("Base.UTF16String", s)
        error("file contains Base.UTF16String, must be converted and re-saved with JLD 0.9 or less")
    end
    _julia_type(s)
end

### Dummy types used for converting attribute strings to Julia types
mutable struct UnsupportedType; end
mutable struct UnconvertedType; end


struct JldDatatype
    dtype::HDF5.Datatype
    index::Int
end
sizeof(T::JldDatatype) = sizeof(T.dtype)

struct JldWriteSession
    persist::Vector{Any} # To hold objects that should not be garbage-collected
    h5ref::Base.IdDict{Any,Any}        # To hold mapping from Object/Array -> HDF5.Reference

    JldWriteSession() = new(Any[], Base.IdDict{Any,Any}())
end

# The Julia Data file type
# Purpose of the nrefs field:
# length(group) only returns the number of _completed_ items in a group. Since
# we'll write recursively, we need to keep track of the number of reference
# objects _started_.
mutable struct JldFile <: HDF5.H5DataStore
    plain::HDF5.File
    version::VersionNumber
    toclose::Bool
    writeheader::Bool
    mmaparrays::Bool
    compatible::Bool
    compress::Bool
    h5jltype::Dict{Int,Type}
    jlh5type::Dict{Type,JldDatatype}
    jlref::Dict{HDF5.Reference,WeakRef}
    truncatemodules::Vector{String}
    gref # Group references; can't annotate type here due to circularity
    nrefs::Int

    function JldFile(plain::HDF5.File, version::VersionNumber=version_current, toclose::Bool=true,
                     writeheader::Bool=false, mmaparrays::Bool=false,
                     compatible::Bool=false, compress::Bool=false)
        f = new(plain, version, toclose, writeheader, mmaparrays, compatible,
                compress & !mmaparrays,
                Dict{HDF5.Datatype,Type}(), Dict{Type,HDF5.Datatype}(),
                Dict{HDF5.Reference,WeakRef}(), String[])
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

iscompatible(f::JldFile) = f.compatible
iscompatible(g::JldGroup) = g.file.compatible

iscompressed(f::JldFile) = f.compress
iscompressed(g::JldGroup) = g.file.compress

ismmapped(f::JldFile) = f.mmaparrays
ismmapped(d::JldGroup) = d.file.mmaparrays

struct PointerException <: Exception; end
show(io::IO, ::PointerException) = print(io, "cannot write a pointer to JLD file")

struct TypeMismatchException <: Exception
    typename::String
end
show(io::IO, e::TypeMismatchException) =
    print(io, "stored type $(e.typename) does not match currently loaded type")

include("jld_types.jl")

file(x::JldFile) = x
file(x::Union{JldGroup, JldDataset}) = x.file

function close(f::JldFile)
    if f.toclose
        # Close types
        for x in values(f.jlh5type)
            close(x.dtype)
        end

        # Close reference group
        isdefined(f, :gref) && close(f.gref)

        # Ensure that all other datasets, groups, and datatypes are closed (ref #176)
        for obj_id in HDF5.API.h5f_get_obj_ids(f.plain.id, HDF5.API.H5F_OBJ_DATASET | HDF5.API.H5F_OBJ_GROUP | HDF5.API.H5F_OBJ_DATATYPE)
            HDF5.API.h5o_close(obj_id)
        end

        # Close file
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

function jldopen(filename::AbstractString, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool; mmaparrays::Bool=false, compatible::Bool=false, compress::Bool=false)
    local fj
    if ff && !wr
        error("Cannot append to a write-only file")
    end
    if !cr && !isfile(filename)
        error("File ", filename, " cannot be found")
    end
    version = version_current
    fapl = HDF5.FileAccessProperties()
    # HDF5.h5p_set_libver_bounds(fapl, HDF5.H5F_LIBVER_18, HDF5.H5F_LIBVER_18)
    try
        fapl.fclose_degree = :strong
        if cr && (tr || !isfile(filename))
            # We're truncating, so we don't have to check the format of an existing file
            # Set the user block to 512 bytes, to save room for the header
            fcpl = HDF5.FileCreateProperties()
            local f
            try
                fcpl.userblock = 512
                f = HDF5.API.h5f_create(filename, HDF5.API.H5F_ACC_TRUNC, fcpl, fapl)
            finally
                close(fcpl)
            end
            fj = JldFile(HDF5.File(f, filename, false), version, true, true, mmaparrays, compatible, compress)
            # Record creator information. Don't use any fancy types here,
            # because we want to be able to read this even when formats change.
            write(fj, _joinpath(pathcreator, "JULIA_MAJOR"), VERSION.major)
            write(fj, _joinpath(pathcreator, "JULIA_MINOR"), VERSION.minor)
            write(fj, _joinpath(pathcreator, "JULIA_PATCH"), VERSION.patch)
            if !isempty(VERSION.prerelease)
                write(fj, _joinpath(pathcreator, "JULIA_PRERELEASE"), [VERSION.prerelease...])
            end
            if !isempty(VERSION.build)
                write(fj, _joinpath(pathcreator, "JULIA_BUILD"), [VERSION.build...])
            end
            write(fj, _joinpath(pathcreator, "WORD_SIZE"), Sys.WORD_SIZE)
            write(fj, _joinpath(pathcreator, "ENDIAN_BOM"), ENDIAN_BOM)
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
                version = VersionNumber(unsafe_string(pointer(magic) + length(magic_base)))
                gcuse(magic)
                if version < v"0.1.0"
                    fj = JLD00.jldopen(filename, rd, wr, cr, tr, ff; mmaparrays=mmaparrays)
                else
                    f = HDF5.API.h5f_open(filename, wr ? HDF5.API.H5F_ACC_RDWR : HDF5.API.H5F_ACC_RDONLY, fapl)
                    fj = JldFile(HDF5.File(f, filename, false), version, true, cr|wr, mmaparrays, compatible, compress)
                    # Load any required files/packages
                    if haskey(fj, pathrequire)
                        r = read(fj, pathrequire)
                        for fn in r
                            mod = path2modsym(fn)
                            Core.eval(Main, Expr(:import, Expr(:., mod)))
                        end
                    end
                end
            else
                if HDF5.ishdf5(filename)
                    println("$filename is an HDF5 file, but it is not a recognized Julia data file. Opening anyway.")
                    # inverse of logic in jldopen() below
                    mode =  rd & !wr & !cr & !tr & !ff ? "r" :
                            rd &  wr & !cr & !tr & !ff ? "r+" :
                           !rd &  wr &  cr &  tr & !ff ? "w" : error("invalid mode")
                    fj = JldFile(h5open(filename, mode), version_current, true, false, mmaparrays, compatible, compress)
                else
                    error("$filename does not seem to be a Julia data or HDF5 file")
                end
            end
        end
    finally
        close(fapl)
    end
    return fj
end

function jldopen(fname::AbstractString, mode::AbstractString="r"; mmaparrays::Bool=false, compatible::Bool=false, compress::Bool=false)
    mode == "r"  ? jldopen(fname, true , false, false, false, false, mmaparrays=mmaparrays, compatible=compatible, compress=compress) :
    mode == "r+" ? jldopen(fname, true , true , false, false, false, mmaparrays=mmaparrays, compatible=compatible, compress=compress) :
    mode == "w"  ? jldopen(fname, false, true , true , true , false, mmaparrays=mmaparrays, compatible=compatible, compress=compress) :
#     mode == "w+" ? jldopen(fname, true , true , true , true , false) :
#     mode == "a"  ? jldopen(fname, false, true , true , false, true ) :
#     mode == "a+" ? jldopen(fname, true , true , true , false, true ) :
    error("invalid open mode: ", mode)
end

function jldopen(f::Function, args...; kws...)
    jld = jldopen(args...; kws...)
    try
        f(jld)
    finally
        close(jld)
    end
end

function creator(file::JldFile, key::AbstractString)
    if file.version < v"0.1.1"
        return nothing
    end
    if key == "VERSION"
        path_prerelease = _joinpath(pathcreator, "JULIA_PRERELEASE")
        prerelease = haskey(file, path_prerelease) ? tuple(read(file, path_prerelease)...) : ()
        path_build = _joinpath(pathcreator, "JULIA_BUILD")
        build = haskey(file, path_build) ? tuple(read(file, path_build)...) : ()
        return VersionNumber(read(file, _joinpath(pathcreator, "JULIA_MAJOR")),
                             read(file, _joinpath(pathcreator, "JULIA_MINOR")),
                             read(file, _joinpath(pathcreator, "JULIA_PATCH")),
                             prerelease,
                             build)
    elseif key in ("WORD_SIZE", "ENDIAN_BOM")
        return read(file, _joinpath(pathcreator, key))
    else
        error("$key not recognized")
    end
end

function jldobject(obj_id::HDF5.API.hid_t, parent)
    obj_type = HDF5.API.h5i_get_type(obj_id)
    obj_type == HDF5.API.H5I_GROUP ? JldGroup(HDF5.Group(obj_id, file(parent.plain)), file(parent)) :
    obj_type == HDF5.API.H5I_DATATYPE ? HDF5.Datatype(obj_id) :
    obj_type == HDF5.API.H5I_DATASET ? JldDataset(HDF5.Dataset(obj_id, file(parent.plain)), file(parent)) :
    error("Invalid object type for path ", path)
end

getindex(parent::Union{JldFile, JldGroup}, path::String) =
    jldobject(HDF5.API.h5o_open(parent.plain.id, path, HDF5.API.H5P_DEFAULT), parent)

function getindex(parent::Union{JldFile, JldGroup, JldDataset}, r::HDF5.Reference)
    r == HDF5.Reference() && error("Reference is null")
    obj_id = HDF5.API.h5r_dereference(parent.plain.id, HDF5.API.H5P_DEFAULT, HDF5.API.H5R_OBJECT, r)
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
eltype(p::JldDataset) = eltype(p.plain)
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
    readas(val)
end
read(parent::Union{JldFile,JldGroup}, name::Symbol) = read(parent, String(string(name)))

function read(obj::JldGroup)
    nms = keys(obj)
    val = Dict{String, Any}()
    for nm in nms
        val[nm] = read(obj[nm])
    end
    return val
end

function read(obj::JldDataset)
    dtype = datatype(obj.plain)
    dspace_id = HDF5.API.h5d_get_space(obj.plain)
    extent_type = HDF5.API.h5s_get_simple_extent_type(dspace_id)
    try
        if extent_type == HDF5.API.H5S_SCALAR
            # Scalar value
            return read_scalar(obj, dtype, jldatatype(file(obj), dtype))
        elseif extent_type == HDF5.API.H5S_SIMPLE
            return read_array(obj, dtype, dspace_id, HDF5.API.H5S_ALL)
        elseif extent_type == HDF5.API.H5S_NULL
            # Empty array
            if HDF5.API.h5t_get_class(dtype) == HDF5.API.H5T_REFERENCE
                T = refarray_eltype(obj)
            else
                T = jldatatype(file(obj), dtype)
            end
            if haskey(obj, "dims")
                dims = read_attribute(obj.plain, "dims")
                return Array{T}(undef, dims...)
            else
                return T[]
            end
        end
    finally
        HDF5.API.h5s_close(dspace_id)
    end
end

## Scalars
function read_scalar(obj::JldDataset, dtype::HDF5.Datatype, ::Type{T}) where {T<:BitsKindOrString}
    if T === Bool && obj.file.version < v"0.1.3"
        return read_scalar_default(obj, dtype, Bool)
    end
    return read(obj.plain, T)
end
read_scalar(obj::JldDataset, dtype::HDF5.Datatype, T::Type) = read_scalar_default(obj, dtype, T)
function read_scalar_default(obj::JldDataset, dtype::HDF5.Datatype, T::Type)
    buf = Vector{UInt8}(undef, sizeof(dtype))
    read_dataset(obj.plain, dtype, buf)
    sc = readas(jlconvert(T, file(obj), pointer(buf)))
    gcuse(buf)
    sc
end

## Arrays

# Read an array
function read_array(obj::JldDataset, dtype::HDF5.Datatype, dspace_id::HDF5.API.hid_t, dsel_id::HDF5.API.hid_t,
                    dims::Tuple{Vararg{Int}} = (Int.(reverse!(HDF5.API.h5s_get_simple_extent_dims(dspace_id)[1]))...,))
    if HDF5.API.h5t_get_class(dtype) == HDF5.API.H5T_REFERENCE
        val = read_refs(obj, refarray_eltype(obj), dspace_id, dsel_id, dims)
    else
        val = read_vals(obj, dtype, jldatatype(file(obj), dtype), dspace_id, dsel_id, dims)
    end
    readas(val)
end

# Arrays of basic HDF5 kinds
function read_vals(obj::JldDataset, dtype::HDF5.Datatype, T::Union{Type{S}, Type{Complex{S}}},
                   dspace_id::HDF5.API.hid_t, dsel_id::HDF5.API.hid_t, dims::Tuple{Vararg{Int}}) where {S<:HDF5.BitsType}
    if S === Bool && obj.file.version < v"0.1.3"
        return read_vals_default(obj, dtype, T, dspace_id, dsel_id, dims)
    end
    if obj.file.mmaparrays && HDF5.iscontiguous(obj.plain) && dsel_id == HDF5.API.H5S_ALL
        readmmap(obj.plain, T)
    else
        out = Array{T}(undef, dims)
        HDF5.API.h5d_read(obj.plain.id, dtype.id, dspace_id, dsel_id, HDF5.API.H5P_DEFAULT, out)
        out
    end
end

# Arrays of immutables/bitstypes
function read_vals(obj::JldDataset, dtype::HDF5.Datatype, T::Type, dspace_id::HDF5.API.hid_t,
                   dsel_id::HDF5.API.hid_t, dims::Tuple{Vararg{Int}})
    return read_vals_default(obj, dtype, T, dspace_id, dsel_id, dims)
end

function read_vals_default(obj::JldDataset, dtype::HDF5.Datatype, T::Type, dspace_id::HDF5.API.hid_t,
                           dsel_id::HDF5.API.hid_t, dims::Tuple{Vararg{Int}})
    out = Array{T}(undef, dims)
    # Empty objects don't need to be read at all
    !ismutabletype(T) && sizeof(T) == 0 && return out

    # Read from file
    n = prod(dims)
    h5sz = sizeof(dtype)
    buf = Vector{UInt8}(undef, h5sz*n)
    HDF5.API.h5d_read(obj.plain.id, dtype.id, dspace_id, dsel_id, HDF5.API.H5P_DEFAULT, buf)

    f = file(obj)
    h5offset = pointer(buf)
    if datatype_pointerfree(T) && !ismutabletype(T)
        jloffset = pointer(out)
        jlsz = sizeof(T)

        # Perform conversion in buffer
        for i = 1:n
            jlconvert!(jloffset, T, f, h5offset)
            jloffset += jlsz
            h5offset += h5sz
        end
    else
        # Convert each item individually
        for i = 1:n
            out[i] = jlconvert(T, f, h5offset)
            h5offset += h5sz
        end
    end
    gcuse(buf)
    gcuse(out)
    out
end

# Arrays of references
function read_refs(obj::JldDataset, ::Type{T}, dspace_id::HDF5.API.hid_t, dsel_id::HDF5.API.hid_t,
                   dims::Tuple{Vararg{Int}}) where T
    refs = Array{HDF5.Reference}(undef, dims)
    HDF5.API.h5d_read(obj.plain.id, HDF5.API.H5T_STD_REF_OBJ, dspace_id, dsel_id, HDF5.API.H5P_DEFAULT, refs)

    out = Array{T}(undef, dims)
    f = file(obj)
    for i = 1:length(refs)
        if refs[i] != HDF5.Reference()
            out[i] = read_ref(f, refs[i])
        end
    end
    out
end

# Get element type of a reference array
function refarray_eltype(obj::JldDataset)
    typename = read_attribute(obj.plain, "julia eltype")
    T = julia_type(typename)
    if T == UnsupportedType
        @warn("type $typename not present in workspace; interpreting array as Array{Any}")
        return Any
    end
    return T
end

## Reference
function read_ref(f::JldFile, ref::HDF5.Reference)
    if haskey(f.jlref, ref)
        # Stored as WeakRefs and may no longer exist
        val = f.jlref[ref].value
        val !== nothing && return val
    end

    dset = f[ref]
    data = try
        read(dset)
    finally
        close(dset)
    end

    f.jlref[ref] = WeakRef(data)
    readas(data)
end

### Writing ###

write(parent::Union{JldFile, JldGroup}, name::String,
      data, wsession::JldWriteSession=JldWriteSession(); kargs...) =
    close(_write(parent, name, writeas(data), wsession; kargs...))

# Pick whether to use compact or default storage based on data size
function dset_create_properties(parent, sz::Int, obj, chunk=Int[]; mmap::Bool=false)
    if sz <= 8192 && !ismmapped(parent) && !mmap
        return compact_properties(), false
    end
    if iscompressed(parent) && !isempty(chunk)
        p = HDF5.DatasetCreateProperties()
        p.chunk = chunk
        if iscompatible(parent)
            p.shuffle = true
            p.deflate = 5
        else
            p.blosc = 5
        end
        return p, true
    else
        return HDF5.DatasetCreateProperties(), false
    end
end

# Write "basic" types
function _write(parent::Union{JldFile, JldGroup},
                    name::String,
                    data::Union{T, Array{T}},
                    wsession::JldWriteSession; kargs...) where T<:Union{HDF5.BitsType, String}
    chunk = T <: String ? Int[] : HDF5.heuristic_chunk(data)
    dprop, dprop_close = dset_create_properties(parent, sizeof(data), data, chunk; kargs...)
    dtype = datatype(data)
    dset = HDF5.Dataset(HDF5.API.h5d_create(parent.plain, String(name), dtype, dataspace(data),
                                        HDF5._link_properties(name), dprop,
                                        HDF5.DatasetAccessProperties()), file(parent.plain))
    try
        # Write the attribute
        isa(data, Array) && isempty(data) && write_attribute(dset, "dims", [size(data)...])
        # Write the data
        write_dataset(dset, dtype, data)
    finally
        close(dtype)
        dprop_close && close(dprop)
    end
    dset
end

# General array types
function _write(parent::Union{JldFile, JldGroup},
                path::String, data::Array{T},
                wsession::JldWriteSession; kargs...) where T
    f = file(parent)
    dtype = h5fieldtype(f, T, true)
    buf = h5convert_array(f, data, dtype, wsession)
    dims = convert(Vector{HDF5.API.hsize_t}, [reverse(size(data))...])
    dspace = dataspace(data)
    chunk = HDF5.heuristic_chunk(dtype, size(data))
    dprop, dprop_close = dset_create_properties(parent, sizeof(buf),buf, chunk; kargs...)
    try
        dset = HDF5.Dataset(HDF5.API.h5d_create(parent.plain, path, dtype.dtype, dspace,
                                            HDF5._link_properties(path), dprop,
                                            HDF5.DatasetAccessProperties()), file(parent.plain))
        if dtype == JLD_REF_TYPE
            write_attribute(dset, "julia eltype", full_typename(f, T))
        end
        if isempty(data) && ndims(data) != 1
            write_attribute(dset, "dims", [size(data)...])
        else
            write_dataset(dset, dtype.dtype, buf)
        end
        return dset
    finally
        dprop_close && close(dprop)
        close(dspace)
    end
end

# Dispatch correct method for Array{Union{}}
_write(parent::Union{JldFile, JldGroup}, path::String, data::Array{Union{}},
       wsession::JldWriteSession; kargs...) =
    # Keyword arguments do not currently work
    invoke(_write, Tuple{Union{JldFile, JldGroup},String,Array,JldWriteSession}, parent,
           path, data, wsession; kargs...)

# Convert an array to the format to be written to the HDF5 file, either
# references or values
function h5convert_array(f::JldFile, data::Array,
                         dtype::JldDatatype, wsession::JldWriteSession)
    if dtype == JLD_REF_TYPE
        # For type stability, return as Vector{UInt8}
        refs = Vector{UInt8}(undef, length(data)*sizeof(HDF5.Reference))
        arefs = reinterpret(HDF5.Reference, refs)
        for i = 1:length(data)
            if isassigned(data, i)
                arefs[i] = write_ref(f, data[i], wsession)
            else
                arefs[i] = HDF5.Reference()
            end
        end
        refs
    else
        gen_h5convert(f, eltype(data))
        h5convert_vals(f, data, dtype, wsession)
    end
end

# Hack to ensure that _h5convert_vals isn't compiled before h5convert!
function h5convert_vals(f::JldFile, @nospecialize(data), dtype::JldDatatype,
                        wsession::JldWriteSession)
    _h5convert_vals(f, data, dtype, wsession)
end

# Convert an array of immutables or bitstypes to a buffer representing
# HDF5 compound objects. A separate function so that it is specialized.
@noinline function _h5convert_vals(f::JldFile, data::Array,
                         dtype::JldDatatype, wsession::JldWriteSession)
    sz = HDF5.API.h5t_get_size(dtype)
    n = length(data)
    buf = Vector{UInt8}(undef, sz*n)
    offset = pointer(buf)
    for i = 1:n
        h5convert!(offset, f, data[i], wsession)
        offset += sz
    end
    gcuse(buf)
    buf
end

# Get reference group, creating a new one if necessary
function get_gref(f::JldFile)
    isdefined(f, :gref) && return f.gref::JldGroup

    if !haskey(f, pathrefs)
        gref = f.gref = create_group(f, pathrefs)
    else
        gref = f.gref = f[pathrefs]
    end
    f.nrefs = length(gref)
    gref
end

# Write a reference
function write_ref(parent::JldFile, data, wsession::JldWriteSession)
    REF_NULL = HDF5.Reference()
    # Check whether we have already written this object
    ref = get(wsession.h5ref, data, REF_NULL)
    ref != REF_NULL && return ref

    # Write an new reference
    gref = get_gref(parent)
    name = @sprintf "%08d" (parent.nrefs += 1)
    dset = _write(gref, name, writeas(data), wsession)

    # Add reference to reference list
    ref = HDF5.Reference(HDF5.API.hobj_ref_t(HDF5.API.h5o_get_info1(HDF5.checkvalid(dset)).addr))
    close(dset)
    if !isa(data, Tuple) && ismutable(data)
        wsession.h5ref[data] = ref
    end
    ref
end
write_ref(parent::JldGroup, data, wsession::JldWriteSession) =
    write_ref(file(parent), data, wsession)

# Expressions, drop line numbers
function _write(parent::Union{JldFile, JldGroup},
                name::String, ex::Expr,
                wsession::JldWriteSession; kargs...)
    args = ex.args
    # Discard "line" expressions
    keep = trues(length(args))
    for i = 1:length(args)
        if (isa(args[i], Expr) && args[i].head == :line) || isa(args[i], LineNumberNode)
            keep[i] = false
        end
    end
    newex = Expr(ex.head)
    newex.args = args[keep]
    write_compound(parent, name, newex, wsession)
end

# Generic (tuples, immutables, and compound types)
_write(parent::Union{JldFile, JldGroup}, name::String, s,
      wsession::JldWriteSession; kargs...) =
    write_compound(parent, name, s, wsession)
function write_compound(parent::Union{JldFile, JldGroup}, name::String,
                        s, wsession::JldWriteSession; kargs...)
    T = typeof(s)
    f = file(parent)
    dtype = h5type(f, T, true)
    gen_h5convert(f, T)

    buf = Vector{UInt8}(undef, HDF5.API.h5t_get_size(dtype))
    h5convert!(pointer(buf), file(parent), s, wsession)
    gcuse(buf)

    dspace = HDF5.Dataspace(HDF5.API.h5s_create(HDF5.API.H5S_SCALAR))
    dprop, dprop_close = dset_create_properties(parent, length(buf), buf; kargs...)
    try
        dset = HDF5.Dataset(HDF5.API.h5d_create(parent.plain, name, dtype.dtype, dspace,
                                            HDF5._link_properties(name), dprop,
                                            HDF5.DatasetAccessProperties()), file(parent.plain))
        write_dataset(dset, dtype.dtype, buf)
        return dset
    finally
        dprop_close && close(dprop)
        close(dspace)
    end
end

### Size, length, etc ###
size(dset::JldDataset) = size(dset.plain)
size(dset::JldDataset, d) = size(dset.plain, d)
length(dset::JldDataset) = prod(size(dset))
lastindex(dset::JldDataset) = length(dset)
ndims(dset::JldDataset) = ndims(dset.plain)

### Read/write via getindex/setindex! ###
function getindex(dset::JldDataset, indices::Union{AbstractRange{Int},Integer}...)
    sz = map(length, indices)
    dsel = HDF5.hyperslab(dset.plain, indices...)
    try
        dspace = HDF5._dataspace(sz)
        try
            return read_array(dset, datatype(dset.plain), dspace.id, dsel.id, sz)
        finally
            close(dspace)
        end
    finally
        HDF5.close(dsel)
    end
end

function setindex!(dset::JldDataset, X::AbstractArray{T,N}, indices::Union{AbstractRange{Int},Integer}...) where {T,N}
    f = file(dset)
    sz = map(length, indices)
    dsel_id = HDF5.hyperslab(dset.plain, indices...)
    try
        dtype = datatype(dset.plain)
        try
            # Convert array to writeable buffer
            if HDF5.API.h5t_get_class(dtype) == HDF5.API.H5T_REFERENCE
                written_eltype = refarray_eltype(dset)
                jldtype = JLD_REF_TYPE
            else
                written_eltype = jldatatype(f, dtype)
                jldtype = JldDatatype(dtype, -1)
            end

            buf = h5convert_array(f, convert(Array{written_eltype,N}, X), jldtype,
                                  JldWriteSession())

            dspace = HDF5._dataspace(sz)
            try
                HDF5.API.h5d_write(dset.plain.id, dtype, dspace, dsel_id, HDF5.API.H5P_DEFAULT, buf)
            finally
                close(dspace)
            end
        finally
            close(dtype)
        end
    finally
        HDF5.API.h5s_close(dsel_id)
    end
end
function setindex!(dset::JldDataset, x::Number, indices::Union{AbstractRange{Int},Integer}...)
    setindex!(dset, fill(x, map(length, indices)), indices...)
end

getindex(dset::JldDataset, I::Union{AbstractRange{Int},Integer,Colon}...) = getindex(dset, ntuple(i-> isa(I[i], Colon) ? (1:size(dset,i)) : I[i], length(I))...)
setindex!(dset::JldDataset, x, I::Union{AbstractRange{Int},Integer,Colon}...) = setindex!(dset, x, ntuple(i-> isa(I[i], Colon) ? (1:size(dset,i)) : I[i], length(I))...)

length(x::Union{JldFile, JldGroup}) = length(keys(x))

### Custom serialization

readas(x) = x
writeas(x) = x

# Wrapper for associative keys
# We write this instead of the associative to avoid dependence on the
# Julia hash function
struct AssociativeWrapper{K,V,T<:AbstractDict}
    keys::Vector{K}
    values::Vector{V}
end

readas(x::AssociativeWrapper{K,V,T}) where {K,V,T} = convert(T, x)
function writeas(x::T) where T<:AbstractDict
    K, V = destructure(eltype(x))
    convert(AssociativeWrapper{K,V,T}, x)
end
destructure(::Type{Pair{K,V}}) where {K,V} = K, V  # not inferrable, julia#10880

# Special case for associative, to rehash keys
function convert(::Type{T}, x::AssociativeWrapper{K,V,T}) where {K,V,T<:AbstractDict}
    ret = T()
    keys = x.keys
    values = x.values
    n = length(keys)
    if applicable(sizehint!, ret, n)
        sizehint!(ret, n)
    end
    for i = 1:n
        ret[keys[i]] = values[i]
    end
    ret
end

function convert(::Type{AssociativeWrapper{K,V,T}}, d::AbstractDict) where {K,V,T}
    n = length(d)
    ks = Vector{K}(undef, n)
    vs = Vector{V}(undef, n)
    i = 0
    for (k,v) in d
        ks[i+=1] = k
        vs[i] = v
    end
    AssociativeWrapper{K,V,typeof(d)}(ks, vs)
end

# Special case for SimpleVector
# Wrapper for SimpleVector
struct SimpleVectorWrapper
    elements::Vector
end

# Special case for SimpleVector
readas(x::SimpleVectorWrapper) = Core.svec(x.elements...)
writeas(x::Core.SimpleVector) = SimpleVectorWrapper([x...])

# function to convert string(mod::Module) back to mod::Module
function modname2mod(modname::AbstractString)
    Meta.parse(modname == "Main" ? modname : string("Main.", modname))
end


# Serializer for anonymous functions
# convert functions to lowered ast expressions
function func2expr(fun::Function)
    @assert !isgeneric(fun) "generic functions not supported"
    ast = Base.uncompressed_ast(fun.code)
    Expr(:function, Expr(:tuple, ast.args[1]...), Expr(:block, ast.args[3].args...))
end
struct AnonymousFunctionSerializer
    expr::Expr
    mod::AbstractString
    AnonymousFunctionSerializer(fun::Function) = new(func2expr(fun), string(fun.code.module))
end
readas(ast::AnonymousFunctionSerializer) = eval(modname2mod(ast.mod)).eval(ast.expr)
writeas(fun::Function) = AnonymousFunctionSerializer(fun)

# Serializer for GlobalRef
struct GlobalRefSerializer
    mod::AbstractString
    name::Symbol
    GlobalRefSerializer(g::GlobalRef) = new(string(g.mod), g.name)
end
readas(grs::GlobalRefSerializer) = GlobalRef(eval(modname2mod(grs.mod)), grs.name)
writeas(gr::GlobalRef) = GlobalRefSerializer(gr)

JLD.writeas(data::Base.StackTraces.StackFrame) =
    Base.StackTraces.StackFrame(data.func,
                    data.file,
                    data.line,
                    nothing,
                    data.from_c,
                    data.inlined,
                    data.pointer)

### Converting attribute strings to Julia types

const _where_macrocall = Symbol("@where")
function expand_where_macro(e::Expr)
    e.head = :where
    popfirst!(e.args)
    popfirst!(e.args)  # source location
    return true
end

is_valid_type_ex(s::Symbol) = true
is_valid_type_ex(s::QuoteNode) = true
is_valid_type_ex(s) = isbitstype(typeof(s))
function is_valid_type_ex(e::Expr)
    if e.head === :curly || e.head == :tuple || e.head == :.
        return all(is_valid_type_ex, e.args)
    elseif e.head === :where
        return is_valid_type_ex(e.args[1])
    elseif e.head === :let && length(e.args) == 2
        return is_valid_type_ex(e.args[2]) &&
               is_valid_type_ex(e.args[1].args[2])
    elseif e.head == :call
        f = e.args[1]
        if f isa Expr
            if f.head === :core
                f = f.args[1]
                return f === :Union || f === :TypeVar || f === :UnionAll
            end
        elseif f isa Symbol
            return f === :Union || f === :TypeVar || f === :symbol
        end
    end
    return false
end

const typemap_Core = Dict(
    :Uint8 => :UInt8,
    :Uint16 => :Uint16,
    :Uint32 => :UInt32,
    :Uint64 => :UInt64,
    :Void => Symbol(Nothing)
)

const _typedict = Dict{String,Type}()

function fixtypes(typ)
    whereall = []
    typ = fixtypes(typ, whereall)
    while !isempty(whereall)
        var = pop!(whereall)
        typ = Expr(:let, var, Expr(:call, Expr(:core, :UnionAll), var.args[1], typ))
    end
    return typ
end
fixtypes(typ, whereall) = typ
function fixtypes(typ::Expr, whereall::Vector{Any})
    if typ.head === :macrocall && typ.args[1] === _where_macrocall
        expand_where_macro(typ) # @where => TypeVar format forwards compatibility
    end
    if typ.head === :.
        if length(typ.args) == 2 && typ.args[1] === :Core
            arg = typ.args[2].value
            return Expr(:., :Core, QuoteNode(get(typemap_Core, arg, arg)))
        else
            return typ
        end
    elseif typ == :(Core.Type{TypeVar(:T,Union(Core.Any,Core.Undef))}) || typ == :(Core.Type{TypeVar(:T)})
        # Work around https://github.com/JuliaLang/julia/issues/8226 and the removal of Top
        return :(Core.Type)
    end

    for i = 1:length(typ.args)
        typ.args[i] = fixtypes(typ.args[i], whereall)
    end

    if (typ.head === :call && !isempty(typ.args) &&
        typ.args[1] === :TypeVar) # TypeVar => where format backwards compatibility
        tv = gensym()
        push!(whereall, Expr(:(=), tv, typ))
        return tv
    end

    if (typ.head === :call && !isempty(typ.args) &&
        typ.args[1] === :Union)
        typ = Expr(:curly, typ.args...)
    end

    if typ.head === :tuple
        if !any(x->isa(x,QuoteNode) || isbits(x), typ.args)
            # guess that we have a tuple type represented as a tuple
            typ = Expr(:curly, :Tuple, typ.args...)
        end
    end

    if typ.head === :curly
        # assume literal TypeVar should work like `T{<:S}`
        while !isempty(whereall)
            var = pop!(whereall)
            typ = Expr(:let, var, Expr(:call, Expr(:core, :UnionAll), var.args[1], typ))
        end
    end
    return typ
end

function _julia_type(s::AbstractString)
    typ = get(_typedict, s, UnconvertedType)
    if typ == UnconvertedType
        sp = Meta.parse(s, raise=false)
        if (isa(sp, Expr) && (sp.head == :error || sp.head == :continue || sp.head == :incomplete))
            println("error parsing type string ", s)
            eval(sp)
        end
        typ = julia_type(fixtypes(sp))
        if typ != UnsupportedType
            _typedict[s] = typ
        end
    end
    typ
end

function julia_type(e::Union{Symbol, Expr})
    if is_valid_type_ex(e)
        try # `try` needed to catch undefined symbols
            # `e` should be fully qualified, and thus reachable from Main
            typ = Core.eval(Main, e)
            typ == Type && return Type
            isa(typ, Type) && return typ
        catch
        end
    end
    return UnsupportedType
end

### Converting Julia types to fully qualified names
function full_typename(io::IO, file::JldFile, jltype::Union)
    print(io, "Union(")
    types = Base.uniontypes(jltype)
    if !isempty(types)
        full_typename(io, file, types[1])
        for i = 2:length(types)
            print(io, ',')
            full_typename(io, file, types[i])
        end
    end
    print(io, ')')
end
function full_typename(io::IO, file::JldFile, x::typeof(Union{}))
    print(io, "Union()")
end
function full_typename(io::IO, file::JldFile, x::UnionAll)
    x == Type && return print(io, "Type")
    print(io, "@where(")
    full_typename(io, file, x.body)
    print(io, ',')
    tv = x.var
    if tv.lb === Union{} && tv.ub === Any
        print(io, tv.name)
    elseif tv.lb === Union{}
        print(io, tv.name)
        print(io, "<:")
        full_typename(io, file, tv.ub)
    else
        full_typename(io, file, tv.lb)
        print(io, "<:")
        print(io, tv.name)
        print(io, "<:")
        full_typename(io, file, tv.ub)
    end
    print(io, ')')
end
function full_typename(io::IO, file::JldFile, tv::TypeVar)
    print(io, tv.name)
end
function full_typename(io::IO, ::JldFile, x)
    # Only allow bitstypes that show as AST literals and make sure that they
    # read back exactly as they are saved. Use show here (instead of print) to
    # preserve as many Julian type distinctions as we can e.g., 1 vs 0x01
    # A different implementation will be required to support custom immutables
    # or things as simple as Int16(1).
    s = sprint(show, x)
    if (isbits(x) || x isa Symbol || (x isa Tuple && all(y->(isbits(y) || y isa Symbol), x))) && eval(Meta.parse(s)) === x
        print(io, s)
    else
        error("type parameters with objects of type ", typeof(x), " are currently unsupported")
    end
end
function full_typename(io::IO, ::JldFile, x::Symbol)
    s = string(x)
    if occursin(" ", s)
        # escape spaces
        print_escaped(io, string("symbol(\"", string(x), "\")"), " ")
    else
        print(io, ":", x)
    end
end
function full_typename(io::IO, file::JldFile, jltype::DataType)
    mod = jltype.name.module
    if mod != Main
        mname = string(mod)
        for x in file.truncatemodules
            if startswith(mname, x)
                mname = length(x) == length(mname) ? "" : mname[sizeof(x)+1:end]
                if startswith(mname, ".")
                    mname = mname[2:end]
                end
                break
            end
        end

        if !isempty(mname)
            print(io, mname)
            print(io, '.')
        end
    end

    print(io, jltype.name.name)
    if !isempty(jltype.parameters)
        print(io, '{')
        full_typename(io, file, jltype.parameters[1])
        for i = 2:length(jltype.parameters)
            print(io, ',')
            full_typename(io, file, jltype.parameters[i])
        end
        print(io, '}')
    elseif jltype <: Tuple
        print(io, "{}")
    end
end
function full_typename(file::JldFile, x)
    io = IOBuffer(Vector{UInt8}(undef, 64), read=true, write=true)
    truncate(io, 0)
    full_typename(io, file, x)
    String(take!(io))
end

function truncate_module_path(file::JldFile, mod::Module)
    push!(file.truncatemodules, string(mod))
end

function Base.keys(parent::Union{JldFile, JldGroup})
    n = keys(parent.plain)
    keep = trues(length(n))
    reserved = [pathrefs[2:end], pathtypes[2:end], pathrequire[2:end], pathcreator[2:end]]
    for i = 1:length(n)
        if in(n[i], reserved)
            keep[i] = false
        end
    end
    n[keep]
end

function save_write(f, s, vname, wsession::JldWriteSession)
    if !isa(vname, Function)
        try
            write(f, s, vname)
        catch e
        end
    end
end

"""
```julia
@save "filename.jld" var1 [var2 var3 ...]
@save compress=true "filename.jld" var1 [var2 var3 ...]
```
Save the variables `var1` (and optionally `var2`, `var3`, etc.) to a JLD file "filename.jld".
Optionally you may specify `compress=true` or `compress=false` as the first argument,
to specify whether the resulting `.jld` file should be compressed (default=`false`).
"""
macro save(filename, vars...)
    if isempty(vars)
        # Save all variables in the current module
        writeexprs = Vector{Expr}(undef, 0)
        m = __module__
        for vname in names(m)
            s = string(vname)
            if !occursin(r"^_+[0-9]*$", s) && s != "ans" # skip IJulia history vars
                v = Core.eval(m, vname)
                if !isa(v, Module)
                    push!(writeexprs, :(write(f, $s, $(esc(vname)), wsession)))
                end
            end
        end
    else
        writeexprs = Vector{Expr}(undef, length(vars))
        for i = 1:length(vars)
            writeexprs[i] = :(write(f, $(string(vars[i])), $(esc(vars[i])), wsession))
        end
    end

    quote
        local f = jldopen($(esc(filename)), "w")
        wsession = JldWriteSession()
        try
            $(Expr(:block, writeexprs...))
        finally
            close(f)
        end
    end
end

macro save(opt::Expr, filename, vars...)
    compress = (opt.head === :(=) && first(opt.args) === :compress) ? last(opt.args) : false

    if isempty(vars)
        # Save all variables in the current module
        writeexprs = Vector{Expr}(undef, 0)
        m = __module__
        for vname in names(m)
            s = string(vname)
            if !occursin(r"^_+[0-9]*$", s) && s != "ans" # skip IJulia history vars
                v = Core.eval(m, vname)
                if !isa(v, Module)
                    push!(writeexprs, :(write(f, $s, $(esc(vname)), wsession)))
                end
            end
        end
    else
        writeexprs = Vector{Expr}(undef, length(vars))
        for i = 1:length(vars)
            writeexprs[i] = :(write(f, $(string(vars[i])), $(esc(vars[i])), wsession))
        end
    end

    quote
        local f = jldopen($(esc(filename)), "w", compress=$(esc(compress)))
        wsession = JldWriteSession()
        try
            $(Expr(:block, writeexprs...))
        finally
            close(f)
        end
    end
end

"""
```julia
@load "filename.jld"
@load "filename.jld" var1 [var2 var3 ...]
```
Load the variables `var1`, `var2`, et cetera, contained in the file
`filename.jld` into the current global scope.

If no variable names are specified, all variables from the file will be loaded.

Returns a `Vector` of `Symbol`s corresponding to the names of the loaded
variables.
"""
macro load(filename, vars...)
    if isempty(vars)
        if isa(filename, Expr)
            @warn("""@load-ing a file without specifying the variables to be loaded may produce
                     unexpected behavior unless the file is specified as a string literal. Future
                     versions of JLD will require that the file is specified as a string literal
                     in this case.""")
            filename = eval(__module__, filename)
        end
        # Load all variables in the top level of the file
        readexprs = Expr[]
        vars = Symbol[]
        f = jldopen(filename)
        try
            for n in keys(f)
                obj = f[n]
                try
                    if isa(obj, JldDataset)
                        push!(vars, Symbol(n))
                    end
                finally
                    close(obj)
                end
            end
        finally
            close(f)
        end
    end
    return quote
        f = jldopen($(esc(filename)))
        ($([esc(x) for x in vars]...),) = try
            ($([:(read(f, $(string(x)))) for x in vars]...),)
        finally
            close(f)
        end
        $(Symbol[v for v in vars]) # convert to Array
    end
end

# Save all the key-value pairs in the dict as top-level variables of the JLD
function fileio_save(f::File{format"JLD"}, dict::AbstractDict; compatible::Bool=false, compress::Bool=false)
    jldopen(FileIO.filename(f), "w"; compatible=compatible, compress=compress) do file
        wsession = JldWriteSession()
        for (k,v) in dict
            if !isa(k, AbstractString)
                error("Keys must be strings (the names of variables), got $k")
            end
            write(file, String(k), v, wsession)
        end
    end
end
# Or the names and values may be specified as alternating pairs
function fileio_save(f::File{format"JLD"}, name::AbstractString, value, pairs...; compatible::Bool=false, compress::Bool=false)
    if isodd(length(pairs)) || !isa(pairs[1:2:end], Tuple{Vararg{AbstractString}})
        throw(ArgumentError("arguments must be in name-value pairs"))
    end
    jldopen(FileIO.filename(f), "w"; compatible=compatible, compress=compress) do file
        wsession = JldWriteSession()
        write(file, String(name), value, wsession)
        for i=1:2:length(pairs)
            write(file, String(pairs[i]), pairs[i+1], wsession)
        end
    end
end
fileio_save(f::File{format"JLD"}, value...; kwargs...) = error("Must supply a name for each variable")

# load with just a filename returns a dictionary containing all the variables
function fileio_load(f::File{format"JLD"})
    jldopen(FileIO.filename(f), "r") do file
        Dict{String,Any}([(var, read(file, var)) for var in keys(file)])
    end
end
# When called with explicitly requested variable names, return each one
function fileio_load(f::File{format"JLD"}, varname::AbstractString)
    jldopen(FileIO.filename(f), "r") do file
        read(file, varname)
    end
end
fileio_load(f::File{format"JLD"}, varnames::AbstractString...) = load(f, varnames)
function fileio_load(f::File{format"JLD"}, varnames::Tuple{Vararg{AbstractString}})
    jldopen(FileIO.filename(f), "r") do file
        map((var)->read(file, var), varnames)
    end
end

# As of this version, packages aren't loaded into Main by default, so the root
# module check verifies that packages are still identified as being top level
# even if a binding to them is not present in Main.
_istoplevel(m::Module) = parentmodule(m) == Main || Base.is_root_module(m)

function addrequire(file::JldFile, mod::Module)
    _istoplevel(mod) || error("must be a toplevel module")
    addrequire(file, nameof(mod))
end

function addrequire(file::JldFile, modsym::Symbol)
    has_require = haskey(file, pathrequire)
    modules = has_require ? map(path2modsym, read(file, pathrequire)) : Symbol[]
    push!(modules, modsym)
    has_require && delete_object(file, pathrequire)
    write(file, pathrequire, modules)
end

function addrequire(file::JldFile, filename::AbstractString)
    @warn("\"addrequire(file, filename)\" is deprecated, please use \"addrequire(file, module)\"")
    addrequire(file, path2modsym(filename))
end

# Cope with JLD0.1 format
path2modsym(s::Symbol) = s
function path2modsym(filename::AbstractString)
    bname = basename(filename)
    if endswith(bname, ".jl")
        bname = bname[1:end-3]
    end
    Symbol(bname)
end

export
    addrequire,
    creator,
    ismmappable,
    jldopen,
    delete_object,
    readmmap,
    @load,
    @save,
    load,
    save,
    translate,
    truncate_module_path

const _runtime_properties = Ref{HDF5.Properties}()
compact_properties() = _runtime_properties[]

function __init__()
    prop = HDF5.DatasetCreateProperties(layout=:compact)

    global _runtime_properties
    _runtime_properties[] = prop

    Base.rehash!(_typedict, length(_typedict.keys))
    Base.rehash!(BUILTIN_TYPES.dict, length(BUILTIN_TYPES.dict.keys))

    nothing
end

include("JLD00.jl")

###
### v0.12.0 deprecations
###

export exists
@noinline function exists(p::Union{JldFile, JldGroup, JldDataset}, path::String)
    Base.depwarn("`exists(p, path)` is deprecated, use `haskey(p, path)` instead.", :exists)
    return haskey(p, path)
end

@noinline function Base.names(parent::Union{JldFile, JldGroup})
    Base.depwarn("`names(parent)` is deprecated, use `keys(parent)` instead.", :names)
    return keys(parent)
end

end
