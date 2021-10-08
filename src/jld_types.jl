# Controls whether tuples and non-pointerfree immutables, which Julia
# stores as references, are stored inline in compound types when
# possible. Currently this is problematic because Julia fields of these
# types may be undefined.
const INLINE_TUPLE = false
const INLINE_POINTER_IMMUTABLE = false

const JLD_REF_TYPE = JldDatatype(HDF5.Datatype(HDF5.H5T_STD_REF_OBJ, false), 0)
const BUILTIN_TYPES = Set([Symbol, Type, BigFloat, BigInt])
const JL_TYPENAME_TRANSLATE = Dict{String, String}()
const JLCONVERT_INFO = Dict{Any, Any}()
const H5CONVERT_INFO = Dict{Any, Any}()

const EMPTY_TUPLE_TYPE = Tuple{}
const TypesType = Core.SimpleVector
TupleType{T<:Tuple} = Type{T}
tupletypes(T::TupleType) = T.parameters
typetuple(types) = Tuple{types...}

## Helper functions

translate(oldname::String, newname::String) = JL_TYPENAME_TRANSLATE[oldname] = newname

# Holds information about the mapping between a Julia and HDF5 type
struct JldTypeInfo
    dtypes::Vector{JldDatatype}
    offsets::Vector{Int}
    size::Int
end

# Get information about the HDF5 types corresponding to Julia types
function JldTypeInfo(parent::JldFile, types::TypesType, commit::Bool)
    dtypes = Vector{JldDatatype}(undef, length(types))
    offsets = Vector{Int}(undef, length(types))
    offset = 0
    for i = 1:length(types)
        dtype = dtypes[i] = h5fieldtype(parent, types[i], commit)
        offsets[i] = offset
        offset += HDF5.h5t_get_size(dtype)
    end
    JldTypeInfo(dtypes, offsets, offset)
end
JldTypeInfo(parent::JldFile, @nospecialize(T), commit::Bool) =
    JldTypeInfo(parent, Base.unwrap_unionall(T).types, commit)

# Write an HDF5 datatype to the file
function commit_datatype(parent::JldFile, dtype::HDF5.Datatype, @nospecialize(T))
    pparent = parent.plain
    if !haskey(pparent, pathtypes)
        gtypes = create_group(pparent, pathtypes)
    else
        gtypes = pparent[pathtypes]
    end

    id = length(gtypes)+1
    try
        HDF5.commit_datatype(gtypes, @sprintf("%08d", id), dtype)
    finally
        close(gtypes)
    end
    write_attribute(dtype, name_type_attr, full_typename(parent, T))

    # Store in map
    parent.jlh5type[T] = JldDatatype(dtype, id)
end

# If parent is nothing, we are creating the datatype in memory for
# validation, so don't commit it
commit_datatype(parent::Nothing, dtype::HDF5.Datatype, @nospecialize(T)) =
    JldDatatype(dtype, -1)

# The HDF5 library loses track of relationships among committed types
# after the file is saved. We mangle the names by appending a
# sequential identifier so that we can recover these relationships
# later.
mangle_name(jtype::JldDatatype, jlname) =
    jtype.index <= 0 ? string(jlname, "_") : string(jlname, "_", jtype.index)
Base.convert(::Type{HDF5.hid_t}, x::JldDatatype) = x.dtype.id

## Serialization of datatypes to JLD
##
## h5fieldtype - gets the JldDatatype corresponding to a given
## Julia type, when the Julia type is stored as an element of an HDF5
## compound type or array. This is the only function that can operate
## on non-concrete types.
##
## h5type - gets the JldDatatype corresponding to an object of the
## given Julia type. For pointerfree types, this is usually the same as
## the h5fieldtype.
##
## h5convert! - converts data from Julia to HDF5 in a buffer. Most
## methods are dynamically generated by gen_h5convert, but methods for
## special built-in types are predefined.
##
## jlconvert - converts data from HDF5 to a Julia object.
##
## jlconvert! - converts data from HDF5 to Julia in a buffer. This is
## only applicable in cases where fields of that type may not be stored
## as references (e.g., not plain types).

## Special types
##
## To create a special serialization of a datatype, one should:
##
## - Define a method of h5fieldtype that dispatches to h5type
## - Define a method of h5type that constructs the type
## - Define h5convert! and jlconvert
## - If the type is an immutable, define jlconvert!
## - Add the type to BUILTIN_TYPES

## HDF5 bits kinds

# This construction prevents these methods from getting called on type unions
const BitsKindTypes = Union{map(x->Type{x}, Base.uniontypes(HDF5.BitsType))...}

h5fieldtype(parent::JldFile, T::BitsKindTypes, ::Bool) =
    h5type(parent, T, false)

h5type(::JldFile, T::BitsKindTypes, ::Bool) =
    JldDatatype(HDF5.Datatype(HDF5.hdf5_type_id(T), false), 0)

h5convert!(out::Ptr, ::JldFile, x::T, ::JldWriteSession) where {T<:HDF5.BitsType} =
    unsafe_store!(convert(Ptr{T}, out), x)

_jlconvert_bits(::Type{T}, ptr::Ptr) where {T} = unsafe_load(convert(Ptr{T}, ptr))
_jlconvert_bits!(out::Ptr, ::Type{T}, ptr::Ptr) where {T} =
    (unsafe_store!(convert(Ptr{T}, out), unsafe_load(convert(Ptr{T}, ptr))); nothing)

jlconvert(T::BitsKindTypes, ::JldFile, ptr::Ptr) = _jlconvert_bits(T, ptr)
jlconvert!(out::Ptr, T::BitsKindTypes, ::JldFile, ptr::Ptr) = _jlconvert_bits!(out, T, ptr)

function jlconvert(T::Type{Char}, f::JldFile, ptr::Ptr)
    c = unsafe_load(convert(Ptr{Char}, ptr))
    if f.version < v"0.1.2"
        c = Char(reinterpret(UInt32, c))
    end
    return c
end
function jlconvert!(out::Ptr, T::Type{Char}, f::JldFile, ptr::Ptr)
    unsafe_store!(convert(Ptr{Char}, out), jlconvert(Char, f, ptr))
end

## Nothing

const NothingType = Type{Nothing}

jlconvert(T::NothingType, ::JldFile, ptr::Ptr) = nothing
jlconvert!(out::Ptr, T::NothingType, ::JldFile, ptr::Ptr) = (unsafe_store!(convert(Ptr{T}, out), nothing); nothing)
h5convert!(out::Ptr, ::JldFile, x::Nothing, ::JldWriteSession) = nothing

## Strings

h5fieldtype(parent::JldFile, ::Type{T}, ::Bool) where {T<:String} =
    h5type(parent, T, false)

# Stored as variable-length strings
function h5type(::JldFile, ::Type{T}, ::Bool) where T<:String
    type_id = HDF5.h5t_copy(HDF5.hdf5_type_id(T))
    HDF5.h5t_set_size(type_id, HDF5.H5T_VARIABLE)
    HDF5.h5t_set_cset(type_id, HDF5.cset(T))
    JldDatatype(HDF5.Datatype(type_id, false), 0)
end

# no-inline needed to ensure gc-root for x
@noinline h5convert!(out::Ptr, ::JldFile, x::String, ::JldWriteSession) =
    unsafe_store!(convert(Ptr{Ptr{UInt8}}, out), pointer(x))

function jlconvert(T::Union{Type{String}}, ::JldFile, ptr::Ptr)
    strptr = unsafe_load(convert(Ptr{Ptr{UInt8}}, ptr))
    str = unsafe_string(strptr)
    Libc.free(strptr)
    str
end

## Symbols

h5fieldtype(parent::JldFile, ::Type{Symbol}, commit::Bool) =
    h5type(parent, Symbol, commit)

# Stored as a compound type that contains a variable length string
function h5type(parent::JldFile, ::Type{Symbol}, commit::Bool)
    haskey(parent.jlh5type, Symbol) && return parent.jlh5type[Symbol]
    id = HDF5.h5t_create(HDF5.H5T_COMPOUND, 8)
    HDF5.h5t_insert(id, "symbol_", 0, h5fieldtype(parent, String, commit))
    dtype = HDF5.Datatype(id, parent.plain)
    commit ? commit_datatype(parent, dtype, Symbol) : JldDatatype(dtype, -1)
end

function h5convert!(out::Ptr, file::JldFile, x::Symbol, wsession::JldWriteSession)
    str = string(x)
    push!(wsession.persist, str)
    h5convert!(out, file, str, wsession)
end

jlconvert(::Type{Symbol}, file::JldFile, ptr::Ptr) = Symbol(jlconvert(String, file, ptr))


## BigInts and BigFloats

h5fieldtype(parent::JldFile, T::Union{Type{BigInt}, Type{BigFloat}}, commit::Bool) =
    h5type(parent, T, commit)

# Stored as a compound type that contains a variable length string
function h5type(parent::JldFile, T::Union{Type{BigInt}, Type{BigFloat}}, commit::Bool)
    haskey(parent.jlh5type, T) && return parent.jlh5type[T]
    id = HDF5.h5t_create(HDF5.H5T_COMPOUND, 8)
    HDF5.h5t_insert(id, "data_", 0, h5fieldtype(parent, String, commit))
    dtype = HDF5.Datatype(id, parent.plain)
    commit ? commit_datatype(parent, dtype, T) : JldDatatype(dtype, -1)
end

function h5convert!(out::Ptr, file::JldFile, x::BigInt, wsession::JldWriteSession)
    str = string(x, base=62)
    push!(wsession.persist, str)
    h5convert!(out, file, str, wsession)
end
function h5convert!(out::Ptr, file::JldFile, x::BigFloat, wsession::JldWriteSession)
    str = string(x)
    push!(wsession.persist, str)
    h5convert!(out, file, str, wsession)
end

jlconvert(::Type{BigInt}, file::JldFile, ptr::Ptr) =
    parse(BigInt, jlconvert(String, file, ptr), base = 62)
jlconvert(::Type{BigFloat}, file::JldFile, ptr::Ptr) =
    parse(BigFloat, jlconvert(String, file, ptr))

## Types

h5fieldtype(parent::JldFile, ::Type{T}, commit::Bool) where {T<:Type} =
    h5type(parent, Type, commit)

# Stored as a compound type that contains a variable length string
function h5type(parent::JldFile, ::Type{T}, commit::Bool) where T<:Type
    haskey(parent.jlh5type, Type) && return parent.jlh5type[Type]
    id = HDF5.h5t_create(HDF5.H5T_COMPOUND, 8)
    HDF5.h5t_insert(id, "typename_", 0, h5fieldtype(parent, String, commit))
    dtype = HDF5.Datatype(id, parent.plain)
    out = commit ? commit_datatype(parent, dtype, Type) : JldDatatype(dtype, -1)
end

function h5convert!(out::Ptr, file::JldFile, x::Type, wsession::JldWriteSession)
    str = full_typename(file, x)
    push!(wsession.persist, str)
    h5convert!(out, file, str, wsession)
end

jlconvert(::Type{T}, file::JldFile, ptr::Ptr) where {T<:Type} =
    julia_type(jlconvert(String, file, ptr))

## Union{}

h5fieldtype(parent::JldFile, ::Type{Union{}}, ::Bool) = JLD_REF_TYPE

## Arrays

# These show up as having T.size == 0, hence the need for
# specialization.
h5fieldtype(parent::JldFile, ::Type{Array{T,N}}, ::Bool) where {T,N} = JLD_REF_TYPE

## User-defined types
##
## Similar to special types, but h5convert!/jl_convert are dynamically
## generated.

## Tuples

if INLINE_TUPLE
    h5fieldtype(parent::JldFile, T::TupleType, commit::Bool) =
        isconcretetype(T) ? h5type(parent, T, commit) : JLD_REF_TYPE
else
    h5fieldtype(parent::JldFile, T::TupleType, ::Bool) = JLD_REF_TYPE
end

function h5type(parent::JldFile, T::TupleType, commit::Bool)
    haskey(parent.jlh5type, T) && return parent.jlh5type[T]
    # Tuples should always be concretely typed, unless we're
    # reconstructing a tuple, in which case commit will be false
    !commit || isconcretetype(T) || error("unexpected non-concrete type $T")

    typeinfo = JldTypeInfo(parent, T, commit)
    if isopaque(T)
        id = HDF5.h5t_create(HDF5.H5T_OPAQUE, opaquesize(T))
    else
        id = HDF5.h5t_create(HDF5.H5T_COMPOUND, typeinfo.size)
    end
    for i = 1:length(typeinfo.offsets)
        fielddtype = typeinfo.dtypes[i]
        HDF5.h5t_insert(id, mangle_name(fielddtype, i), typeinfo.offsets[i], fielddtype)
    end

    dtype = HDF5.Datatype(id, parent.plain)
    if commit
        jlddtype = commit_datatype(parent, dtype, T)
        if T == EMPTY_TUPLE_TYPE
            # to allow recovery of empty tuples, which HDF5 does not allow
            write_attribute(dtype, "empty", UInt8(1))
        end
        jlddtype
    else
        JldDatatype(dtype, -1)
    end
end

## All other objects

# For cases not defined above: If the type is mutable and non-empty,
# this is a reference. If the type is immutable, this is a type itself.
if INLINE_POINTER_IMMUTABLE
    h5fieldtype(parent::JldFile, @nospecialize(T), commit::Bool) =
        isconcretetype(T) && (!((VERSION >= v"1.7.0-DEV.1279") ? (Base.ismutabletype(T)) : (T.mutable)) || T.size == 0) ? h5type(parent, T, commit) : JLD_REF_TYPE
else
            h5fieldtype(parent::JldFile, @nospecialize(T), commit::Bool) =
        isconcretetype(T) && (!((VERSION >= v"1.7.0-DEV.1279") ? (Base.ismutabletype(T)) : (T.mutable)) || T.size == 0) && datatype_pointerfree(T) ? h5type(parent, T, commit) : JLD_REF_TYPE
end

function h5type(parent::JldFile, T::Type{Bool}, commit::Bool)
    if parent.version < v"0.1.3"
        h5type_default(parent, T, commit)
    else
        JldDatatype(HDF5.Datatype(HDF5.hdf5_type_id(T), false), 0)
    end
end

h5type(parent::JldFile, @nospecialize(T), commit::Bool) = h5type_default(parent, T, commit)

function h5type_default(parent::JldFile, @nospecialize(T), commit::Bool)
    !isa(T, DataType) && unknown_type_err(T)
    T = T::DataType

    haskey(parent.jlh5type, T) && return parent.jlh5type[T]
    isconcretetype(T) || error("unexpected non-concrete type ", T)

    if isopaque(T)
        # Empty type or non-basic bitstype
        id = HDF5.h5t_create(HDF5.H5T_OPAQUE, opaquesize(T))
    else
        # Compound type
        typeinfo = JldTypeInfo(parent, T.types, commit)
        id = HDF5.h5t_create(HDF5.H5T_COMPOUND, typeinfo.size)
        for i = 1:length(typeinfo.offsets)
            fielddtype = typeinfo.dtypes[i]
            HDF5.h5t_insert(id, mangle_name(fielddtype, fieldnames(T)[i]), typeinfo.offsets[i], fielddtype)
        end
    end

    dtype = HDF5.Datatype(id, parent.plain)
    if commit
        jlddtype = commit_datatype(parent, dtype, T)
        if T.size == 0
            # to allow recovery of empty types, which HDF5 does not allow
            write_attribute(dtype, "empty", UInt8(1))
        end
        jlddtype
    else
        JldDatatype(dtype, -1)
    end
end

# Normal objects
function _gen_jlconvert_type(typeinfo::JldTypeInfo, @nospecialize(T))
    ex = Expr(:block)
    args = ex.args
    push!(args, :(out = ccall(:jl_new_struct_uninit, Ref{T}, (Any,), T)))
    for i = 1:length(typeinfo.dtypes)
        h5offset = typeinfo.offsets[i]

        if HDF5.h5t_get_class(typeinfo.dtypes[i]) == HDF5.H5T_REFERENCE
            push!(args, quote
                ref = unsafe_load(convert(Ptr{HDF5.Reference}, ptr)+$h5offset)
                if ref != HDF5.Reference()
                    out.$(fieldnames(T)[i]) = convert($(T.types[i]), read_ref(file, ref))
                end
            end)
        else
            push!(args, :(out.$(fieldnames(T)[i]) = jlconvert($(T.types[i]), file, ptr+$h5offset)))
        end
    end
    push!(args, :(return out))
    return ex
end

function _gen_jlconvert_type!(typeinfo::JldTypeInfo, @nospecialize(T))
    error("unimplemented")
end

# Immutables
function _gen_jlconvert_immutable(typeinfo::JldTypeInfo, @nospecialize(T))
    ex = Expr(:block)
    args = ex.args
    jloffsets = map(idx->fieldoffset(T, idx), 1:fieldcount(T))
    if isbitstype(T)
        push!(args, :(out = Ref{T}()))
        push!(args, :(jlconvert!(unsafe_convert(Ptr{T}, out), T, file, ptr)))
        push!(args, :(return out[]))
    else
        nf = length(typeinfo.dtypes)
        push!(args, :(fieldvals = Vector{Any}(undef, $nf)))
        push!(args, :(ninit = 0))
        for i = 1:nf
            h5offset = typeinfo.offsets[i]
            jloffset = jloffsets[i]
            obj = gensym("obj")
            if isa(T.types[i], TupleType) && isbitstype(T.types[i])
                # We continue to store tuples as references for the sake of
                # backwards compatibility, but on 0.4 they are now stored
                # inline
                push!(args, quote
                    ref = unsafe_load(convert(Ptr{HDF5.Reference}, ptr)+$h5offset)
                    if ref == HDF5.Reference()
                        @warn("""A pointerfree tuple field was undefined.
                                 This is not supported in Julia 0.4 and the corresponding tuple will be uninitialized.""")
                    else
                        fieldvals[$i] = convert($(T.types[i]), read_ref(file, ref))
                        ninit += 1
                    end
                end)
            elseif HDF5.h5t_get_class(typeinfo.dtypes[i]) == HDF5.H5T_REFERENCE
                push!(args, quote
                    ref = unsafe_load(convert(Ptr{HDF5.Reference}, ptr)+$h5offset)
                    if ref != HDF5.Reference()
                        fieldvals[$i] = convert($(T.types[i]), read_ref(file, ref))
                        ninit += 1
                    end
                end)
            else
                push!(args, :(fieldvals[$i] = jlconvert($(T.types[i]), file, ptr+$h5offset)))
                push!(args, :(ninit += 1))
            end
        end
        push!(args, :(out = ccall(:jl_new_structv, Ref{T}, (Any,Ptr{Cvoid},UInt32), T, fieldvals, ninit)))
        push!(args, :(return out))
    end
    return ex
end

function _gen_jlconvert_immutable!(typeinfo::JldTypeInfo, @nospecialize(T))
    ex = Expr(:block)
    args = ex.args
    jloffsets = map(idx->fieldoffset(T, idx), 1:fieldcount(T))
    if isbitstype(T)
        for i = 1:length(typeinfo.dtypes)
            h5offset = typeinfo.offsets[i]
            jloffset = jloffsets[i]

            if isa(T.types[i], TupleType) && isbitstype(T.types[i])
                # We continue to store tuples as references for the sake of
                # backwards compatibility, but on 0.4 they are now stored
                # inline
                push!(args, quote
                    ref = unsafe_load(convert(Ptr{HDF5.Reference}, ptr)+$h5offset)
                    if ref == HDF5.Reference()
                        @warn("""A pointerfree tuple field was undefined.
                                 This is not supported in Julia 0.4 and the corresponding tuple will be uninitialized.""")
                    else
                        unsafe_store!(convert(Ptr{$(T.types[i])}, out)+$jloffset, read_ref(file, ref))
                    end
                end)
            elseif HDF5.h5t_get_class(typeinfo.dtypes[i]) == HDF5.H5T_REFERENCE
                error("reference encountered in pointerfree immutable; this is a bug")
            else
                push!(args, :(jlconvert!(out+$jloffset, $(T.types[i]), file, ptr+$h5offset)))
            end
        end
    else
        error("unimplemented")
    end
    push!(args, nothing)
    return ex
end

function _gen_jlconvert_tuple(typeinfo::JldTypeInfo, @nospecialize(T))
    # custom converter (rather than falling back to immutable datatype handling)
    # as we need to make sure to return a valid Tuple concrete type,
    # but some parameters of T may be Any
    ex = Expr(:block)
    args = ex.args
    tup = Expr(:tuple)
    tupargs = tup.args
    types = tupletypes(T)
    for i = 1:length(typeinfo.dtypes)
        h5offset = typeinfo.offsets[i]
        field = Symbol(string("field", i))

        if HDF5.h5t_get_class(typeinfo.dtypes[i]) == HDF5.H5T_REFERENCE
            push!(args, :($field = read_ref(file, unsafe_load(convert(Ptr{HDF5.Reference}, ptr)+$h5offset))))
        else
            push!(args, :($field = jlconvert($(types[i]), file, ptr+$h5offset)))
        end
        push!(tupargs, field)
    end
    push!(args, :(return $tup))
    return ex
end


function gen_jlconvert(typeinfo::JldTypeInfo, @nospecialize(T))
    T === Nothing && return
    # TODO: this is probably invalid, so try to do this differently
    JLCONVERT_INFO[T] = typeinfo
    nothing
end

function gen_jlconvert(@nospecialize(T))
    typeinfo = JLCONVERT_INFO[T]::JldTypeInfo
    if isa(T, TupleType)
        return _gen_jlconvert_tuple(typeinfo, T)
    elseif isempty(fieldnames(T))
            if T.size == 0 && !((VERSION >= v"1.7.0-DEV.1279") ? (Base.ismutabletype(T)) : (T.mutable))
                return T.instance
            else
               return :(_jlconvert_bits(T, ptr))
            end
        
    elseif T.size == 0
        return :(ccall(:jl_new_struct_uninit, Ref{T}, (Any,), T))
    elseif ((VERSION >= v"1.7.0-DEV.1279") ? (Base.ismutabletype(T)) : (T.mutable))
        return _gen_jlconvert_type(typeinfo, T)
    else
        return _gen_jlconvert_immutable(typeinfo, T)
    end
end

function gen_jlconvert!(@nospecialize(T))
    typeinfo = JLCONVERT_INFO[T]::JldTypeInfo
    if isa(T, TupleType)
        error("unimplemented")
    elseif isempty(fieldnames(T))
        if T.size == 0
            if !((VERSION >= v"1.7.0-DEV.1279") ? (Base.ismutabletype(T)) : (T.mutable))
                return nothing
            else
                return :(unsafe_store!(convert(Ptr{Any}, out), $(T.instance)); nothing)
            end
        else
           return :(_jlconvert_bits!(out, T, ptr); nothing)
        end
    elseif T.size == 0
        return nothing
    elseif ((VERSION >= v"1.7.0-DEV.1279") ? (Base.ismutabletype(T)) : (T.mutable))
        return _gen_jlconvert_type!(typeinfo, T)
    else
        return _gen_jlconvert_immutable!(typeinfo, T)
    end
end

@generated function jlconvert!(out::Ptr, ::Type{T}, file::JldFile, ptr::Ptr) where T
    return gen_jlconvert!(T)
end

@generated function jlconvert(::Type{T}, file::JldFile, ptr::Ptr) where T
    return gen_jlconvert(T)
end

## Common functions for all non-special types (including gen_h5convert)

# Whether this datatype should be stored as opaque
isopaque(t::TupleType) = t == EMPTY_TUPLE_TYPE
# isopaque(t::DataType) = isempty(fieldnames(t))
isopaque(t::DataType) = isa(t, TupleType) ? t == EMPTY_TUPLE_TYPE : isempty(fieldnames(t))

# The size of this datatype in the HDF5 file (if opaque)
opaquesize(t::TupleType) = 1
opaquesize(t::DataType) = max(1, t.size)

# Whether a type that is stored inline in HDF5 should be stored as a
# reference in Julia. This will only be called such that it returns
# true for some unions of special types defined above, unless either
# INLINE_TUPLE or INLINE_POINTER_IMMUTABLE is true.
uses_reference(T::DataType) = !datatype_pointerfree(T)
uses_reference(::TupleType) = true
uses_reference(::Union) = true

unknown_type_err(T) =
    error("""$T is not of a type supported by JLD
             Please report this error at https://github.com/JuliaIO/HDF5.jl""")

const BUILTIN_H5_types = Union{Nothing, Type, String, HDF5.BitsType, Symbol, BigInt, BigFloat}
function gen_h5convert(parent::JldFile, @nospecialize(T))
    T <: BUILTIN_H5_types && return
    # TODO: this is probably invalid, so try to do this differently
    haskey(H5CONVERT_INFO, T) && return
    H5CONVERT_INFO[T] = parent

    dtype = parent.jlh5type[T].dtype
    istuple = isa(T, TupleType)

    if isopaque(T)
        return
    end

    if istuple
        types = tupletypes(T)
    else
        types = (T::DataType).types
    end

    n = HDF5.h5t_get_nmembers(dtype.id)
    for i = 1:n
        if HDF5.h5t_get_member_class(dtype.id, i-1) != HDF5.H5T_REFERENCE
            gen_h5convert(parent, types[i])
        end
    end
    nothing
end

# There is no point in specializing this
function _gen_h5convert!(@nospecialize(T))
    parent = H5CONVERT_INFO[T]::JldFile
    dtype = parent.jlh5type[T].dtype
    istuple = isa(T, TupleType)

    if isopaque(T)
        if T.size == 0
            return nothing
        else
            return :(unsafe_store!(convert(Ptr{T}, out), x))
        end
    end

    if istuple
        types = tupletypes(T)
    else
        types = (T::DataType).types
    end

    getindex_fn = istuple ? (:getindex) : (:getfield)
    ex = Expr(:block)
    args = ex.args
    n = HDF5.h5t_get_nmembers(dtype.id)
    for i = 1:n
        offset = HDF5.h5t_get_member_offset(dtype.id, i-1)
        if HDF5.h5t_get_member_class(dtype.id, i-1) == HDF5.H5T_REFERENCE
            if istuple
                push!(args, :(unsafe_store!(convert(Ptr{HDF5.Reference}, out)+$offset,
                                            write_ref(file, $getindex_fn(x, $i), wsession))))
            else
                push!(args, quote
                    if isdefined(x, $i)
                        ref = write_ref(file, $getindex_fn(x, $i), wsession)
                    else
                        ref = HDF5.Reference()
                    end
                    unsafe_store!(convert(Ptr{HDF5.Reference}, out)+$offset, ref)
                end)
            end
        else
            push!(args, :(h5convert!(out+$offset, file, $getindex_fn(x, $i), wsession)))
        end
    end
    push!(args, nothing)
    return ex
end

@generated function h5convert!(out::Ptr, file::JldFile, x::T, wsession::JldWriteSession) where T
    return _gen_h5convert!(T)
end

## Find the corresponding Julia type for a given HDF5 type

# Type mapping function. Given an HDF5.Datatype, find (or construct) the
# corresponding Julia type.
function jldatatype(parent::JldFile, dtype::HDF5.Datatype)
    class_id = HDF5.h5t_get_class(dtype.id)
    if class_id == HDF5.H5T_STRING
        cset = HDF5.h5t_get_cset(dtype.id)
        if cset == HDF5.H5T_CSET_ASCII
            return String
        elseif cset == HDF5.H5T_CSET_UTF8
            return String
        else
            error("character set ", cset, " not recognized")
        end
    elseif class_id == HDF5.H5T_INTEGER || class_id == HDF5.H5T_FLOAT || class_id == HDF5.H5T_BITFIELD
        # This can be a performance hotspot
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_DOUBLE) > 0 && return Float64
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_INT64) > 0 && return Int64
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_FLOAT) > 0 && return Float32
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_INT32) > 0 && return Int32
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_UINT8) > 0 && return UInt8
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_UINT64) > 0 && return UInt64
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_UINT32) > 0 && return UInt32
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_INT8) > 0 && return Int8
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_INT16) > 0 && return Int16
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_UINT16) > 0 && return UInt16
        HDF5.h5t_equal(dtype.id, HDF5.H5T_NATIVE_B8) > 0 && return Bool
        error("unrecognized integer or float type")
    elseif class_id == HDF5.H5T_BITFIELD
        Bool
    elseif class_id == HDF5.H5T_COMPOUND || class_id == HDF5.H5T_OPAQUE
        addr = object_info(dtype).addr
        haskey(parent.h5jltype, addr) && return parent.h5jltype[addr]

        typename = read_attribute(dtype, name_type_attr)
        typename = get(JL_TYPENAME_TRANSLATE, typename, typename)
        T = julia_type(typename)
        if T == UnsupportedType
            @warn("type $typename not present in workspace; reconstructing")
            T = reconstruct_type(parent, dtype, typename)
        end

        T = readas(T)

        if !(T in BUILTIN_TYPES)
            # Call jldatatype on dependent types to validate them and
            # define jlconvert
            if class_id == HDF5.H5T_COMPOUND
                for i = 0:HDF5.h5t_get_nmembers(dtype.id)-1
                    member_name = HDF5.h5t_get_member_name(dtype.id, i)
                    idx = first(something(findlast("_", member_name), 0:-1))
                    if idx != sizeof(member_name)
                        member_dtype = open_datatype(parent.plain, string(pathtypes, '/', lpad(member_name[idx+1:end], 8, '0')))
                        jldatatype(parent, member_dtype)
                    end
                end
            end

            gen_jlconvert(JldTypeInfo(parent, T, false), T)
        end

        # Verify that types match
        newtype = h5type(parent, T, false).dtype
        if T !== Expr  # in julia >= 0.7 Expr has 1 fewer fields; the trailing field can be ignored.
            if dtype == newtype
                # If the types compare equal, continue
            else
                # Otherwise the types may be compatible yet still naively compare as
                # unequal. This seems to be related to H5T_REFERENCE members in the compound
                # datatype (starting with libhdf5 v1.12) comparing unequal if one half has
                # been committed while the other is still transient.
                #
                # Both dtype and newtype may be committed (former) or contain committed member
                # types loaded from the cache (latter), so copy both for comparison.
                dtype′  = HDF5.Datatype(HDF5.h5t_copy(dtype.id))
                newtype = HDF5.Datatype(HDF5.h5t_copy(newtype.id))
                dtype′ == newtype || throw(TypeMismatchException(typename))
            end
        end

        # Store type in type index
        index = typeindex(parent, addr)
        parent.jlh5type[T] = JldDatatype(dtype, index)
        parent.h5jltype[addr] = T
        T
    else
        error("unrecognized HDF5 datatype class ", class_id)
    end
end

# Create a Julia type based on the HDF5.Datatype from the file. Used
# when the type is no longer available.
function reconstruct_type(parent::JldFile, dtype::HDF5.Datatype, savedname::AbstractString)
    name = gensym(savedname)
    class_id = HDF5.h5t_get_class(dtype.id)
    if class_id == HDF5.H5T_OPAQUE
        if haskey(dtype, "empty")
            @eval (struct $name; end; $name)
        else
            sz = Int(HDF5.h5t_get_size(dtype.id))*8
            @eval (primitive type $name $sz end; $name)
        end
    else
        # Figure out field names and types
        nfields = HDF5.h5t_get_nmembers(dtype.id)
        fieldnames = Vector{Symbol}(undef, nfields)
        fieldtypes = Vector{Type}(undef, nfields)
        for i = 1:nfields
            membername = HDF5.h5t_get_member_name(dtype.id, i-1)
            idx = first(something(findlast("_", membername), 0:-1))
            fieldname = fieldnames[i] = Symbol(membername[1:prevind(membername,idx)])

            if idx != sizeof(membername)
                # There is something past the underscore in the HDF5 field
                # name, so the type is stored in file
                memberdtype = open_datatype(parent.plain, string(pathtypes, '/', lpad(membername[idx+1:end], 8, '0')))
                fieldtypes[i] = jldatatype(parent, memberdtype)
            else
                memberclass = HDF5.h5t_get_member_class(dtype.id, i-1)
                if memberclass == HDF5.H5T_REFERENCE
                    # Field is a reference, so use Any
                    fieldtypes[i] = Any
                else
                    # Type is built-in
                    memberdtype = HDF5.Datatype(HDF5.h5t_get_member_type(dtype.id, i-1), parent.plain)
                    fieldtypes[i] = jldatatype(parent, memberdtype)
                end
            end
        end

        if startswith(savedname, "(") || startswith(savedname, "Core.Tuple{")
            # We're reconstructing a tuple
            typetuple(fieldtypes)
        else
            # We're reconstructing some other type
            @eval begin
                struct $name
                    $([:($(fieldnames[i])::$(fieldtypes[i])) for i = 1:nfields]...)
                end
                $name
            end
        end
    end
end

# Get the index of a type in the types group. This could be cached, but
# it's already many times faster than calling H5Iget_name with a lot of
# data in the file, and it only needs to be called once per type.
# Revisit if this ever turns out to be a bottleneck.
function typeindex(parent::JldFile, addr::HDF5.haddr_t)
    gtypes = parent.plain[pathtypes]
    i = 1
    for x in gtypes
        if object_info(x).addr == addr
            return i
        end
        i += 1
    end
end
