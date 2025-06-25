# NOTE SimpleAtom could be improved by designing better constructors.
# Current ones have some tests for consitency or avoid them completely.

"""
    SimpleAtom(id::AtomsBase.AtomId, r::SVector; kwargs...)
    SimpleAtom(id::AtomsBase.AtomId, r::SVector, v::SVector; kwargs...)
    SimpleAtom(sa::SimpleAtom; kwargs...

Stucture to represent a single atom with species, position, and optional velocity and other properties.

The `SimpleAtom` is bitstype when the data is bitstypes.

Mass is set to the mass of the species by default, but can be overridden with a keyword argument.

# Examples
```julia
SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å")
SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å", [0.1, 0.0, 0.0]u"Å/s"; mass = 16.0u"u", charge = -1.0u"q")
SimpleAtom(ChemicalSpecies(:H), [0.0, 0.0, 0.0]u"Å")
SimpleAtom( :O => [1.0, 0.0, 0.0]u"Å" )
```
"""
struct SimpleAtom{D, TD, TP}
    data::TD
    function SimpleAtom(data::NamedTuple)
        @argcheck haskey(data, :species) && isa(data.species, ChemicalSpecies)
        @argcheck haskey(data, :position) && eltype(data.position) <: Unitful.Length
        @argcheck ! ( haskey(data, :velocity) && length(data.velocity) != length(data.position) )
        new{length(data.position), typeof(data), eltype(data.position)}(data)
    end
    function SimpleAtom(spc::ChemicalSpecies, r::SVector{D, TP}) where{D, TP<:Unitful.Length}
        tmp = ( species=spc, position=r)
        new{D, typeof(tmp), TP}(tmp)
    end
    function SimpleAtom(
        spc::ChemicalSpecies,
        r::SVector{D, TP},
        v::SVector{D, <:Unitful.Velocity};
        kwargs...
    ) where{D, TP<:Unitful.Length}
        tmp = ( species=spc, position=r, velocity=v, kwargs... )
        new{D, typeof(tmp), TP}(tmp)
    end
end

function SimpleAtom(
    spc::ChemicalSpecies, 
    r::AbstractVector{<:Unitful.Length}; 
    kwargs...
)   
    if haskey(kwargs, :velocity)
        if length(kwargs[:velocity]) != length(r) && !isa(v, SVector)
            throw( ArgumentError("Position and velocity have same length") )
        end
        v = SVector( kwargs[:velocity]... )
        filter!(x->x!=:velocity, keys(kwargs))
        tmp = ( species=spc, position=SVector(r...), velocity=v, kwargs... )
    end
    tmp = ( species=spc, position=SVector(r...), kwargs... )
    return SimpleAtom(tmp)
end

function SimpleAtom(at::Union{AtomsBase.Atom, AtomsBase.AtomView})
    spc = species(at)
    r = position(at)
    if haskey(at, :velocity)
        v = velocity(at)
        properties = NamedTuple( k=>at[k] for k in keys(at) if ! in(k, (:species, :position, :velocity)))
        return SimpleAtom(spc, SVector(r...), SVector(v...); properties...)
    end
    properties = NamedTuple( k=>at[k] for k in keys(at) if ! in(k, (:species, :position)))
    return SimpleAtom(spc, r; properties...)
end

function SimpleAtom(sa::SimpleAtom; kwargs...)
    length(kwargs) == 0 && return sa
    tmp = Dict{Symbol, Any}( pairs(sa) )
    foreach( pairs(kwargs) ) do (k,v)
        tmp[k] = v
    end
    return SimpleAtom( NamedTuple( tmp ) )
end

function SimpleAtom(id::AtomsBase.AtomId, r::AbstractVector{<:Unitful.Length}; kwargs...)
    return SimpleAtom(ChemicalSpecies(id), r; kwargs...)
end


function SimpleAtom(
        id::AtomsBase.AtomId,
        r::AbstractVector{<:Unitful.Length},
        v::AbstractVector{<:Unitful.Velocity};
        kwargs...
    )
    return SimpleAtom(ChemicalSpecies(id), SVector(r...), SVector(v...); kwargs...)
end

function SimpleAtom(pr::Pair; kwargs...)
    return SimpleAtom(pr[1], pr[2]; kwargs...)
end

AtomsBase.n_dimensions(::SimpleAtom{D, TD}) where {D, TD} = D

AtomsBase.velocity(sa::SimpleAtom) = haskey(sa, :velocity) ? sa.data.velocity : missing
AtomsBase.position(sa::SimpleAtom) = sa.data.position
AtomsBase.mass(sa::SimpleAtom)     = haskey(sa, :mass) ? sa.data.mass : AtomsBase.mass(AtomsBase.species(sa))
AtomsBase.species(sa::SimpleAtom)  = sa.data.species

AtomsBase.atom_name(sa::SimpleAtom)     = atom_name(species(sa))
AtomsBase.atomic_symbol(sa::SimpleAtom) = atomic_symbol(species(sa))
AtomsBase.atomic_number(sa::SimpleAtom) = atomic_number(species(sa))
AtomsBase.element(sa::SimpleAtom)       = element(species(sa))

Base.getindex(sa::SimpleAtom, x::Symbol) = x == :mass ? mass(sa) : sa.data[x]
Base.haskey(sa::SimpleAtom, x::Symbol)   = Base.haskey(sa.data, x)
Base.keys(sa::SimpleAtom)                = Base.keys(sa.data)
Base.pairs(sa::SimpleAtom)               = pairs(sa.data)

Base.show(io::IO, sa::SimpleAtom) = AtomsBase.show_atom(io, sa)
Base.show(io::IO, mime::MIME"text/plain", sa::SimpleAtom) = AtomsBase.show_atom(io, mime, sa)


## Vector properties

const AtomsVector{D,TD, TP} = AbstractVector{SimpleAtom{D, TD, TP}}

AtomsBase.atomkeys(sys::AtomsVector) = keys(sys[1])
AtomsBase.hasatomkey(sys::AtomsVector, x::Symbol) = haskey(sys[1], x)

AtomsBase.cell(::AtomsVector{D, TD, TP}) where{D, TD, TP} = IsolatedCell(D, TP)
AtomsBase.cell_vectors(sys::AtomsVector) = cell_vectors(cell(sys))
AtomsBase.periodicity(sys::AtomsVector) = periodicity(cell(sys))

AtomsBase.mass(sys::AtomsVector, i) = map(j->sys[j][:mass], i)
AtomsBase.mass(sys::AtomsVector, ::Colon) = map(j->sys[j][:mass], eachindex(sys))

AtomsBase.position(sys::AtomsVector, i) = map(j->sys[j][:position], i)
AtomsBase.position(sys::AtomsVector, ::Colon) = map(j->sys[j][:position], eachindex(sys))

AtomsBase.species(sys::AtomsVector, i) = map(j->sys[j][:species], i)
AtomsBase.species(sys::AtomsVector, ::Colon) = map(j->sys[j][:species], eachindex(sys))

AtomsBase.velocity(sys::AtomsVector, i) = map(j->sys[j][:velocity], i)
AtomsBase.velocity(sys::AtomsVector, ::Colon) = map(j->sys[j][:velocity], eachindex(sys))
