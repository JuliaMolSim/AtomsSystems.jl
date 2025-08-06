

"""
    AtomicPropertySystem{D, LU, TB, TP, SM}

A system that contains atomic properties in addition to the basic attributes of atoms. 

This system is designed to work with `AbstractSimpleSystem` and allows for the storage of additional properties for each atom, such as custom mass, charge, or any other custom property.

It is recommeded to use `generic_system` to create instances of this system.

# Type Parameters
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit of length for positions.
- `TB`: Type of the base system, either `SimpleSystem` or `SimpleVelocitySystem`.
- `TP`: Type of the atomic properties (`NamedTuple`).
- `SM`: Boolean indicating if the system has a custom mass property.

# Fields
- `base_system`: The base system, which is either a `SimpleSystem` or a `SimpleVelocitySystem`.
- `atom_properties`: A vector of `TB` containing the properties of each atom.
"""
mutable struct AtomicPropertySystem{D, LU, TB, TP, SM} <: AbstractIsolatedSystem{D, LU}
    base_system::TB
    atom_properties::Vector{TP}
    function AtomicPropertySystem(
        sys::AbstractSimpleSystem{D, LU}, 
        properties::AbstractVector{<:NamedTuple}
    ) where {D, LU}
        @argcheck length(sys) == length(properties)
        new{D, LU, typeof(sys), eltype(properties), haskey(properties[1], :mass)}(sys, properties)
    end
end


AtomicPropertySystem(sys::AtomicPropertySystem) = deepcopy(sys)
function _AtomicPropertySystem(sys::AtomicPropertySystem, i)
    tmp = SimpleVelocitySystem(sys.base_system, i)
    return AtomicPropertySystem( tmp, sys.atom_properties[i] )
end
AtomicPropertySystem(sys::AtomicPropertySystem, i) = _AtomicPropertySystem( AtomicPropertySystem(sys), i )
AtomicPropertySystem(sys::AtomicPropertySystem, i::Colon) = _AtomicPropertySystem( AtomicPropertySystem(sys), i)
AtomicPropertySystem(sys::AtomicPropertySystem, i::BitVector) = _AtomicPropertySystem( AtomicPropertySystem(sys), i)

function AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector}, i)
    prop_names = [ k for k in atomkeys(sys) if ! in(k, (:species, :position, :velocity)) ]
    if length(prop_names) == 0
        return SimpleVelocitySystem(sys, i)
    end
    tmp = map( sys[i] ) do at
        NamedTuple( k=>at[k] for k in prop_names )
    end
    return AtomicPropertySystem( SimpleVelocitySystem(sys, i), tmp )
end

AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector}, ::Colon) = AtomicPropertySystem(sys, 1:length(sys))
AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector})          = AtomicPropertySystem(sys, :)
AtomicPropertySystem(sys::AbstractSimpleSystem, ::Colon)               = SimpleVelocitySystem(sys, :)
AtomicPropertySystem(sys::AbstractSimpleSystem, i)                     = SimpleVelocitySystem(sys, i)
AtomicPropertySystem(sys::AbstractSimpleSystem, i::BitVector)          = SimpleVelocitySystem(sys, i)
AtomicPropertySystem(sys::AbstractSimpleSystem)                        = SimpleVelocitySystem(sys)

AtomicPropertySystem(sys::AbstractSimpleSystem, spc::ChemicalSpecies...)  = SimpleVelocitySystem(sys, spc...)

function AtomicPropertySystem(ss::Union{AbstractSystem,AtomsVector}, spc::ChemicalSpecies...)
    i = map( x -> x in spc, species(ss, :) )
    return AtomicPropertySystem(ss, i)
end

function AtomicPropertySystem(sys::AtomicPropertySystem, spc::ChemicalSpecies...)
    i = map( x -> x in spc, species(sys, :) )
    return AtomicPropertySystem(sys, i)
end

function AtomicPropertySystem(sys::Union{AbstractSystem, AtomsVector}, i::BitVector)
    @argcheck length(sys) == length(i)
    j = (1:length(sys))[i]
    return AtomicPropertySystem(sys, j)
end


function Base.getindex(sys::AtomicPropertySystem, i::Int)
    tmp = sys.base_system[i]
    SimpleAtom( (; tmp.data..., sys.atom_properties[i]...) )
end

Base.getindex(sys::AtomicPropertySystem, x::Symbol) = sys.base_system[x]

function AtomsBase.atomkeys(sys::AbstractIsolatedSystem)
    base_keys = AtomsBase.atomkeys(sys.base_system)
    property_keys = _property_keys(sys)
    # there is a risk of key collisions, but this is not a problem in practice
    return (base_keys..., property_keys...)
end

_property_keys(sys::AbstractIsolatedSystem) = keys(sys.atom_properties[1]) 

function AtomsBase.mass(sys::AtomicPropertySystem{D, LU, TB, TP, true}, i::Int) where {D, LU, TB, TP}
    return getindex(sys.atom_properties[i], :mass)
end
function AtomsBase.mass(sys::AtomicPropertySystem{D, LU, TB, TP, true}, i) where {D, LU, TB, TP}
    return [ mass(sys, j)  for j in i ]
end
function AtomsBase.mass(sys::AtomicPropertySystem{D, LU, TB, TP, true}, ::Colon) where {D, LU, TB, TP}
    return mass(sys, 1:length(sys))
end

function AtomsBase.mass(sys::AtomicPropertySystem{D, LU, TB, TP, false}, i) where {D, LU, TB, TP}
    return mass(sys.base_system, i)
end

function Base.append!(sys1::T, sys2::T) where{T<:AtomicPropertySystem}
    Base.append!(sys1.base_system, sys2.base_system)
    Base.append!(sys1.atom_properties, sys2.atom_properties)
    return sys1
end

function Base.deleteat!(sys::AtomicPropertySystem, i)
    Base.deleteat!(sys.base_system, i)
    Base.deleteat!(sys.atom_properties, i)
    return sys
end