
"""
    SimpleSystem{D, LU, TP} <: AbstractSimpleSystem{D, LU}

A simple system that holds only species and positions of atoms, without any additional properties.

SimpleSystem has `IsolatedCell` as its cell.

No intended to be called directly. Instead use `generic_system` to build systems.

# Type Parameters
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit of length for positions (LU=unit(TP)).
- `TP`: Type of position vector, in form of `SVector{D, TP}` where `TP` is a unitful length type.

# Fields
- `species`: A vector of `ChemicalSpecies` representing the species of each atom.
- `position`: A vector of position vectors, each of type `SVector{D, TP}`.
"""
mutable struct SimpleSystem{D, LU, TP} <: AbstractSimpleSystem{D, LU}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    function SimpleSystem(species::AbstractVector{ChemicalSpecies}, r::AbstractVector{SVector{D, TP}}) where{D,TP<:Unitful.Length}
        @argcheck length(species) == length(r)
        LU = unit(TP)
        new{D, LU, TP}(species, r)
    end
end

function SimpleSystem(species::AbstractVector{<:AtomsBase.AtomId}, r::AbstractVector{<:AbstractVector})
    tmp = map( x -> SVector(x...), r )
    spc = map( x -> ChemicalSpecies(x), species )
    return SimpleSystem(spc, tmp)
end

function SimpleSystem(species::AtomsBase.AtomId, pos::AbstractVector{<:Unitful.Length})
    return SimpleSystem([species], [pos])
end

SimpleSystem(sys::SimpleSystem) = deepcopy(sys)

function SimpleSystem(ss::Union{AbstractSystem,AtomsVector}, i)
    return SimpleSystem(species(ss, i), position(ss, i))
end

SimpleSystem(sys::Union{AbstractSystem, AtomsVector}) = SimpleSystem(sys, :)

function SimpleSystem(ss::Union{AbstractSystem,AtomsVector}, spc::ChemicalSpecies...)
    i = map( x -> x in spc, species(ss, :) )
    return SimpleSystem(ss, i)
end


Base.getindex(ss::SimpleSystem, i::Int) = SimpleAtom(ss.species[i], ss.position[i])
function Base.getindex(ss::AbstractSimpleSystem, x::Symbol)
    if x === :cell_vectors
        return cell_vectors(ss)
    elseif x === :periodicity
        return periodicity(ss)
    else
        throw( KeyError(x) )
    end
end


Base.length(ss::AbstractSimpleSystem) = length(ss.species)

AtomsBase.atomkeys(::SimpleSystem) = (:position, :species)
AtomsBase.cell(::SimpleSystem{D, LU, TP}) where{D, LU, TP} = IsolatedCell(D, TP)


function Base.append!(sys1::SimpleSystem{D, LU, TP}, sys2::SimpleSystem{D, LU, TP}) where{D, LU, TP}
    Base.append!(sys1.position, sys2.position)
    Base.append!(sys1.species, sys2.species)
    return sys1
end

function Base.deleteat!(sys::SimpleSystem, i)
    Base.deleteat!(sys.species, i)
    Base.deleteat!(sys.position, i)
    return sys
end


##
"""
    SimpleVelocitySystem{D, LU, UV, TP, TV} <: AbstractSimpleSystem{D, LU}

System that holds species, positions, and velocities of atoms.

# Type Parameters
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit of length for positions (LU=unit(TP)).
- `UV`: Unit of velocity for velocities (UV=unit(TV)).
- `TP`: Type of position vector, in form of `SVector{D, TP}` where `TP` is a unitful length type.
- `TV`: Type of velocity vector, in form of `SVector{D, TV}` where `TV` is a unitful velocity type.

No intended to be called directly. Instead use `generic_system` to build systems.

# Fields
- `species`: A vector of `ChemicalSpecies` representing the species of each atom.
- `position`: A vector of position vectors, each of type `SVector{D, TP}`.
- `velocity`: A vector of velocity vectors, each of type `SVector{D, TV}`.
"""
mutable struct SimpleVelocitySystem{D, LU, UV, TP, TV} <: AbstractSimpleSystem{D, LU}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    velocity::Vector{SVector{D, TV}}
    function SimpleVelocitySystem(
        species::AbstractVector{ChemicalSpecies}, 
        r::AbstractVector{SVector{D, TP}}, 
        v::AbstractVector{SVector{D, TV}}
    )   where {D, TP<:Unitful.Length, TV<:Unitful.Velocity}
        @argcheck length(species) == length(r) == length(v)
        LU = unit(TP)
        UV = unit(TV)
        new{D, LU, UV, TP, TV}(species, r, v)
    end
end

function SimpleVelocitySystem(
        species::AbstractVector{<:AtomsBase.AtomId}, 
        r::AbstractVector{<:AbstractVector}, 
        v::AbstractVector{<:AbstractVector}
    )
    tmp_r = map( x -> SVector(x...), r )
    tmp_v = map( x -> SVector(x...), v )
    spc = map( x -> ChemicalSpecies(x), species )
    return SimpleVelocitySystem(spc, tmp_r, tmp_v)
end


function SimpleVelocitySystem(sys::Union{AbstractSystem, AtomsVector}, i)
    if hasatomkey(sys, :velocity)
        scp = species(sys, i)
        pos = position(sys, i)
        vel = velocity(sys, i)
        return SimpleVelocitySystem(scp, pos, vel)
    else
        return SimpleSystem(sys, i)
    end
end

SimpleVelocitySystem(sys::Union{AbstractSystem, AtomsVector}) = SimpleVelocitySystem(sys, :)
SimpleVelocitySystem(sys::SimpleVelocitySystem) = deepcopy(sys)

function SimpleVelocitySystem(
    species::AtomsBase.AtomId, 
    pos::AbstractVector{<:Unitful.Length}, 
    vel::AbstractVector{<:Unitful.Velocity}
)
    return SimpleVelocitySystem([species], [pos], [vel])
end

function SimpleVelocitySystem(ss::Union{AbstractSystem,AtomsVector}, spc::ChemicalSpecies...)
    i = map( x -> x in spc, species(ss, :) )
    return SimpleVelocitySystem(ss, i)
end

Base.getindex(ss::SimpleVelocitySystem, i::Int) = SimpleAtom(ss.species[i], ss.position[i], ss.velocity[i])

AtomsBase.atomkeys(::SimpleVelocitySystem) = (:position, :velocity, :species)
AtomsBase.cell(::SimpleVelocitySystem{D, LU, UV, TP, TV}) where{D, LU,UV, TP, TV} = IsolatedCell(D, TP)
AtomsBase.velocity(sys::SimpleVelocitySystem, i::Int) = sys.velocity[i]
AtomsBase.velocity(sys::SimpleVelocitySystem, i) = view(sys.velocity, i)
AtomsBase.velocity(sys::SimpleVelocitySystem, ::Colon) = sys.velocity

AtomsBase.set_velocity!(sys::SimpleVelocitySystem, i, v) = setindex!(sys.velocity, v, i)

function Base.append!(
    sys1::SimpleVelocitySystem{D, LU, UV, TP, TV}, 
    sys2::SimpleVelocitySystem{D, LU, UV, TP, TV}
) where{D, LU, UV, TP, TV}
    Base.append!(sys1.position, sys2.position)
    Base.append!(sys1.species, sys2.species)
    Base.append!(sys1.velocity, sys2.velocity)
    return sys1
end


function Base.deleteat!(sys::SimpleVelocitySystem, i)
    Base.deleteat!(sys.species, i)
    Base.deleteat!(sys.position, i)
    Base.deleteat!(sys.velocity, i)
    return sys
end