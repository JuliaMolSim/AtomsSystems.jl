
mutable struct SimpleTrajectory{D, LU, TP} <: AbstractSimpleTrajectory{D, LU, TP}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    n_atoms::Int
    function SimpleTrajectory(
        species::AbstractVector{ChemicalSpecies},
        positions::AbstractVector{SVector{D, TP}}
    ) where {D, TP<:Unitful.Length}
        @argcheck length(species) > 0 "Species vector must not be empty"
        n_atoms = length(species)
        @argcheck length(positions) == n_atoms "Number of positions must match number of species"
        LU = unit(TP)
        new{D, LU, TP}(deepcopy(species), deepcopy(positions), n_atoms)
    end
end


function SimpleTrajectory(sys)
    return SimpleTrajectory(
        species(sys, :),
        position(sys, :)
    )
end

function SimpleTrajectory(traj::AbstractVector{<:AbstractSystem})
    first_frame = traj[begin]
    tmp = SimpleTrajectory(
        species(first_frame, :),
        position(first_frame, :)
    )
    for frame in Iterators.drop(traj, 1)
        append!(tmp, frame)
    end
    return tmp
end


function Base.append!(traj::SimpleTrajectory{D}, pos::AbstractVector{SVector{D, TP}}) where {D, TP<:Unitful.Length}
    @argcheck length(pos) % n_atoms(traj) == 0 "Position vector must have a length that is a multiple of the number of atoms in the trajectory"
    append!(traj.position, pos)
    return traj
end

function Base.append!(traj::SimpleTrajectory{D}, sys::AbstractSystem{D}) where {D}
    @argcheck length(sys) == n_atoms(traj) "System must have the same number of atoms as the trajectory"
    @argcheck all(species(sys, :) .=== traj.species) "Species must match the trajectory species"
    return Base.append!(traj, position(sys, :))    
end


function Base.eltype(::SimpleTrajectory{D, LU, TP}) where {D, LU, TP}
    return AtomsSystems.SimpleSystemView{D, LU, TP, Tuple{UnitRange{Int}}, true}
end


function Base.getindex(traj::SimpleTrajectory, i::Int)
    @argcheck i >= 1 && i <= length(traj) BoundsError(traj, i)
    k = traj.n_atoms * (i-1) +1  : traj.n_atoms * i
    return AtomsSystems.SimpleSystemView(view(traj.species, 1:traj.n_atoms), view(traj.position, k))
end

Base.size(traj::AbstractSimpleTrajectory) = (Int(length(traj.position)//traj.n_atoms), )


Base.show(io::IO, trj::SimpleTrajectory) =
    print(io, "SimpleTrajectory with ", length(trj), " frames of ", n_atoms(trj), " atoms")

##

mutable struct SimpleVelocityTrajectory{D, LU, TP, TV} <: AbstractSimpleTrajectory{D, LU, TP}
    species::Vector{ChemicalSpecies}
    position::Vector{SVector{D, TP}}
    velocity::Vector{SVector{D, TV}}
    n_atoms::Int
    function SimpleVelocityTrajectory(
        species::AbstractVector{ChemicalSpecies}, 
        positions::AbstractVector{SVector{D, TP}}, 
        velocities::AbstractVector{SVector{D, TV}}
        ) where {D, TP<:Unitful.Length, TV<:Unitful.Velocity}
        @argcheck length(species) > 0 "Species vector must not be empty"
        n_atoms = length(species)
        @argcheck length(positions) == n_atoms "Number of positions must match number of species"
        @argcheck length(velocities) == n_atoms "Number of velocities must match number of species"
        LU = unit(TP)
        new{D, LU, TP, TV}(deepcopy(species), deepcopy(positions), deepcopy(velocities), n_atoms)
    end
end


function SimpleVelocityTrajectory(sys::AbstractSystem)
    if hasatomkey(sys, :velocity)
        return SimpleVelocityTrajectory(
            species(sys, :),
            position(sys, :),
            velocity(sys, :)
        )
    end
    return SimpleTrajectory(sys)
end

function SimpleVelocityTrajectory(traj::AbstractVector{<:AbstractSystem})
    first_frame = traj[begin]
    if hasatomkey(first_frame, :velocity)
        tmp = SimpleVelocityTrajectory(
            species(first_frame, :),
            position(first_frame, :),
            velocity(first_frame, :)
        )
        for frame in Iterators.drop(traj, 1)
            append!(tmp, frame)
        end
        return tmp
    end
    return SimpleTrajectory(traj)
end

function Base.append!(
    traj::SimpleVelocityTrajectory{D},
    pos::AbstractVector{SVector{D, TP}},
    vel::AbstractVector{SVector{D, TV}}
) where {D, TP<:Unitful.Length, TV<:Unitful.Velocity}
    @argcheck length(pos) % n_atoms(traj) == 0 "Position vector must have a length that is a multiple of the number of atoms in the trajectory"
    @argcheck length(vel) == length(pos) "Velocity vector must have the same length as position vector"
    append!(traj.position, pos)
    append!(traj.velocity, vel)
    return traj
end

function Base.append!(traj::SimpleVelocityTrajectory{D}, sys::AbstractSystem{D}) where {D}
    @argcheck length(sys) == n_atoms(traj) "System must have the same number of atoms as the trajectory"
    @argcheck all(species(sys, :) .=== traj.species) "Species must match the trajectory species"
    return Base.append!(traj, position(sys, :), velocity(sys, :))    
end

function Base.eltype(::SimpleVelocityTrajectory{D, LU, TP, TV}) where {D, LU, TP, TV}
    return AtomsSystems.SimpleVelocitySystemView{D, LU, TP, TV, Tuple{UnitRange{Int}}, true}
end

function Base.getindex(traj::SimpleVelocityTrajectory, i::Int)
    @argcheck i >= 1 && i <= length(traj) BoundsError(traj, i)
    k = traj.n_atoms * (i-1) +1  : traj.n_atoms * i
    return AtomsSystems.SimpleVelocitySystemView(
        view(traj.species, 1:traj.n_atoms), 
        view(traj.position, k), 
        view(traj.velocity, k)
    )
end
    
        
Base.show(io::IO, trj::SimpleVelocityTrajectory) =
    print(io, "SimpleVelocityTrajectory with ", length(trj), " frames of ", n_atoms(trj), " atoms")