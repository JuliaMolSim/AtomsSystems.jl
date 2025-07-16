
"""
    VariableVolumeTrajectory{D, LU, TP, TB} <: AbstractCellTrajectory{D, LU, TP}

A trajectory that allows volume changes.
This is done by storing cell for each frame separately.

Cells are `PeriodicCell`. For isolates cell trajectories the return type is either
`SimpleTrajectory` or `SimpleVelocityTrajectory`.

Chemfiles extension provides a way to read trajectories from files.

# Type Parameters
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit of length for positions.
- `TP`: Type of positions `Unitful.Length`.
- `TB`: Type of the base trajectory, which is a `SimpleVelocityTrajectory` or `SimpleTrajectory`.

# Fields
- `base_trajectory`: The base trajectory, which is either a `SimpleVelocityTrajectory` or a `SimpleTrajectory`.
- `cell`: A vector of periodic cells, one for each frame

# Example
```julia
traj = VariableVolumeTrajectory(first_frame)
push!(traj, second_frame)

# Chemfiles extension
using Chemfiles
traj = ConstantVolumeTrajectory("path/to/trajectory")
```
"""
mutable struct VariableVolumeTrajectory{D, LU, TP, TB} <: AbstractCellTrajectory{D, LU, TP}
    base_trajectory::TB
    cell::Vector{PeriodicCell{D, TP}}
    function VariableVolumeTrajectory(
        traj::AbstractSimpleTrajectory{D, LU, TP},
        cell::PeriodicCell{D, TP}
    ) where {D, LU, TP}
        new{D, LU, TP, typeof(traj)}(
            traj,
            [cell]
        )
    end
end

function VariableVolumeTrajectory(sys::AbstractSystem)
    base = SimpleVelocityTrajectory(sys)
    if ! (cell(sys) isa PeriodicCell)
        return base
    end
    return VariableVolumeTrajectory(base, cell(sys))
end

function VariableVolumeTrajectory(systems::AbstractVector{<:AbstractSystem})
    tmp = VariableVolumeTrajectory(systems[begin])
    for frame in Iterators.drop(systems, 1)
        push!(tmp, frame)
    end
    return tmp
    
end


Base.size(sys::AbstractCellTrajectory) = size(sys.base_trajectory)

function Base.getindex(sys::VariableVolumeTrajectory, i::Int)
    @argcheck i >= 1 && i <= length(sys) BoundsError(sys, i)
    return AtomsSystems.CellSystemView(
        sys.base_trajectory[i],
        sys.cell[i]
    )
end

function Base.push!(traj::VariableVolumeTrajectory{D, LU, TP}, sys::AbstractSystem{D}) where {D, LU, TP}
    @argcheck length(sys) == n_atoms(traj) "System must have the same number of atoms as the trajectory"
    if cell(sys) isa PeriodicCell{D, TP}
        push!(traj.base_trajectory, sys)
        push!(traj.cell, cell(sys))
    else
        error("System must have a periodic cell to push to a VariableVolumeTrajectory") 
    end
    return traj
end

function Base.append!(
        traj::VariableVolumeTrajectory{D,LU,TP,TB},
        pos::AbstractVector{SVector{D, TP}},
        vel::AbstractVector{SVector{D, TV}},
        cell::PeriodicCell{D, TP}
    ) where {D, LU, TP<:Unitful.Length, TV<:Unitful.Velocity, TB<:SimpleVelocityTrajectory}
    append!(traj.base_trajectory, pos, vel)
    push!(traj.cell, cell)
    return traj
end

function Base.append!(
        traj::VariableVolumeTrajectory{D,LU,TP,TB},
        pos::AbstractVector{SVector{D, TP}},
        cell::PeriodicCell{D, TP}
    ) where {D, LU, TP, TB<:SimpleTrajectory}
    append!(traj.base_trajectory, pos)
    push!(traj.cell, cell)
    return traj
end


function Base.eltype(traj::VariableVolumeTrajectory{D, LU, TP, TB}) where {D, LU, TP, TB}
    return AtomsSystems.CellSystemView{D, LU, eltype(traj.base_trajectory), PeriodicCell{D, TP}}   
end


Base.show(io::IO, trj::VariableVolumeTrajectory) =
    print(io, "VariableVolumeTrajectory with ", length(trj), " frames, each with ", n_atoms(trj), " atoms")


##

"""
    ConstantVolumeTrajectory{D, LU, TP, TB} <: AbstractCellTrajectory{D, LU, TP}

A trajectory with constant cell. Use `VariableVolumeTrajectory` for trajectories with changing cell.

This should be used only if you know that the cell does not change during the trajectory.
It is a bit more efficient than `VariableVolumeTrajectory` because it only stores once cell for the whole trajectory.

# Type Parameters
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit of length for positions.
- `TP`: Type of positions `Unitful.Length`.
- `TB`: Type of the base trajectory, which is a `SimpleVelocityTrajectory` or `SimpleTrajectory`.

# Fields
- `base_trajectory`: The base trajectory, which is either a `SimpleVelocityTrajectory` or a `SimpleTrajectory`.
- `cell`: The periodic cell that defines the boundaries and periodicity for the trajectory.

# Example
```julia
traj = ConstantVolumeTrajectory(first_frame)
push!(traj, second_frame) 

# Chemfiles extension
using Chemfiles
traj = ConstantVolumeTrajectory("path/to/trajectory")
```
"""
mutable struct ConstantVolumeTrajectory{D, LU, TP, TB} <: AbstractCellTrajectory{D, LU, TP}
    base_trajectory::TB
    cell::PeriodicCell{D, TP}
    function ConstantVolumeTrajectory(
        traj::AbstractSimpleTrajectory{D, LU, TP},
        cell::PeriodicCell{D, TP}
    ) where {D, LU, TP}
        new{D, LU, TP, typeof(traj)}(
            traj,
            cell
        )
    end
end

ConstantVolumeTrajectory(traj::AbstractSimpleTrajectory{D}, ::IsolatedCell{D}) where {D} = traj

function ConstantVolumeTrajectory(sys)
    base = SimpleVelocityTrajectory(sys)
    if ! (cell(sys) isa PeriodicCell)
        return base
    end
    return ConstantVolumeTrajectory(base, cell(sys))
end

function ConstantVolumeTrajectory(traj::AbstractVector{<:AbstractSystem})
    frame = traj[begin]
    tmp = ConstantVolumeTrajectory(
        SimpleVelocityTrajectory(traj),
        cell(frame)
    )
    return tmp
end

function Base.getindex(sys::ConstantVolumeTrajectory, i::Int)
    @argcheck i >= 1 && i <= length(sys) BoundsError(sys, i)
    return AtomsSystems.CellSystemView(
        sys.base_trajectory[i],
        sys.cell
    )   
end


function Base.push!(traj::ConstantVolumeTrajectory{D, LU, TP}, sys::AbstractSystem{D}) where {D, LU, TP}
    @argcheck length(sys) == n_atoms(traj) "System must have the same number of atoms as the trajectory"
    if cell(sys) == traj.cell
        push!(traj.base_trajectory, sys)
    else
        error("System cell must match the trajectory cell to push") 
    end
    return traj
end

function Base.eltype(traj::ConstantVolumeTrajectory{D, LU, TP, TB}) where {D, LU, TP, TB}
    return AtomsSystems.CellSystemView{D, LU, eltype(traj.base_trajectory), PeriodicCell{D, TP}}   
end

Base.show(io::IO, trj::ConstantVolumeTrajectory) =
    print(io, "ConstantVolumeTrajectory with ", length(trj), " frames, each with ", n_atoms(trj), " atoms")