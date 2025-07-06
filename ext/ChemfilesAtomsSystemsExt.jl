module ChemfilesAtomsSystemsExt

using AtomsBase
using AtomsSystems
using AtomsSystems.AtomsTrajectories
using Chemfiles
using Unitful
using StaticArrays


function AtomsSystems.SimpleSystem(frame::Chemfiles.Frame)
    pos = positions(frame)
    r = reinterpret(reshape, SVector{3, Float64}, pos) * u"Å"
    spc = map( frame ) do a
        ChemicalSpecies( atomic_number(a) )
    end
    AtomsSystems.SimpleSystem(spc, r)
end

function AtomsSystems.SimpleVelocitySystem(frame::Chemfiles.Frame)
    pos = positions(frame)
    r = reinterpret(reshape, SVector{3, Float64}, pos) * u"Å"
    spc = map( frame ) do a
        ChemicalSpecies( atomic_number(a) )
    end
    if has_velocities(frame)
        vel = velocities(frame)
        v = reinterpret(reshape, SVector{3, Float64}, vel) * u"Å/ps"
        return AtomsSystems.SimpleVelocitySystem(spc, r, v)
    end
    return AtomsSystems.SimpleSystem(spc, r)
end

function _load_cell(frame::Chemfiles.Frame)
    ccell = Chemfiles.UnitCell(frame)
    cell_shape = Chemfiles.shape(ccell)
    if cell_shape == Chemfiles.Infinite
        return IsolatedCell(3)
    else
        cell_vectors = collect(eachrow(Chemfiles.matrix(ccell)))u"Å"
        return PeriodicCell(; cell_vectors, periodicity=(true, true, true))
    end 
end

function AtomsSystems.CellSystem(frame::Chemfiles.Frame)
    sys = AtomsSystems.SimpleVelocitySystem(frame)
    cell = _load_cell(frame)
    return AtomsSystems.CellSystem(sys, cell)
end


function AtomsSystems.AtomsTrajectories.SimpleTrajectory(traj::Chemfiles.Trajectory)
    # SimpleTrajectory has constant cell, so we can use the first frame to initialize the system
    first_frame = Chemfiles.read_step(traj, 0)
    sys = AtomsSystems.SimpleSystem(first_frame)
    ntraj = AtomsSystems.AtomsTrajectories.SimpleTrajectory(sys)
    for frame in Iterators.drop(traj,1)
        pos = positions(frame) * u"Å"
        r = reinterpret(reshape, SVector{3, Float64}, pos) * u"Å"
        append!(ntraj, r)
    end
    return ntraj 
end


function AtomsSystems.AtomsTrajectories.SimpleVelocityTrajectory(traj::Chemfiles.Trajectory)
    # SimpleVelocityTrajectory has constant cell, so we can use the first frame to initialize the system
    first_frame = Chemfiles.read_step(traj, 0)
    if has_velocities(first_frame)
        sys = AtomsSystems.SimpleVelocitySystem(first_frame)
        ntraj = AtomsSystems.AtomsTrajectories.SimpleVelocityTrajectory(sys)
        for frame in Iterators.drop(traj, 1)
            pos = positions(frame) * u"Å"
            vel = velocities(frame) * u"Å/ps"
            r = reinterpret(reshape, SVector{3, Float64}, pos) * u"Å"
            v = reinterpret(reshape, SVector{3, Float64}, vel) * u"Å/ps"
            append!(ntraj, r, v)
        end
        return ntraj
    end
    return AtomsSystems.AtomsTrajectories.SimpleTrajectory(traj)
end

function AtomsSystems.AtomsTrajectories.ConstantVolumeTrajectory(traj::Chemfiles.Trajectory)
    # ConstantVolumeTrajectory has constant cell, so we can use the first frame to initialize the system
    first_frame = Chemfiles.read_step(traj, 0)
    tcell = _load_cell(first_frame)
    btraj = AtomsSystems.AtomsTrajectories.SimpleVelocityTrajectory(traj)
    return AtomsSystems.AtomsTrajectories.ConstantVolumeTrajectory(btraj, tcell)
end

function AtomsSystems.AtomsTrajectories.SimpleTrajectory(fname::AbstractString)
    traj = Chemfiles.Trajectory(fname)
    return AtomsSystems.AtomsTrajectories.SimpleTrajectory(traj)
end

function AtomsSystems.AtomsTrajectories.SimpleVelocityTrajectory(fname::AbstractString)
    traj = Chemfiles.Trajectory(fname)
    return AtomsSystems.AtomsTrajectories.SimpleVelocityTrajectory(traj)    
end

function AtomsSystems.AtomsTrajectories.ConstantVolumeTrajectory(fname::AbstractString; species_from=nothing)
    ttraj = Chemfiles.Trajectory(fname)
    traj = AtomsSystems.AtomsTrajectories.ConstantVolumeTrajectory(ttraj)
    if isnothing(species_from)
        return traj
    end
    tmp_traj = Chemfiles.Trajectory(species_from)
    frame = Chemfiles.read(tmp_traj)
    tmp = AtomsSystems.SimpleSystem(frame)
    sys = traj[1]
    if length(tmp) == length(sys)
        AtomsBase.set_species!(sys, :, species(tmp, :))
        return traj
    end
    error("Species vector length does not match the number of atoms in the trajectory")
end


end