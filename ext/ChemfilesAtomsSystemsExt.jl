module ChemfilesAtomsSystemsExt

using AtomsBase
using AtomsSystems
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

end