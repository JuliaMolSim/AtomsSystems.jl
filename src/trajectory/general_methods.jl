

function AtomsSystems.distance(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, frame::Int=1)
    return AtomsSystems.distance(traj[frame], i, j)
end

function AtomsSystems.distance(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, ::Colon)
    map(traj) do sys
        AtomsSystems.distance(sys, i, j)
    end
end

function AtomsSystems.distance(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, frames)
    map( frames ) do frame
        AtomsSystems.distance(traj[frame], i, j)
    end
end

function AtomsSystems.dihedral_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, l::Int, frame::Int)
    return AtomsSystems.dihedral_angle(traj[frame], i, j, k, l)
end

function AtomsSystems.dihedral_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, l::Int, ::Colon)
    map(traj) do sys
        AtomsSystems.dihedral_angle(sys, i, j, k, l)
    end
end

function AtomsSystems.dihedral_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, l::Int, frames)
    map( frames ) do frame
        AtomsSystems.dihedral_angle(traj[frame], i, j, k, l)
    end
end


function AtomsSystems.bond_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, frame::Int)
    return AtomsSystems.bond_angle(traj[frame], i, j, k)
end

function AtomsSystems.bond_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, ::Colon)
    map(traj) do sys
        AtomsSystems.bond_angle(sys, i, j, k)
    end     
end

function AtomsSystems.bond_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, frames)
    map( frames ) do frame
        AtomsSystems.bond_angle(traj[frame], i, j, k)
    end       
end

function AtomsSystems.distance_vector(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, frame::Int)
    return AtomsSystems.distance_vector(traj[frame], i, j)
end

function AtomsSystems.distance_vector(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, ::Colon)
    map(traj) do sys
        AtomsSystems.distance_vector(sys, i, j)
    end
end

function AtomsSystems.distance_vector(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, frames)
    map( frames ) do frame
        AtomsSystems.distance_vector(traj[frame], i, j)
    end
end