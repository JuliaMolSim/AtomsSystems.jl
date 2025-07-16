abstract type AbstractTrajectory{D, LU} <: AbstractVector{AbstractSystem{D}} end
abstract type AbstractSimpleTrajectory{D, LU, TP} <: AbstractTrajectory{D, LU} end
abstract type AbstractCellTrajectory{D, LU, TP} <: AbstractTrajectory{D, LU} end


Base.show(io::IO, ::MIME"text/plain", trj::AbstractTrajectory) = show(io, trj)

n_atoms(trj::AbstractTrajectory) = length(trj[1])
@inline n_atoms(trj::AbstractSimpleTrajectory) = length(trj.species)
n_atoms(trj::AbstractCellTrajectory) = n_atoms(trj.base_trajectory)