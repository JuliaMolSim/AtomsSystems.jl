module AtomsTrajectories

using ArgCheck
using AtomsBase
using LinearAlgebra: norm
using StaticArrays
using Unitful
import ..AtomsSystems

export SimpleTrajectory
export SimpleVelocityTrajectory
export ConstantVolumeTrajectory
export VariableVolumeTrajectory
export AbstractTrajectory
export AbstractSimpleTrajectory


include("abstract_trajectory.jl")
include("simple_trajectory.jl")
include("cell_trajectory.jl")

end
