# Trajectory Submodule

!!! note ""
    Trajectory submodule is still work in progress. While it works, there is a lot of room to improve.

Experimental trajectory submodule adds support for trajectory structures that are subtypes of `AbstractVector{AbstractSystem}`.
These are more optimized than simply storing systems in a vector.

The optimizations come from the fact that data structures are only allocated once for the trajectory not for each frame.


## Trajectory types

The basic principle is same as with systems. The structures are layered with additional properties.

- [`SimpleTrajectory`](@ref) - species and position only
- [`SimpleVelocityTrajectory`](@ref) - species, position and velocity
- [`ConstantVolumeTrajectory`](@ref) - constant cell trajectory
- [`VariableVolumeTrajectory`](@ref) - cell can change for each frame

All trajectory types expect that number of atoms and species are constants.




## Examples

You can create trajectories by giving vector of systems as and input

```@example traj
using AtomsBase
using AtomsSystems
using AtomsSystems.AtomsTrajectories

first_frame = generic_system"""
    H 0 0 0
    O 1 0 0
"""

second_frame = generic_system"""
    H 0   0 0
    O 1.1 0 0
"""

traj = VariableVolumeTrajectory([first_frame, second_frame])
```

or by creating trajectory by giving first frame and then pushing to the system 

```@example traj
traj = VariableVolumeTrajectory(first_frame)

push!(traj, second_frame)
```

To use you can just index it

```@repel traj
traj[2]
```