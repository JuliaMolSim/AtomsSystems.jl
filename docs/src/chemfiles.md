# Chemfiles Extension

AtomsSystems - Chemfiles extension allows you to load files with [Chemfiles](https://github.com/chemfiles/Chemfiles.jl).

This is mainly intended to be implemented in [AtomsIO](https://github.com/mfherbst/AtomsIO.jl), but the interface is there, if you want to use it.


```julia
using AtomsSystems
using Chemfiles


# read a frame
sys = Trajectory("example.xyz") do trajectory
    frame = read(trajectory, 0) # adjust to read a different frame
    CellSystem(frame)
end
```

## Trajectories with Chemfiles

Trajectories have an extension that allows them to be loaded directly from a file

```julia
using AtomsSystems
using AtomsSystems.AtomsTrajectories
using Chemfiles

# Works with all Chemfiles supported trajectories
traj = VariableVolumeTrajectory("trajectory file")

# If you know the trajectory has constant volume
traj = ConstantVolumeTrajectory("trajectory file") 
```

You can also give keyword argument `species_from` that allows loading species information form a different file

```julia
traj = VariableVolumeTrajectory("trajectory file"; species_from="file with species information")
```