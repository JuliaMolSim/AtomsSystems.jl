# Chemfiles Extension

AtomsSystems - Chemfiles extension allows you to load files with [Chemfiles](https://github.com/chemfiles/Chemfiles.jl).

This is mainly intended to be implemented in [AtomsIO](https://github.com/mfherbst/AtomsIO.jl), but the interface is there, if you want to use it.


```julia
using AtomsSystems
using Chemfiles


# read a frame
sys = Trajectory("example.xyz") do trajectory
    frame = read(trajectory)
    CellSystem(frame)
end


# Read the whole trajectory
# Note this can be slow!
traj = Trajectory("example.xyz") do trajectory
    map(trajectory) do frame
        SimpleSystem(frame)
        # or SimpleVelocitySystem(frame)
        # or CellSystem(frame)
    end
end
```