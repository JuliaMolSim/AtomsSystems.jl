# AtomsSystems.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaMolSim.github.io/AtomsSystems.jl/dev)
[![Build Status](https://github.com/JuliaMolSim/AtomsSystems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaMolSim/AtomsSystems.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JuliaMolSim/AtomsSystems.jl/graph/badge.svg?token=QPK831PYGJ)](https://codecov.io/gh/JuliaMolSim/AtomsSystems.jl)




## Get Started

This will get you going, but please look for the [documentation](https://JuliaMolSim.github.io/AtomsSystems.jl/dev) and docstrings for more information.

```julia
using AtomsBase
using AtomsSystems
using Unitful

# Make a water molecule
sys = generic_system"""
    O     -2.1   0.6    0.0
    H     -1.4   0.4    0.6
    H     -1.8   1.3   -0.6
"""


# Add cell to the system (by default periodic)
sys = generic_system(sys; cell_vectors = [
    [5., 0., 0.]u"Å",
    [0., 5., 0.]u"Å",
    [0., 0., 5.]u"Å"],
)

# Set periodicity conditions
sys = generic_system(sys; periodicity=(false, true, false))

# Add additional data
sys = generic_system(sys; energy=1.2u"eV", label="my water")
```

### Build Systems with Atomic Features

```julia
# First build a vector of atoms
# Syntax is SimpleAtom(species, pos, [vel]; kwords...)
# or SimpleAtom(species, x, y, z; kwords...) e.g. SimpleAtom(:0, -2.1, 0.6, 0.0)
atoms = [
    SimpleAtom(:O, [-2.1, 0.6, 0.0]u"Å"; mass=16.5u"u" )
    SimpleAtom(:H, [-1.4, 0.4, 0.6]u"Å"; mass=2.3u"u"  )
    SimpleAtom(:H, [-1.8, 1.3, -0.6]u"Å"; mass=3.3u"u" )
]

# You can add more features for existing atoms
atoms = [
    SimpleAtom( atoms[1]; velocity=[0.1, 0.2, 0.0]u"Å/fs" , charge=-1.0u"q" )
    SimpleAtom( atoms[2]; velocity=[-0.2, 0.0, 0.1]u"Å/fs", charge=1.0u"q"  )
    SimpleAtom( atoms[3]; velocity=[0.0, -0.1, 0.2]u"Å/fs", charge=0.0u"q"  )
]

# Use vector of atoms to build a system
sys = generic_system(atoms)
```
