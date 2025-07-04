# System Structures

```@setup system
using AtomsBase
using AtomsSystems
using Unitful
```

Under the hood there are several different system structures that differ on what data they store. You only need to call `generic_system` to build all of them.
The basic principle is that the system is optimized to carry only the information you give.

The most simple system is build by giving only species and positions

```@example system
# Use macro to build system with xyz-file style syntax
generic_system"""
    H 0 0 0
    O 1 0 0
"""

# Build based on vector of species and positions
# generic_system(species::AbstractVector, position::AbstractVector{AbstractVector}; kwargs...)
generic_system(
    [:H, :O],
    [ [0.0, 0.0, 0.0]u"Å",  [1.0, 0.0, 0.0]u"Å" ]
)

# Save with vector of atoms
# generic_system(atoms::AbstractVector{SimpeAtom}; kwargs...)
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"), 
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]
)
```

Keyword arguments can be used to add global features 


```@example system
# Add global features at build
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")];
    energy = 10.0u"eV",
    label  = "my sys"
)

# or add global features to an existing system
sys = generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"), 
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]
)

sys = generic_system(sys; energy = 10.0u"eV", label  = "my sys" )
```

To add velocity you can do some of the following

```@example system

# using vectors - geric_system(spc, pos, [vel]; kwargs...)
generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"],
    [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"];
)

# or with atoms
generic_system(
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å", [0.1, 0.0, 0.0]u"Å/s"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å", [0.2, 0.0, 0.0]u"Å/s")
)

```

## Add Cell to the System

`generic_system` can alter cell, if given as a keyword argument. You can give only `cell_vectors`, in this case `periodicity` is set to `(true,true,true)`. 

```@example system
# Add cell to the system by giving only cell vectors
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]; 
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"], 
    #periodicity = (true, true, true) this comes by default, but you can edit it
)
```

You can also add/edit the cell of an existing system by keyword arguments

```@example system
sys = generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"), 
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]
)

sys = generic_system(
    sys;
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"],
)
```


```

# Create a system from an array of pairs
# generic_system(AbstractVector(<:Pair); kwargs...)
sys = generic_system([:H => [0.0, 0.0, 0.0]u"Å", :O => [1.0, 0.0, 0.0]u"Å"])

# Create a system vectors of atom symbols, positions and velocities
# geric_system(spc, pos, [vel]; kwargs...)
sys = generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"],
    [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"];
)
```

## Build System from Other Systems

You can build system from other system and modify system global features

```julia
# Form a copy of old system
# gneric_system(old_sys; kwargs...)
new_sys = generic_system(old_sys)

# Copy system and add a global feature
new_sys = generic_system(old_sys; energy=10.0u"eV")

# Copy system and add/change cell
new_sys = generic_system(
    old_sys; 
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"],
    periodicity = (true, true, true)
)
```

## Create a Subsystem from an Existing System

The created subsystem does not share any data with the system it was build

```julia
# Subsystem with atoms 1:5
# generic_system(sys, subsys_definition...; kwargs...)
sub_sys = generic_system(sys, 1:5)

# Subsystem with all H and O atoms from sys
sub_sys = generic_system(sys, ChemicalSpecies(:H), ChemicalSpecies(:O))

# Add global feature to subsystem
sub_sys = generic_system(sys, 1:5; label="the first 5 atoms")
```

## Subsystem Views

You can create subsystems that share all data with the host system by calling `system_view`

```julia
# Subsystem with atoms 1:5
# system_view(sys, subsys_definition...)
syb_sys = system_view(sys, 1:5)

# Subsystem with all H and O atoms from sys
sub_sys = system_view(sys, ChemicalSpecies(:H), ChemicalSpecies(:O))
```

Any changes you make to `system_view` structures is made to host system and vise versa.

Note that `system_view` does not see the global features of the host system.


## Changing Systems

AtomsBase defines funtions to modify structures, the following list is supported

- `set_position!(system, i, x)` - all structures
- `set_velocity!(system, i, v)` - all structures that have velocity
- `set_species!(system, i, spc)` - all structures
- `set_cell!(system, cell)` - only for structures with `PeriodicCell` and to another `PeriodicCell` with same dimension. System view structures do not support cell update.
- `set_cell_vectors!(system, bb)` - same as for `set_cell`
- `set_periodicity!(cell, pbc)` - same as for `set_cell`
- `append!(system1, system2)` - if systems have same information fields (e.g. both have velocity), same cell and dont have global features.

**Example**

```julia
sys = generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"],
    [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"];
)

AtomsBase.set_position!(sys, 2, [0.0, 1.0, 0.0]u"Å")
AtomsBase.set_velocity!(sys, 2, [0.0, 0.2, 0.0]u"Å/s")
AtomsBase.set_species!(sys, 1, ChemicalSpecies(:N))
```