# System Structures

Under the hood there are several different system structures that differ on what data they store. But you only need to call `generic_system` to build all of them.
Depending on what information you provide you get different structure

```julia
# Build based on vector of atoms
# generic_system(atoms::AbstractVector{SimpeAtom}; kwargs...)
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"), 
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]
)

# Same but added key for energy
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")];
    energy = 10.0u"eV"
)

# Add cell to the system
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]; 
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"], 
    periodicity = (true, true, true)
)

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