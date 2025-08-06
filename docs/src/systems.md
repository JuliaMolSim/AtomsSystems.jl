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
# generic_system(species, position; kwargs...)
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
# generic_system(sys; kwargs...) - style
sys = generic_system(sys; energy = 10.0u"eV", label  = "my sys" )
```

To add velocity you can do some of the following

```@example system

# using vectors - generic_system(spc, pos, [vel]; kwargs...)
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

# or add velocity to existing system
sys = generic_system"""
    H 0 0 0
    O 1 0 0
"""
# using generic_system(sys, vel) call
sys = generic_system(sys, [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"])
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

You can also add a cell to an existing system by keyword arguments

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




## Create a Subsystem from an Existing System

You can give selection parameters for `generic_system` command to create subsystems. The created systems are allocated, so they **do not share** common data with the original system.

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

`system_view` allows you to create subsystems that **share** all data with the original system.


```julia
# Subsystem with atoms 1:5
# system_view(sys, subsys_definition...)
syb_sys = system_view(sys, 1:5)

# Subsystem with all H and O atoms from sys
sub_sys = system_view(sys, ChemicalSpecies(:H), ChemicalSpecies(:O))
```

Any changes you make to `system_view` structures is made to the host system and vise versa.

Note that `system_view` does not see the global features of the host system.


## Changing Systems

AtomsBase defines funtions to modify structures, the following list is supported:

- `set_position!(system, i, x)` - all structures
- `set_velocity!(system, i, v)` - all structures that have velocity
- `set_species!(system, i, spc)` - all structures
- `set_cell!(system, cell)` - only for structures with `PeriodicCell` and to another `PeriodicCell` with same dimension. System view structures do not support cell update.
- `set_cell_vectors!(system, bb)` - same as for `set_cell`
- `set_periodicity!(cell, pbc)` - same as for `set_cell`
- `append!(system1, system2)` - if systems have same information fields (e.g. both have velocity), same cell and dont have global features.
- `deleteat!(system, i)` - all systems excluding system views.

**Example**

```@example system
sys = generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"],
    [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"];
)

AtomsBase.set_position!(sys, 2, [0.0, 1.0, 0.0]u"Å")
AtomsBase.set_velocity!(sys, 2, [0.0, 0.2, 0.0]u"Å/s")
AtomsBase.set_species!(sys, 1, ChemicalSpecies(:N))
```

## System Structures Explained

There are following abstract types defined

- `AbstractCompositeSystem{D, LU}` - supertype for all structures and a subtype of `AbstractSystem{D}`
- `AbstractIsolatedSystem{D, LU}` - supertype for system with isolated cell (no cell information)
- `AbstractSimpleSystem{D, LU}` - lowest level of systems that store data for species, position and velocity

### SimpleSystems

The simplest types that are subtypes of `AbstractSimpleSystem` are

- [`SimpleSystem`](@ref) - holds only species and positions
- [`SimpleVelocitySystem`](@ref) - holds species, position and velocity

Every system structure is either directly one of these or holds a one. 

### CompositeSystems

Composite systems are system structures that hold a system, called `base_system`, and some additional information to that system.

For example [`CellSystem`](@ref) has following definition

```julia
mutable struct CellSystem{D, LU, TB, TC} <: AbstractCompositeSystem{D, LU}
    base_system::TB
    cell::TC
    function CellSystem(sys::AbstractIsolatedSystem{D, LU}, cell::PeriodicCell{D,T}) where {D, LU, T}
        new{D, LU, typeof(sys), typeof(cell)}(sys, cell)
    end
end
```

that adds `PeriodicCell` to a system that does not have cell.

[`AtomicPropertySystem`](@ref) is another composite system. It takes in a `AbstractSimpleSystem` and adds atomic properties to it. That is properties like charge, custom mass, etc.

Finally there is `GeneralSystem` that adds global features. Note that `GeneralSystem` should never be called directly. Call instead `generic_system`.


## Using Different Systems Directly Control What Information is Stored

You can use different systems to control what information is stored. This allows you to drop features you don't need.

Lets build a system what has a lots of features using `AtomsBaseTesting`

```@example system_example
using AtomsBase
using AtomsBaseTesting
using AtomsSystems
using Unitful

ref = make_test_system()
sys = generic_system(ref.system)
```

As you can see it has a lot of features. If you only need species and positions you can use [`SimpleSystem`](@ref)

```@example system_example
SimpleSystem(sys)
```

Similarly using [`SimpleVelocitySystem`](@ref) allows limiting to species, position and velocities

```@example system_example
SimpleVelocitySystem(sys)
```

[`AtomicPropertySystem`](@ref) preserves all atomic properties, but ignores cell and global features

```@example system_example
AtomicPropertySystem(sys)
```

[`CellSystem`](@ref) drops global features, but rest is preserved

```@example system_example
CellSystem(sys)
```