
"""
    SimpleSystemView{D, LU, TP, TI, L}

View subsystem of a system that contains only species and position.

It is recommended to use `system_view` function to create a view of the system
instead of creating this type directly.

# Parameter Definitions
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit length of the system.
- `TP`: Type of position vector (e.g., `SVector{D, T}`).
- `TI`: Type index for SubArray used for the view.
- `L`:  True when SybArray uses linear indexing, false otherwise.

# Example
```julia
# crete view of a system with atoms 2 to 30
subsys = SimpleSystemView(sys, 2:30)
```
"""
mutable struct SimpleSystemView{D, LU, TP, TI, L} <: AbstractSimpleSystem{D, LU}
    species::SubArray{ChemicalSpecies, 1, Vector{ChemicalSpecies}, TI, L}
    position::SubArray{SVector{D,TP}, 1, Vector{SVector{D,TP}}, TI, L}
    function SimpleSystemView(sys::AbstractCompositeSystem{D, LU}, i) where {D, LU}
        vspc = view( species(sys, :), i)
        vpos = view( position(sys, :), i)
        T = typeof(vspc)
        TI = T.parameters[4]
        L  = T.parameters[5]
        TP = (eltype ∘ eltype)(vpos) 
        new{D, LU, TP, TI, L}(vspc, vpos)
    end
    function SimpleSystemView(
            spc::SubArray{ChemicalSpecies, 1, Vector{ChemicalSpecies}, TI, L},
            pos::SubArray{SVector{D,TP}, 1, Vector{SVector{D,TP}}, TI, L}
        ) where {D, TP, TI, L}
        @argcheck length(spc) == length(pos)
        LU = unit(TP)
        new{D, LU, TP, TI, L}(spc, pos)
    end
    function SimpleSystemView(sys::SimpleSystemView{D, LU, TP, TI, L}, i) where {D, LU, TP, TI, L}
        vspc = view(sys.species, i)
        vpos = view(sys.position, i)
        T = typeof(vspc)
        nTI = T.parameters[4]
        nL  = T.parameters[5]
        new{D, LU, TP, nTI, nL}(vspc, vpos)  
    end
end

"""
    SimpleVelocitySystemView{D, LU, TP, TV, TI, L}

View subsystem of a system that contains species, position, and velocity.

It is recommended to use `system_view` function to create a view of the system
instead of creating this type directly.

# Parameter Definitions
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit length of the system.
- `TP`: Type of position vector (e.g., `SVector{D, T}`).
- `TV`: Type of velocity vector (e.g., `SVector{D, T}`).
- `TI`: Type index for SubArray used for the view.
- `L`:  True when SybArray uses linear indexing, false otherwise.

# Example
```julia
# Create view of a system with atoms 2 to 30
subsys = SimpleVelocitySystemView(sys, 2:30)
````
"""
mutable struct SimpleVelocitySystemView{D, LU, TP, TV, TI, L} <: AbstractSimpleSystem{D, LU}
    species::SubArray{ChemicalSpecies, 1, Vector{ChemicalSpecies}, TI, L}
    position::SubArray{SVector{D,TP}, 1, Vector{SVector{D,TP}}, TI, L}
    velocity::SubArray{SVector{D,TV}, 1, Vector{SVector{D,TV}}, TI, L}
    function SimpleVelocitySystemView(sys::AbstractCompositeSystem{D, LU}, i) where {D, LU}
        vspc = view( species(sys, :), i)
        vpos = view( position(sys, :), i)
        vvel = view( velocity(sys, :), i)
        T = typeof(vspc)
        TI = T.parameters[4]
        L  = T.parameters[5]
        TP = (eltype ∘ eltype)(vpos)
        TV = (eltype ∘ eltype)(vvel)    
        new{D, LU, TP, TV, TI, L}(vspc, vpos, vvel)
    end
    function SimpleVelocitySystemView(
            spc::SubArray{ChemicalSpecies, 1, Vector{ChemicalSpecies}, TI, L},
            pos::SubArray{SVector{D,TP}, 1, Vector{SVector{D,TP}}, TI, L},
            vel::SubArray{SVector{D,TV}, 1, Vector{SVector{D,TV}}, TI, L}
        ) where {D, TP, TV, TI, L}
        @argcheck length(spc) == length(pos) == length(vel)
        LU = unit(TP)
        new{D, LU, TP, TV, TI, L}(spc, pos, vel) 
    end
    function SimpleVelocitySystemView(sys::SimpleVelocitySystemView{D, LU, TP, TV, TI, L}, i) where {D, LU, TP, TV, TI, L}
        vspc = view(sys.species, i)
        vpos = view(sys.position, i)
        vvel = view(sys.velocity, i)
        T = typeof(vspc)
        nTI = T.parameters[4]
        nL  = T.parameters[5]
        new{D, LU, TP, TV, nTI, nL}(vspc, vpos, vvel)
    end
end

"""
    AtomicPropertySystemView{D, LU, TB, TP, TI, SM, L}

View subsystem of a system that contains species, position, and additional atomic properties.

It is recommended to use `system_view` function to create a view of the system
instead of creating this type directly.

# Parameter Definitions
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit length of the system.
- `TB`: Type of the base system (e.g., `SimpleSystemView`).
- `TP`: Type of the atomic properties (e.g., `SVector{D, T}`).
- `TI`: Type index for SubArray used for the view.
- `SM`: Type of the atomic properties (e.g., `SVector{D, T}`).
- `L`:  True when SybArray uses linear indexing, false otherwise.
"""
mutable struct AtomicPropertySystemView{D, LU, TB, TP, TI, SM, L} <: AbstractIsolatedSystem{D, LU}
    base_system::TB
    atom_properties::SubArray{TP, 1, Vector{TP}, TI, L}
    function AtomicPropertySystemView(sys::AtomicPropertySystem{D, LU, TT, TP, SM}, i) where {D, LU, TT, TP, SM}
        base = system_view(sys.base_system, i)
        atom_properties = view(sys.atom_properties, i)
        T = typeof(atom_properties)
        TI = T.parameters[4]
        L = T.parameters[5]
        new{D, LU, typeof(base), TP, TI, SM, L}(base, atom_properties)
    end
    function AtomicPropertySystemView(
            sys::Union{SimpleSystemView{D, LU}, SimpleVelocitySystemView{D, LU}},
            atom_properties::SubArray{TP, 1, Vector{TP}, TI, L}
        ) where {D, LU, TP, TI, L}
        @argcheck length(sys) == length(atom_properties)
        SM = haskey(atom_properties[1], :mass) && atom_properties[1][:mass] isa Unitful.Mass
        new{D, LU, typeof(sys), TP, TI, SM, L}(sys, atom_properties)
    end
    function AtomicPropertySystemView(
            sys::AtomicPropertySystemView{D, LU, TB, TP, TIX, SM, LX},
            i
        ) where {D, LU, TB, TP, TIX, SM, LX}
        base = system_view(sys.base_system, i)
        atom_properties = view(sys.atom_properties, i)
        T = typeof(atom_properties)
        TI = T.parameters[4]
        L = T.parameters[5]
        new{D, LU, typeof(base), TP, TI, SM, L}(base, atom_properties)
    end
end

"""
    CellSystemView{D, LU, TB, TC}

View subsystem of a system that contains a base system and a cell.

It is recommended to use `system_view` function to create a view of the system
instead of creating this type directly.

# Parameter Definitions
- `D`: Dimension of the system (e.g., 2 or 3).
- `LU`: Unit length of the system.
- `TB`: Type of the base system (e.g., `SimpleSystemView`).
- `TC`: Type of the cell (has to be subtype of `PeriodicCell`).
"""
mutable struct CellSystemView{D, LU, TB, TC} <: AbstractCompositeSystem{D, LU}
    base_system::TB
    cell::TC
    function CellSystemView(sys::CellSystem{D, LU, TT, TC}, i) where {D, LU, TT, TC}
        base = system_view(sys.base_system, i)
        new{D, LU, typeof(base), TC}(base, sys.cell)
    end
    function CellSystemView(
        ssys::Union{SimpleSystemView{D, LU}, SimpleVelocitySystemView{D, LU}, AtomicPropertySystemView{D, LU}},
        cell::PeriodicCell{D}
    ) where {D, LU}
        new{D, LU, typeof(ssys), typeof(cell)}(ssys, cell)
    end
    function CellSystemView(sys::CellSystemView{D, LU, TB, TC}, i) where {D, LU, TB, TC}
        base = system_view(sys.base_system, i)
        new{D, LU, typeof(base), TC}(base, sys.cell)        
    end
end

##

Base.getindex(sys::SimpleSystemView, i::Int) = SimpleAtom( sys.species[i],  sys.position[i] )
Base.getindex(sys::SimpleVelocitySystemView, i::Int) = SimpleAtom( sys.species[i],  sys.position[i], sys.velocity[i] )

function Base.getindex(sys::AtomicPropertySystemView, i::Int)
    tmp = sys.base_system[i]
    SimpleAtom( (; tmp.data..., sys.atom_properties[i]...) )
end

Base.getindex(sys::AtomicPropertySystemView, x::Symbol) = sys.base_system[x]

function Base.getindex(sys::CellSystemView, x::Symbol)
    if x == :cell_vectors
        cell_vectors(cell(sys))
    elseif x == :periodicity
        return periodicity(cell(sys))
    else
        return sys.base_system[x]
    end
end

AtomsBase.atomkeys(::SimpleSystemView) = (:position, :species)
AtomsBase.atomkeys(::SimpleVelocitySystemView) = (:position, :species, :velocity)

AtomsBase.cell(::SimpleSystemView{D, LU, TP}) where {D, LU, TP} = IsolatedCell(D, TP)
AtomsBase.cell(::SimpleVelocitySystemView{D, LU, TP}) where {D, LU, TP} = IsolatedCell(D, TP)
AtomsBase.cell(sys::CellSystemView) = sys.cell

AtomsBase.velocity(sys::SimpleVelocitySystemView, i::Int)  = sys.velocity[i]
AtomsBase.velocity(sys::SimpleVelocitySystemView, i)       = view(sys.velocity, i)
AtomsBase.velocity(sys::SimpleVelocitySystemView, ::Colon) = sys.velocity

AtomsBase.set_velocity!(sys::SimpleVelocitySystemView, i, x) = setindex!(sys.velocity, x, i)




"""
    system_view(sys, i)
    system_view(sys, spc...)

Create a view of a system for the given indices or species.
    
This does not create a copy of the system, but rather a view of the existing system
and thus does not allocate memory (except for the view itself).

# Example
```julia
# Create a view of a system with atoms 2 to 30
subsys = system_view(sys, 2:30)

# Create a view of a system with species Si and Al
subsys = system_view(sys, ChemicalSpecies(:Si), ChemicalSpecies(:Al))
```
"""
system_view(sys::Union{SimpleSystem, SimpleSystemView}, i) = SimpleSystemView(sys, i)
system_view(sys::Union{SimpleVelocitySystem, SimpleVelocitySystemView}, i) = SimpleVelocitySystemView(sys, i)
system_view(sys::Union{AtomicPropertySystem, AtomicPropertySystemView}, i) = AtomicPropertySystemView(sys, i)
system_view(sys::Union{CellSystem, CellSystemView}, i) = CellSystemView(sys, i)

system_view(sys::GeneralSystem, i) = CellSystemView(sys.base_system, i)

system_view(sys::Union{SimpleSystem, SimpleSystemView}, spc::ChemicalSpecies...) = _system_view(sys, spc...)
system_view(sys::Union{SimpleVelocitySystem, SimpleVelocitySystemView}, spc::ChemicalSpecies...) = _system_view(sys, spc...)
system_view(sys::Union{AtomicPropertySystem, AtomicPropertySystemView}, spc::ChemicalSpecies...) = _system_view(sys, spc...)
system_view(sys::Union{CellSystem, CellSystemView}, spc::ChemicalSpecies...) = _system_view(sys, spc...) 

system_view(sys::GeneralSystem, spc::ChemicalSpecies...) = _system_view(sys, spc...)   
    
# Helper function to fight abiquity in the system_view function
function _system_view(sys::AbstractCompositeSystem, spc::ChemicalSpecies...)
    i = species(sys, :) .∈ Ref(spc)
    if count(i) == 0
       throw(ArgumentError("No species found in the system matching the provided species."))
    end
    return system_view(sys, i)
end