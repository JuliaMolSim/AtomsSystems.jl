"""
    GeneralSystem{D, LU, TB} <: AbstractCompositeSystem{D, LU}

This is the most generic system type that can hold any system properties.

To build this use always the `generic_system` function. This type is meant to be used
only on the background.
"""
mutable struct GeneralSystem{D, LU, TB} <: AbstractCompositeSystem{D, LU}
    base_system::TB
    properties::Dict{Symbol, Any}
    function GeneralSystem(
        sys::Union{CellSystem{D,LU}, AbstractIsolatedSystem{D,LU}};
        kwargs...
    ) where {D,LU}
        props = Dict{Symbol, Any}( k=>v for (k,v) in pairs(kwargs) if ! in(k, (:cell, :cell_vectors, :periodicity)) )
        if length(props) == 0 
            return sys
        end
        new{D, LU, typeof(sys)}(sys, props)
    end
end

GeneralSystem(sys::GeneralSystem) = deepcopy(sys)


function Base.getindex(sys::GeneralSystem, x::Symbol)
    if x in (:cell_vectors, :periodicity)
        return sys.base_system[x]
    else
        return sys.properties[x]
    end
end

function Base.keys(sys::GeneralSystem)
    tmp = keys(sys.properties)
    return (:cell_vectors, :periodicity, tmp...)
end

Base.haskey(sys::GeneralSystem, x::Symbol) = in(x, keys(sys))


function AtomsBase.set_cell_vectors!(sys::GeneralSystem, x)
    sys = AtomsBase.set_cell_vectors!(sys.base_system, x)
    return sys
end

function AtomsBase.set_periodicity!(sys::GeneralSystem, x)
    sys = AtomsBase.set_periodicity!(sys.base_system, x)
    return sys
end

function AtomsBase.set_cell!(sys::GeneralSystem, x)
    sys = AtomsBase.set_cell!(sys.base_system, x)
    return sys
end

function AtomsBase.set_cell!(sys::GeneralSystem, x::IsolatedCell)
    tmp = AtomsBase.set_cell!(sys.base_system, x)
    return GeneralSystem(tmp; sys.properties...)
end


## Generic builders

function _parse_cell_vectors(cv::Vararg{AbstractVector{<:Unitful.Length}})
    # make sure we have a square matrix
    tmp = reduce(hcat, cv)
    if isa(tmp, SMatrix)
        return Tuple( x for x in eachcol(tmp) )
    end
    if size(tmp, 1) == size(tmp, 2)
        D = size(tmp, 1)
        tmp = SMatrix{D,D}( tmp )
        return Tuple( x for x in eachcol(tmp) )
    end
    throw(ArgumentError("Cell vectors must be square matrix"))
end

"""
    generic_system(sys::AbstractSystem; kwargs...)

Create a system from an existing system, allowing to update the cell or other system wide properties.

The new system is a copy of the original system, it does not share any data with the original system.

# Examples
```julia
# form a come of a system
sys = generic_system(some_system)

# form a copy and add energy
sys = generic_system(some_system; energy = 10.0u"eV")

# Create a system with updated cell vectors and periodicity
sys = generic_system(some_system; cell_vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"Å", periodicity = (true, true, true))
my_cell = PeriodicCell([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"Å", (true, true, true))
sys = generic_system(some_system; cell = my_cell)

# Copy of a system with changed cell vectors
sys = generic_system(some_system; cell_vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"Å")

# Copy some_system and take cell from another system
sys = generic_system(some_system; cell = cell(another_system))
```
"""
function generic_system(sys::AbstractSystem; kwargs...)
    # Figure out do we need to update the cell
    new_cell = nothing
    if haskey(kwargs, :cell_vectors) && haskey(kwargs, :periodicity)
        cv = _parse_cell_vectors(kwargs[:cell_vectors]...)
        new_cell = PeriodicCell( cell_vectors=cv, periodicity=kwargs[:periodicity] )
    elseif haskey(kwargs, :cell_vectors) && isa(cell(sys), PeriodicCell)
        cv = _parse_cell_vectors(kwargs[:cell_vectors]...)
        new_cell = PeriodicCell( cell_vectors=cv, periodicity=periodicity(sys) )
    elseif haskey(kwargs, :periodicity) && isa(cell(sys), PeriodicCell)
        new_cell = PeriodicCell( cell_vectors=cell_vectors(sys), periodicity=kwargs[:periodicity] )
    elseif haskey(kwargs, :cell)
        new_cell = kwargs[:cell]
    elseif isa(cell(sys), IsolatedCell) &&
            ( haskey(kwargs, :cell_vectors) && !haskey(kwargs, :periodicity))
        cm = _parse_cell_vectors(kwargs[:cell_vectors]...)
        new_cell = PeriodicCell(cm, Tuple( true for _ in 1:length(cm) ) )
    elseif isa(cell(sys), IsolatedCell) &&
            ( haskey(kwargs, :periodicity) && !haskey(kwargs, :cell_vectors))
            throw(ArgumentError("Cannot set periodicity without cell vectors"))
    end

    if isnothing(new_cell)
        tsys = CellSystem(sys)
        if length(kwargs) == 0 && length(keys(sys)) == 2 # only :cell_vectors and :periodicity
            return tsys
        elseif length(kwargs) == 0 && length(keys(sys)) > 2
            tmp = NamedTuple( k=>sys[k] for k in keys(sys) )
            return GeneralSystem(tsys; tmp...)
        elseif length(kwargs) > 0 && length(keys(sys)) == 2
            return GeneralSystem(tsys; kwargs...)
        else
            tmp = Dict{Symbol, Any}( k=>sys[k] for k in keys(sys) )
            foreach( pairs(kwargs) ) do (k,v)
                tmp[k] = v
            end
            return GeneralSystem(tsys; tmp...)
        end
    end

    tsys = CellSystem( AtomicPropertySystem(sys), new_cell )
    if length(kwargs) > 0 && length(keys(sys)) == 2 # only :cell_vectors and :periodicity with old system
        return GeneralSystem(tsys; kwargs...)
    else
        tmp = Dict{Symbol, Any}( k=>sys[k] for k in keys(sys) )
        foreach( pairs(kwargs) ) do (k,v)
            tmp[k] = v
        end
        return GeneralSystem(tsys; tmp...)
    end
end

"""
    generic_system(vec_atoms::AbstractVector{<:AtomsBase.Atom}; kwargs...)

Create a system from an array of atoms, with optional keyword arguments to set system properties.

If cell is not specified, the new system will have an isolated cell.

# Examples
```julia
# Create a system from an array of atoms
sys = generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"), 
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]
)

# Create a system and set energy
sys = generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")];
    energy = 10.0u"eV"
)

# Create a system and set cell vectors and periodicity
sys = generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]; 
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"], 
    periodicity = (true, true, true)
)
```
"""
function generic_system(sys::AtomsVector; kwargs...)
    tmp = AtomicPropertySystem(sys)
    if length(kwargs) > 0
        return generic_system(tmp; kwargs...)
    end
    return tmp
end

function generic_system(sys::AbstractVector{<:Union{AtomsBase.Atom, AtomsBase.AtomView}}; kwargs...)
    tmp = SimpleAtom.(sys)
    return generic_system(tmp; kwargs...)
end

"""
    generic_system(species, r; kwargs...)
    generic_system(species, r, v; kwargs...)

Create a system from an array of species and positions, with optional keyword arguments to set system properties.

If cell is not specified, the new system will have an isolated cell.

# Examples
```julia
# Create a system from species and positions
sys = generic_system([:H, :O], [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"])

# Create a system and set energy
sys = generic_system([:H, :O], [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"]; energy = 10.0u"eV")

# Create a system and set cell vectors and periodicity
sys = generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"]; 
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"], 
    periodicity = (true, true, true)
)

# Create a system with velocities
sys = generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"],
    [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"];
    energy = 10.0u"eV"
)
```
"""
function generic_system(
    species::AbstractVector{<:AtomsBase.AtomId},
    r::AbstractVector{<:AbstractVector};
    kwargs...
)   
    spc = ChemicalSpecies.(species)
    tmp = SimpleSystem(spc, r)
    if length(kwargs) > 0
        return generic_system(tmp; kwargs...)
    end
    return tmp
end

 
function generic_system(
    species::AbstractVector{<:AtomsBase.AtomId}, 
    r::AbstractVector{<:AbstractVector},
    v::AbstractVector{<:AbstractVector}; 
    kwargs...
)   
    spc = ChemicalSpecies.(species)
    tmp = SimpleVelocitySystem(spc, r, v)
    if length(kwargs) > 0
        return generic_system(tmp; kwargs...)
    end
    return tmp
end

"""
    generic_system(prs::AbstractVector{<:Pair}; kwargs...)

Create a system from an array of pairs, where each pair consists of a species and a position.
The new system is a copy of the original system, allowing to set additional properties via keyword arguments.

If cell is not specified, the new system will have an isolated cell.

# Examples
```julia
# Create a system from an array of pairs
sys = generic_system([:H => [0.0, 0.0, 0.0]u"Å", :O => [1.0, 0.0, 0.0]u"Å"])

# Create a system and set energy
sys = generic_system([:H => [0.0, 0.0, 0.0]u"Å", :O => [1.0, 0.0, 0.0]u"Å"]; energy = 10.0u"eV")

# Create a system with cell vectors and periodicity
sys = generic_system([
    :H => [0.0, 0.0, 0.0]u"Å", 
    :O => [1.0, 0.0, 0.0]u"Å"];
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"], 
    periodicity = (true, true, true)
)
```
"""
function generic_system(prs::AbstractVector{<:Pair}; kwargs...)
    return generic_system( SimpleAtom.(prs); kwargs... )
end


"""
    generic_system(sys::AbstractSystem, i; kwargs...)
    generic_system(sys::AbstractSystem, vspc::ChemicalSpecies...; kwargs...)

Create a subsystem from an existing system and a selection of atoms, allowing to set additional properties via keyword arguments.

The new subsystem allocates its own data and does not share any data with the original system.

# Examples
```julia
# Create a subsystem with atoms 1:5 from an existing system
sys = generic_system(some_system, 1:5)

# Create a sybsys with specific species from an existing system
sys = generic_system(some_system, ChemicalSpecies(:H), ChemicalSpecies(:O))

# Create a subsystem and set energy
sys = generic_system(some_system, 1:5; energy = 10.0u"eV")
```
"""
function generic_system(sys::AbstractSystem, i; kwargs...)
    return generic_system(CellSystem(sys, i); kwargs...)
end

function generic_system(sys::AbstractSystem, vspc::ChemicalSpecies...; kwargs...)
    tmp = CellSystem(sys, vspc...)
    if length(kwargs) > 0
        return GeneralSystem(tmp; kwargs...)
    end
    return tmp
end

"""
    generic_system(sys::AbstractSystem{D}, v::AbstractVector{<:AbstractVector{<:Unitful.Velocity}}) where {D}

Create a system from an existing system and a vector of velocities.

This allows you to add velocities to a system while keeping the original system's properties.
"""
function generic_system(sys::AbstractSystem{D}, v::AbstractVector{<:AbstractVector{<:Unitful.Velocity}}) where {D}
    @argcheck length(v) == length(sys) "The length of the velocity vector must match the number of atoms in the system"
    @argcheck all( length.(v) .== D ) "All velocity vectors must have the same dimension as the system"
    
    tmp = AtomicPropertySystem(sys)
    vtmp = SimpleVelocitySystem(species(tmp, :), position(tmp, :), v)
    ntmp = _combine_helper(vtmp, tmp)
    if cell(sys) isa IsolatedCell
        # Don't pass cell vectors or periodicity
        t = filter( x -> !in(x, (:cell_vectors, :periodicity)), keys(sys) )
        if isempty(t)
            return ntmp
        end
        return generic_system(ntmp; (k => sys[k] for k in t)...)
    end
    return generic_system(ntmp; pairs(sys)...)
end

_combine_helper(bsys::AbstractSimpleSystem, sys2::AtomicPropertySystem) = AtomicPropertySystem( bsys, sys2.atom_properties )
_combine_helper(bsys::AbstractSimpleSystem, ::AbstractSimpleSystem) = bsys



function generic_system(atoms::SimpleAtom...; kwargs...)
    return generic_system(collect(atoms); kwargs...)
end


##

macro generic_system_str(input::AbstractString)
    atoms = []
    for line in split(input, '\n')
        line = strip(line)
        parts = split(line)
        if isempty(line) || line[1] == '#'
            continue  # Skip empty lines or comments
        end
        if length(parts) < 2
            throw(ArgumentError("Each line must contain at least a species and a position"))
        end
        spc = (ChemicalSpecies ∘ Symbol)(parts[1])
        pos = [ parse(Float64, word) for word in parts[2:end] ]
        push!(atoms, (species=spc, position=pos))
    end

    return quote
        generic_system([SimpleAtom(spc, pos...) for (spc, pos) in $(atoms)])
    end
end

