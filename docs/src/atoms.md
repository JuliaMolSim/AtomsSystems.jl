# Atom Structure

New atom structure `SimpleAtom` is a bits type structure that aims to be a simple small structure.

You can create `SimpleAtom` in one of the following ways

```julia
SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å")
SimpleAtom( 1, [0.0, 0.0, 0.0]u"Å") # same as above
SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å", [0.1, 0.0, 0.0]u"Å/s"; mass = 16.0u"u", charge = -1.0u"q")
SimpleAtom(ChemicalSpecies(:H), [0.0, 0.0, 0.0]u"Å")
SimpleAtom( :O => [1.0, 0.0, 0.0]u"Å" )
```


Comparison to AtomsBase `Atom`
```julia-repl
julia> ab_atom = AtomsBase.Atom( :O, [1.0, 0.0, 0.0]u"Å" )
Atom(O, Z = 8, m = 15.999 u):
    position          : [1,0,0]u"Å"
    species           : O

julia> sa = SimpleAtom( :O, [1.0, 0.0, 0.0]u"Å" )
SimpleAtom(O, Z = 8, m = 15.999 u):
    position          : [1,0,0]u"Å"
    species           : O

julia> Base.summarysize(ab_atom)
456

julia> Base.summarysize(sa)
32

julia> isbits(ab_atom)
false

julia> isbits(sa)
true
```