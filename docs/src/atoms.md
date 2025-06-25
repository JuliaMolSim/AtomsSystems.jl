# Atom Structure

New atom structure `SimpleAtom` is a bits type structure that aims to be a simple small structure.

You can create `SimpleAtom` in one of the following ways

```@repl
SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å")
SimpleAtom( 1, [0.0, 0.0, 0.0]u"Å") # same as above
SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å", [0.1, 0.0, 0.0]u"Å/s"; mass = 16.0u"u", charge = -1.0u"q")
SimpleAtom(ChemicalSpecies(:H), [0.0, 0.0, 0.0]u"Å")
SimpleAtom( :O => [1.0, 0.0, 0.0]u"Å" )
```

You can add extra atomkeys to an existing atom, by creating a new `SimpleAtom` and adding a keyword argument

```@example
sa = SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å")
SimpleAtom(sa; charge = 1.0u"q" ) # same as sa but with added charge
```


Compared to AtomsBase `Atom` [`SimpleAtom`](@ref) is smaller in size and is also bitstype 

```@repl
ab_atom = AtomsBase.Atom( :O, [1.0, 0.0, 0.0]u"Å" )
sa = SimpleAtom( :O, [1.0, 0.0, 0.0]u"Å" )
Base.summarysize(ab_atom)
Base.summarysize(sa)
isbits(ab_atom)
isbits(sa)
```


## Vector of SimpleAtoms

Vector os SimpleAtoms has basic AtomsBase interface implemented

```@repl
va = [ SimpleAtom(i, i * ones(3)u"Å") for i in 1:5 ]

species(va, 3)

position(va, 4)

mass(va, :)

hasatomkey(va, :velocity)

cell(va)
```