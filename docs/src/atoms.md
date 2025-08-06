# Atom Structure

```@setup atom
using AtomsBase
using AtomsSystems
using Unitful
```

New atom structure `SimpleAtom` is a bits type structure that aims to be a simple small structure.

The main advantage of `SimpleAtom` over the default atom type of AtomsBase is that it is bitstype.
The second advantage is that it does not have `Dict` to hold custom data, thus making it a lot more smaller data structure.  


## Create SimpleAtoms

You can create `SimpleAtom` in one of the following ways


```@repl atom
SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å")
SimpleAtom( 1, [0.0, 0.0, 0.0]u"Å") # same as above
SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å", [0.1, 0.0, 0.0]u"Å/s"; mass = 16.0u"u", charge = -1.0u"q")
SimpleAtom(ChemicalSpecies(:H), [0.0, 0.0, 0.0]u"Å")
SimpleAtom( :O => [1.0, 0.0, 0.0]u"Å" )
```

You can add extra atomkeys to an existing atom, by creating a new `SimpleAtom` and adding a keyword argument

```@repel atom
sa = SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å")
SimpleAtom(sa; charge = 1.0u"q" ) # add charge
```


## Vector of SimpleAtoms

Vector os SimpleAtoms has basic AtomsBase interface implemented

```@repl atom
va = [ SimpleAtom(i, i * ones(3)u"Å") for i in 1:5 ]

species(va, 3)

position(va, 4)

mass(va, :)

hasatomkey(va, :velocity)

cell(va)
```