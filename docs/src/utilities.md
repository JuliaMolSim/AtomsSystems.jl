# Quality of Life Extensions to AtomsBase

**Methods to change positions**

```julia
center_of_mass(sys)

# translate the whole system by r
translate_system!(sys, r)

# translate a copy of system by r
translate_system(sys, r)

# rotate the system
using Rotations
rot = rand(RotMatrix{3})
rotate_system!(sys, rot)

# rotate a copy of the system
rotate_system(sys, rot)
```

**Add system together or repeat them**

```julia
# make system that has sys2 added to sys1, keep sys1 and sys2 as they are
add_systems(sys1, sys2)

# repeat system along cell vectors
# repeat system 3 times along all cell vectors
repeat(sys, 3)

# repeat 2 times on the first cell vector, 3 times on th esecond cell vector
# and 4 times along the third cell vector
repeat(sys, (2,3,4))
```

**Methods to get information from systems**

```julia
# distance of atoms i and j as a vector
distance_vector(sys, i , j)

# distance of atoms i and j
distance(sys, i, j)

# bond angle of atom i, j and k (j->i vs j->k)
bond_angle(sys, i, j, k)

# dihedral angle of atoms i, j, k and m
dihedral_angle(sys, i, j, k, m)
```

## Fractional Coordinate Methods

```julia
# get inverse cell of the system as matrix
inv_cell(sys)

# get cell_vectors as matrix
cell_matrix(sys)

# fractional coordinates of atom(s) i
fractional_coordinates(sys, i)

# fractional coordinates as matrix for atom(s) i
fractional_coordinates_as_matrix(sys, i)

# wrap atoms inside the cell
wrap_coordinates!(sys)
```