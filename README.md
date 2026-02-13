# QuasiCrystal.jl

[![docs: dev](https://img.shields.io/badge/docs-stable-blue.svg)](https://codes.sota-shimozono/QuasiCrystal.jl/stable/)
[![Julia](https://img.shields.io/badge/julia-v1.10+-9558b2.svg)](https://julialang.org)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![codecov](https://codecov.io/gh/sotashimozono/Lattices.jl/graph/badge.svg?token=6E7VZ9MJMK)](https://codecov.io/gh/sotashimozono/Lattices.jl)
[![Build Status](https://github.com/sotashimozono/QuasiCrystal.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sotashimozono/QuasiCrystal.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

In this module, QuasiCrystals structures are available including,

- Fibonacci Lattice
- Penrose Tile
- Amman Beenker Lattice

## Features

### Unified AbstractLattice Interface

QuasiCrystal.jl provides a unified interface (`AbstractLattice{D}`) that allows both periodic lattices and aperiodic quasicrystals to be accessed consistently. This design enables applications like Lattice2DMonteCarlo to work uniformly with either structure type through consistent access to:

- **Sites/positions**: Access vertex positions in the structure
- **Bonds/connections**: Query nearest-neighbor relationships
- **Indexing**: Uniform site and bond indexing

### Key Types

- `AbstractLattice{D}`: Base abstract type for all lattice-like structures in D dimensions
- `AbstractQuasicrystal{D}`: Abstract type for quasicrystals (inherits from AbstractLattice)
- `AbstractTopology{D}`: Abstract type for periodic lattice topologies
- `UnitCell{D,T}`: Unit cell definition for periodic lattices
- `Lattice{Topology,T,B,I}`: Concrete type for periodic lattices with unit cell structure
- `QuasicrystalData{D,T,TileType}`: Data structure for generated quasicrystal patterns
- `Bond`: Represents connections between sites
- `Connection`: Connection rules for unit cell-based lattices

### Common Interface Methods

```julia
# Works for both quasicrystals and periodic lattices
get_positions(lattice)           # Get site positions
get_bonds(lattice)               # Get bond list
get_nearest_neighbors(lattice)   # Get neighbor indices for each site
num_sites(lattice)               # Total number of sites
num_bonds(lattice)               # Total number of bonds

# Build nearest-neighbor connectivity
build_nearest_neighbor_bonds!(qc_data; cutoff=2.0)
```

## Example

### Quasicrystal Generation and Analysis

```julia
using QuasiCrystal

# Generate a Fibonacci lattice (1D quasicrystal)
fibonacci_qc = generate_fibonacci_projection(20)

# Build nearest neighbor bonds
build_nearest_neighbor_bonds!(fibonacci_qc; cutoff=2.0)

# Access properties via unified interface
println("Number of sites: ", num_sites(fibonacci_qc))
println("Number of bonds: ", num_bonds(fibonacci_qc))

positions = get_positions(fibonacci_qc)
bonds = get_bonds(fibonacci_qc)
neighbors = get_nearest_neighbors(fibonacci_qc)

# Generate a 2D Penrose tiling
penrose_qc = generate_penrose_projection(3.0)
build_nearest_neighbor_bonds!(penrose_qc; cutoff=1.5)
```

### Periodic Lattice Unit Cells

```julia
# Get unit cell for square lattice
uc = get_unit_cell(Square)

# Access unit cell properties
println("Basis vectors: ", uc.basis)
println("Sublattice positions: ", uc.sublattice_positions)
println("Connections: ", uc.connections)
```

See `example/lattice_interface.jl` for a comprehensive demonstration of the unified interface.

## Additional Examples

Tight-binding model is available in `example/tight_binding.jl`
