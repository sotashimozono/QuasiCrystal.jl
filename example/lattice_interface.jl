"""
Example demonstrating the unified AbstractLattice interface.

This example shows how both periodic lattices and aperiodic quasicrystals
can be accessed through a common interface, enabling applications like
Lattice2DMonteCarlo to work uniformly with either structure type.
"""

using QuasiCrystal
using LinearAlgebra

println("=" ^ 70)
println("Unified AbstractLattice Interface Example")
println("=" ^ 70)

# ============================================================================
# 1. Creating a Quasicrystal (Fibonacci Lattice)
# ============================================================================

println("\n1. Creating a 1D Fibonacci Quasicrystal")
println("-" ^ 70)

# Generate a Fibonacci lattice using projection method
fibonacci_qc = generate_fibonacci_projection(20)

println("Generated Fibonacci lattice:")
println("  Number of sites: ", num_sites(fibonacci_qc))
println("  Dimension: 1D")
println("  Generation method: ", typeof(fibonacci_qc.generation_method))

# ============================================================================
# 2. Building Bonds and Nearest Neighbors
# ============================================================================

println("\n2. Building Bonds and Nearest Neighbors")
println("-" ^ 70)

# Build nearest neighbor bonds with a distance cutoff
build_nearest_neighbor_bonds!(fibonacci_qc; cutoff=2.0)

println("After building bonds:")
println("  Number of bonds: ", num_bonds(fibonacci_qc))
println(
    "  Average neighbors per site: ",
    sum(length.(get_nearest_neighbors(fibonacci_qc))) / num_sites(fibonacci_qc),
)

# ============================================================================
# 3. Accessing Lattice Properties via Unified Interface
# ============================================================================

println("\n3. Accessing Properties via Unified Interface")
println("-" ^ 70)

# These interface methods work for both quasicrystals and periodic lattices
positions = get_positions(fibonacci_qc)
bonds = get_bonds(fibonacci_qc)
neighbors = get_nearest_neighbors(fibonacci_qc)

println("Using unified interface methods:")
println("  get_positions() returns: ", typeof(positions))
println("  get_bonds() returns: ", typeof(bonds))
println("  get_nearest_neighbors() returns: ", typeof(neighbors))

# ============================================================================
# 4. Examining Bond Structure
# ============================================================================

println("\n4. Examining Bond Structure")
println("-" ^ 70)

if num_bonds(fibonacci_qc) > 0
    println("First 5 bonds:")
    for (i, bond) in enumerate(bonds[1:min(5, num_bonds(fibonacci_qc))])
        pos_i = positions[bond.src][1]
        pos_j = positions[bond.dst][1]
        println("  Bond $i: site $(bond.src) -> $(bond.dst)")
        println("    Positions: $pos_i -> $pos_j")
        println("    Distance: $(norm(bond.vector))")
    end
end

# ============================================================================
# 5. Creating a 2D Quasicrystal (Penrose Tiling)
# ============================================================================

println("\n5. Creating a 2D Penrose Quasicrystal")
println("-" ^ 70)

penrose_qc = generate_penrose_projection(3.0)

println("Generated Penrose tiling:")
println("  Number of sites: ", num_sites(penrose_qc))
println("  Dimension: 2D")

# Build bonds for 2D structure
build_nearest_neighbor_bonds!(penrose_qc; cutoff=1.5)

println("After building bonds:")
println("  Number of bonds: ", num_bonds(penrose_qc))
println("  Coordination numbers (neighbors per site):")

# Analyze coordination numbers
coord_numbers = [length(nn) for nn in get_nearest_neighbors(penrose_qc)]
println("    Min: ", minimum(coord_numbers))
println("    Max: ", maximum(coord_numbers))
println("    Mean: ", round(sum(coord_numbers) / length(coord_numbers); digits=2))

# ============================================================================
# 7. Type Hierarchy
# ============================================================================

println("\n7. Type Hierarchy")
println("-" ^ 70)

println("AbstractLattice hierarchy:")
println("  FibonacciLattice <: AbstractQuasicrystal{1} <: AbstractLattice{1}")
println("  PenroseP3 <: AbstractQuasicrystal{2} <: AbstractLattice{2}")
println("  AmmannBeenker <: AbstractQuasicrystal{2} <: AbstractLattice{2}")
println("  Square <: AbstractTopology{2} <: AbstractLattice{2}")
println("  UnitCell{2,Float64} <: AbstractLattice{2}")

println("\nQuasicrystalData instances work with the unified interface:")
println("  typeof(fibonacci_qc) = ", typeof(fibonacci_qc))
println("  typeof(penrose_qc) = ", typeof(penrose_qc))

# ============================================================================
# 8. Application Example: Computing Adjacency Information
# ============================================================================

println("\n8. Application Example: Computing Adjacency Information")
println("-" ^ 70)

function compute_adjacency_stats(lattice_data)
    """
    Example function that works with any AbstractLattice-compatible structure.
    This demonstrates how applications like Lattice2DMonteCarlo can uniformly
    access bond and site information.
    """
    n_sites = num_sites(lattice_data)
    n_bonds = num_bonds(lattice_data)
    neighbors = get_nearest_neighbors(lattice_data)

    # Compute statistics
    coord_numbers = [length(nn) for nn in neighbors]
    avg_coord = sum(coord_numbers) / n_sites

    return Dict(
        :n_sites => n_sites,
        :n_bonds => n_bonds,
        :avg_coordination => avg_coord,
        :min_coordination => minimum(coord_numbers),
        :max_coordination => maximum(coord_numbers),
    )
end

println("Computing adjacency statistics for Fibonacci lattice:")
fib_stats = compute_adjacency_stats(fibonacci_qc)
for (key, value) in fib_stats
    println("  $key: $value")
end

println("\nComputing adjacency statistics for Penrose tiling:")
penrose_stats = compute_adjacency_stats(penrose_qc)
for (key, value) in penrose_stats
    println("  $key: $value")
end

println("\n" * "=" ^ 70)
println("Example completed successfully!")
println("=" ^ 70)
