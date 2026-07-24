using QuasiCrystal
using LatticeCore: bond_center
using StaticArrays
using Test

const _CB_PHI = (1 + sqrt(5)) / 2

# undirected (i,j) edge set of a lattice's bonds
_edges(d) = Set((min(b.i, b.j), max(b.i, b.j)) for b in d.bonds)

@testset "build_chain_bonds! wires the chain in position order" begin
    d = generate_fibonacci_substitution(4)
    build_chain_bonds!(d)
    @test length(d.bonds) == num_sites(d) - 1
    # Sorted by position, consecutive sites are bonded.
    order = sortperm(d.positions; by=p -> p[1])
    @test _edges(d) == Set(
        (min(order[k], order[k + 1]), max(order[k], order[k + 1])) for
        k in 1:(num_sites(d) - 1)
    )
    # Every interior site has coordination 2, ends have 1.
    coords = sort(length.(d.nearest_neighbors))
    @test coords[1] == 1 && coords[end] == 2
end

@testset "chain bonds are inflation-covariant (scale-invariant adjacency)" begin
    inf = InfiniteQuasicrystal(FibonacciLattice())
    base = build_chain_bonds!(materialize(inf; n_points=60))
    up = build_chain_bonds!(materialize(inflate(inf); n_points=60))

    # Same combinatorial bond graph...
    @test _edges(base) == _edges(up)
    @test length(base.bonds) == length(up.bonds)
    # ...and every bond vector scaled by exactly λ = ϕ.
    bvec(d) = sort([b.vector[1] for b in d.bonds])
    @test bvec(up) ≈ _CB_PHI .* bvec(base)
    # Bond centres scale by λ too (bond-centring is covariant).
    @test sort([bond_center(up, b)[1] for b in up.bonds]) ≈ _CB_PHI .* sort([bond_center(base, b)[1] for b in base.bonds])
end

@testset "a fixed cutoff is NOT inflation-covariant (contrast)" begin
    # The distance-cutoff builder, with the cutoff held fixed, produces a
    # different bond graph once the chain is inflated — the very failure
    # build_chain_bonds! avoids.
    inf = InfiniteQuasicrystal(FibonacciLattice())
    base = build_nearest_neighbor_bonds!(materialize(inf; n_points=60); cutoff=1.2)
    up = build_nearest_neighbor_bonds!(materialize(inflate(inf); n_points=60); cutoff=1.2)
    @test length(base.bonds) != length(up.bonds)
    # Whereas chain bonds keep the same count.
    cbase = build_chain_bonds!(materialize(inf; n_points=60))
    cup = build_chain_bonds!(materialize(inflate(inf); n_points=60))
    @test length(cbase.bonds) == length(cup.bonds)
end
