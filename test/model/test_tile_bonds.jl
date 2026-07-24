using QuasiCrystal
using LatticeCore: plaquettes
using LinearAlgebra
using Test

const _SUBST_2D = (
    ("penrose", generate_penrose_substitution), ("ab", generate_ammann_beenker_substitution)
)

_edges(d) = Set((min(b.i, b.j), max(b.i, b.j)) for b in d.bonds)
_lengths(d) = [norm(b.vector) for b in d.bonds]

@testset "build_tile_bonds! builds the unit-edge tiling graph" begin
    for (name, gen) in _SUBST_2D
        d = gen(3)
        build_tile_bonds!(d)

        @test !isempty(d.bonds)
        # Every bond is a tile edge — unit length for the rhombus tilings.
        @test all(≈(1.0; atol=1e-8), _lengths(d))
        # Undirected, deduplicated: one Bond per edge, and the adjacency
        # lists agree with the bond set.
        @test length(d.bonds) == length(_edges(d))
        deg = sum(length, d.nearest_neighbors)
        @test deg == 2 * length(d.bonds)          # each edge contributes to two sites
        # No self-loops, indices in range.
        @test all(b -> b.i != b.j, d.bonds)
        @test all(b -> 1 <= b.i <= num_sites(d) && 1 <= b.j <= num_sites(d), d.bonds)
    end
end

@testset "build_tile_bonds! is parameter-free across inflation generations" begin
    # The tiling-intrinsic graph stays the unit-edge network at every
    # generation — no cutoff to retune — unlike a fixed-distance builder.
    for (name, gen) in _SUBST_2D
        for g in 2:4
            d = build_tile_bonds!(gen(g))
            @test all(≈(1.0; atol=1e-8), _lengths(d))
            @test !isempty(d.bonds)
        end
    end
end

@testset "build_tile_bonds! independently reproduces the tile-edge set" begin
    d = generate_penrose_substitution(3)
    build_tile_bonds!(d)
    # Recompute the undirected tile-boundary edges directly.
    expected = Set{Tuple{Int,Int}}()
    for p in plaquettes(d)
        vs = p.vertices
        n = length(vs)
        for k in 1:n
            i, j = vs[k], vs[mod1(k + 1, n)]
            i == j && continue
            push!(expected, (min(i, j), max(i, j)))
        end
    end
    @test _edges(d) == expected
end

@testset "build_tile_bonds! refuses a tiling-less patch" begin
    # Projection-generated patches carry no tiles.
    proj = generate_penrose_projection(6.0)
    @test isempty(proj.tiles)
    @test_throws ArgumentError build_tile_bonds!(proj)

    # A Fibonacci chain (1D) also has no tiles — use build_chain_bonds!.
    chain = generate_fibonacci_substitution(4)
    @test_throws ArgumentError build_tile_bonds!(chain)
end
