using QuasiCrystal
using QuasiCrystal: Square, RhombusAB, tile_parentage
using LatticeCore: plaquettes
using LinearAlgebra
using Test

# diagonals of a reconstructed tile (vertices in cyclic order v,a,w,b)
_diags(t) = (norm(t.vertices[1] - t.vertices[3]), norm(t.vertices[2] - t.vertices[4]))
function _edges_ok(t)
    return all(
        abs(norm(t.vertices[i] - t.vertices[mod1(i + 1, 4)]) - 1.0) < 1e-6 for i in 1:4
    )
end

@testset "AB projection reconstructs a valid square+rhombus tiling" begin
    d = generate_ammann_beenker_projection(6.0)
    @test !isempty(d.tiles)

    sq = filter(t -> t.type isa Square, d.tiles)
    rh = filter(t -> t.type isa RhombusAB, d.tiles)
    @test !isempty(sq) && !isempty(rh)
    @test length(sq) + length(rh) == length(d.tiles)

    # every tile has unit edges
    @test all(_edges_ok, d.tiles)

    # square diagonals are both √2; rhombus diagonals are 2sin(π/8) and 2cos(π/8)
    for t in sq
        e, f = _diags(t)
        @test e ≈ sqrt(2) atol = 1e-6
        @test f ≈ sqrt(2) atol = 1e-6
    end
    short, long = 2sin(π / 8), 2cos(π / 8)
    for t in rh
        e, f = sort(collect(_diags(t)))
        @test e ≈ short atol = 1e-6
        @test f ≈ long atol = 1e-6
    end
end

@testset "build_tile_bonds! now works on the AB projection" begin
    d = generate_ammann_beenker_projection(6.0)
    build_tile_bonds!(d)
    @test !isempty(d.bonds)
    @test all(b -> abs(norm(b.vector) - 1.0) < 1e-6, d.bonds)
    # undirected dedup
    edges = Set((min(b.i, b.j), max(b.i, b.j)) for b in d.bonds)
    @test length(edges) == length(d.bonds)
    # plaquette centring is now available for AB (materialised from tiles)
    @test length(plaquettes(d)) == length(d.tiles)
end

@testset "AB tile_parentage: exact partition by coarse-tile containment" begin
    # AB is not a stone inflation (fine tiles straddle coarse tile
    # boundaries), but every fine-tile centre lies in exactly one coarse
    # tile, so centre containment against the k-inflated tiling gives an
    # exact tile parentage — the AB analogue of Penrose's tile_parentage.
    for R in (6.0, 8.0), k in 1:2
        d = generate_ammann_beenker_projection(R)
        tp = tile_parentage(d, k)
        @test sort(vcat(tp...)) == collect(1:length(d.tiles))   # exact partition
        @test all(!isempty, tp)
    end

    d = generate_ammann_beenker_projection(8.0)
    # coarsening: fewer cells with more deflation
    @test length(tile_parentage(d, 2)) <= length(tile_parentage(d, 1))
    @test_throws ArgumentError tile_parentage(d, 0)                       # k >= 1
    # the AB substitution generator is a placeholder (no :radius) → refused
    @test_throws ArgumentError tile_parentage(generate_ammann_beenker_substitution(2), 1)
end
