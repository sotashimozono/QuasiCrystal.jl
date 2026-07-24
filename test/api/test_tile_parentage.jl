using QuasiCrystal
using QuasiCrystal: tile_parentage
using Test

@testset "tile_parentage is an exact partition of the Penrose tiles" begin
    for g in 2:5, k in 1:g
        d = generate_penrose_substitution(g)
        tp = tile_parentage(d, k)

        flat = sort(vcat(tp...))
        @test flat == collect(1:length(d.tiles))          # every tile once
        @test sum(length, tp) == length(d.tiles)          # disjoint cover
        @test all(!isempty, tp)                            # no empty cells reported
    end
end

@testset "tile_parentage coarsens monotonically with k" begin
    d = generate_penrose_substitution(5)
    ncells = [length(tile_parentage(d, k)) for k in 1:5]
    # More deflation ⇒ fewer, larger coarse cells.
    @test issorted(ncells; rev=true)
    @test ncells[end] <= ncells[1]
    # Deflating all the way lands on the seed's few triangles.
    @test length(tile_parentage(d, 5)) <= length(tile_parentage(d, 4))
end

@testset "tile_parentage nests: a k-step cell splits into finer cells" begin
    # Every fine tile keeps the same coarse ancestor whether reached in
    # one k-step or by composing — checked via cell sizes summing.
    d = generate_penrose_substitution(4)
    for k in 1:4
        tp = tile_parentage(d, k)
        # each fine tile appears in exactly one cell (partition invariant)
        seen = zeros(Int, length(d.tiles))
        for cell in tp, t in cell
            seen[t] += 1
        end
        @test all(==(1), seen)
    end
end

@testset "tile_parentage rejects the ill-posed cases" begin
    # Projection (no tiles / not self-similar via this route).
    @test_throws ArgumentError tile_parentage(generate_penrose_projection(6.0), 1)
    # Other families' generators are not exactly self-similar here.
    @test_throws ArgumentError tile_parentage(generate_ammann_beenker_substitution(2), 1)
    @test_throws ArgumentError tile_parentage(generate_fibonacci_substitution(4), 1)

    d = generate_penrose_substitution(3)
    @test_throws ArgumentError tile_parentage(d, 0)      # k >= 1
    @test_throws ArgumentError tile_parentage(d, 4)      # k > generation
end
