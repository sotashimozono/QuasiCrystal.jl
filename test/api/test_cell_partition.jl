using QuasiCrystal
using LatticeCore: cell_partition, rescale
using Test

@testset "cell_partition: hand-verified Fibonacci example (g=3, k=1)" begin
    # word(3) = L S L L S (5 intervals, 6 vertices). Deflating one
    # generation gives word(2) = L S L (3 tiles), and the substitution
    # L→LS, S→L labels the fine intervals [1,1,2,3,3]:
    #   cell 1 = {1,2}, cell 2 = {3}, cell 3 = {4,5} + final vertex 6.
    d = generate_fibonacci_substitution(3)
    @test num_sites(d) == 6
    @test cell_partition(d, 1) == [[1, 2], [3], [4, 5, 6]]
end

@testset "cell_partition is an exact partition matching the coarse tiling" begin
    for g in 2:8, k in 1:g
        d = generate_fibonacci_substitution(g)
        cells = cell_partition(d, k)

        # Exact partition of all vertices.
        flat = sort(vcat(cells...))
        @test flat == collect(1:num_sites(d))
        @test sum(length, cells) == num_sites(d)

        # One cell per coarse tile = (coarse vertices − 1).
        coarse = rescale(d, -k)
        @test length(cells) == num_sites(coarse) - 1

        # Each cell's vertices are contiguous and the cells are ordered.
        for c in cells
            @test c == first(c):last(c)
        end
        @test issorted(first.(cells))
    end
end

@testset "cell_partition composes with the size sequence" begin
    # Deflating k in one step must group the same way as the coarse
    # lattice's own tile count at each intermediate generation.
    d = generate_fibonacci_substitution(7)
    for k in 1:5
        @test length(cell_partition(d, k)) ==
            num_sites(generate_fibonacci_substitution(7 - k)) - 1
    end
end

@testset "cell_partition rejects the ill-posed cases" begin
    # Projection patches: no substitution parentage.
    proj = generate_fibonacci_projection(50)
    @test_throws ArgumentError cell_partition(proj, 1)

    # Other substitution families (re-centring inflation).
    pen = generate_penrose_substitution(2)
    @test_throws ArgumentError cell_partition(pen, 1)

    d = generate_fibonacci_substitution(3)
    @test_throws ArgumentError cell_partition(d, 0)         # k >= 1
    @test_throws ArgumentError cell_partition(d, 4)         # k > generation
end
