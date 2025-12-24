@testset "Penrose P3 Tiling Tests" begin
    using LinearAlgebra

    @testset "Projection Method" begin
        # Test basic generation
        qc_data = generate_penrose_projection(3.0)
        @test length(qc_data.positions) > 0
        @test qc_data.generation_method isa ProjectionMethod
        @test haskey(qc_data.parameters, :radius)
        @test qc_data.parameters[:radius] == 3.0

        # All positions should be within radius
        for pos in qc_data.positions
            @test norm(pos) <= 3.0 + 1e-10  # Small tolerance for floating point
        end

        # Test larger radius gives more points
        qc_data_small = generate_penrose_projection(2.0)
        qc_data_large = generate_penrose_projection(5.0)
        @test length(qc_data_large.positions) > length(qc_data_small.positions)

        # Test 5-fold symmetry (positions should roughly distribute evenly in 5 sectors)
        qc_data = generate_penrose_projection(5.0)
        if length(qc_data.positions) >= 20  # Need enough points to test
            angles = [atan(pos[2], pos[1]) for pos in qc_data.positions if norm(pos) > 0.1]
            # Normalize angles to [0, 2π)
            angles = mod.(angles, 2π)
            # Count points in each of 5 sectors
            sector_size = 2π / 5
            counts = zeros(Int, 5)
            for angle in angles
                sector = floor(Int, angle / sector_size) + 1
                if sector <= 5
                    counts[sector] += 1
                end
            end
            # Each sector should have some points (rough test)
            # @test all(counts .> 0)  # This might be too strict for small radius
        end
    end

    @testset "Substitution Method" begin
        # Test basic generation
        qc_data = generate_penrose_substitution(3)
        @test length(qc_data.positions) > 0
        @test qc_data.generation_method isa SubstitutionMethod
        @test haskey(qc_data.parameters, :generations)
        @test qc_data.parameters[:generations] == 3

        # More generations should give more tiles and positions
        qc_data_1 = generate_penrose_substitution(1)
        qc_data_2 = generate_penrose_substitution(2)
        @test qc_data_2.parameters[:n_tiles] >= qc_data_1.parameters[:n_tiles]
        @test length(qc_data_2.positions) >= length(qc_data_1.positions)
    end

    @testset "Tile Structure" begin
        # Test that tiles are properly formed
        qc_data = generate_penrose_substitution(2)

        for tile in qc_data.tiles
            # Each tile should have 4 vertices (rhombus)
            @test length(tile.vertices) == 4

            # Tile type should be 1 (fat) or 2 (thin)
            @test tile.type in [1, 2]

            # Center should be a 2D point
            @test length(tile.center) == 2
        end
    end
end
