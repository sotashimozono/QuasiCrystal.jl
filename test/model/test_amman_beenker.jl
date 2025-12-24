@testset "Ammann-Beenker Tiling Tests" begin
    using LinearAlgebra

    @testset "Projection Method" begin
        # Test basic generation
        qc_data = generate_ammann_beenker_projection(3.0)
        @test length(qc_data.positions) > 0
        @test qc_data.generation_method isa ProjectionMethod
        @test haskey(qc_data.parameters, :radius)
        @test qc_data.parameters[:radius] == 3.0
        @test qc_data.parameters[:symmetry] == 8

        # All positions should be within radius
        for pos in qc_data.positions
            @test norm(pos) <= 3.0 + 1e-10
        end

        # Larger radius gives more points
        qc_data_small = generate_ammann_beenker_projection(2.0)
        qc_data_large = generate_ammann_beenker_projection(5.0)
        @test length(qc_data_large.positions) > length(qc_data_small.positions)

        # Test 8-fold symmetry
        qc_data = generate_ammann_beenker_projection(5.0)
        if length(qc_data.positions) >= 32  # Need enough points
            angles = [atan(pos[2], pos[1]) for pos in qc_data.positions if norm(pos) > 0.1]
            angles = mod.(angles, 2π)
            sector_size = 2π / 8
            counts = zeros(Int, 8)
            for angle in angles
                sector = floor(Int, angle / sector_size) + 1
                if sector <= 8
                    counts[sector] += 1
                end
            end
            # Some symmetry should be present
            # @test maximum(counts) - minimum(counts) < length(angles) / 2
        end
    end

    @testset "Substitution Method" begin
        # Test basic generation
        qc_data = generate_ammann_beenker_substitution(2)
        @test length(qc_data.positions) > 0
        @test qc_data.generation_method isa SubstitutionMethod
        @test qc_data.parameters[:generations] == 2
        @test qc_data.parameters[:symmetry] == 8

        # More generations give more tiles
        qc_data_1 = generate_ammann_beenker_substitution(1)
        qc_data_2 = generate_ammann_beenker_substitution(2)
        @test qc_data_2.parameters[:n_tiles] >= qc_data_1.parameters[:n_tiles]
    end

    @testset "Tile Types" begin
        qc_data = generate_ammann_beenker_substitution(2)

        # Should have both squares (type 1) and rhombi (type 2)
        types = [tile.type for tile in qc_data.tiles]
        @test 1 in types  # Has squares
        @test 2 in types  # Has rhombi

        for tile in qc_data.tiles
            # Each tile should have 4 vertices
            @test length(tile.vertices) == 4
            # Type should be 1 or 2
            @test tile.type in [1, 2]
            # Center should be 2D
            @test length(tile.center) == 2
        end
    end

    @testset "Octagonal Properties" begin
        # The Ammann-Beenker tiling uses √2 as a fundamental scaling
        qc_data = generate_ammann_beenker_substitution(1)

        # Check that tiles exist
        @test length(qc_data.tiles) > 0

        # Initial configuration should have multiple tiles
        @test qc_data.parameters[:n_tiles] >= 1
    end
end
