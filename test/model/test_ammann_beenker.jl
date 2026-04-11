@testset "Ammann–Beenker tiling" begin
    @testset "projection method" begin
        qc = generate_ammann_beenker_projection(3.0)
        @test qc isa AbstractLattice{2,Float64}
        @test num_sites(qc) > 0
        @test qc.generation_method isa ProjectionMethod
        @test haskey(qc.parameters, :radius)
        @test qc.parameters[:radius] == 3.0
        @test qc.parameters[:symmetry] == 8

        # Every site is within the requested radius.
        for i in 1:num_sites(qc)
            @test norm(position(qc, i)) <= 3.0 + 1e-10
        end

        # Larger radius → more sites.
        qc_small = generate_ammann_beenker_projection(2.0)
        qc_large = generate_ammann_beenker_projection(5.0)
        @test num_sites(qc_large) > num_sites(qc_small)

        # Positions are SVector{2, Float64}.
        @test position(qc, 1) isa SVector{2,Float64}
    end

    @testset "substitution method (placeholder inflation)" begin
        qc = generate_ammann_beenker_substitution(2)
        @test num_sites(qc) > 0
        @test qc.generation_method isa SubstitutionMethod
        @test qc.parameters[:generations] == 2
        @test qc.parameters[:symmetry] == 8

        # Deeper generations produce at least as many tiles.
        qc_1 = generate_ammann_beenker_substitution(1)
        qc_2 = generate_ammann_beenker_substitution(2)
        @test qc_2.parameters[:n_tiles] >= qc_1.parameters[:n_tiles]
    end

    @testset "LatticeCore traits" begin
        qc = generate_ammann_beenker_projection(4.0)
        @test periodicity(qc) isa Aperiodic
        @test reciprocal_support(qc) isa HasFourierModule
        @test is_finite(qc) == true
    end

    @testset "distance-based bond builder" begin
        qc = generate_ammann_beenker_projection(4.0)
        @test isempty(bonds(qc))

        build_nearest_neighbor_bonds!(qc; cutoff=1.1)
        @test !isempty(bonds(qc))
        for b in bonds(qc)
            @test b isa Bond{2,Float64}
            @test norm(b.vector) < 1.1
        end
    end
end
