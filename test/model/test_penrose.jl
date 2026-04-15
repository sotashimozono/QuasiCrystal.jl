@testset "Penrose P3 tiling" begin
    @testset "projection method" begin
        qc = generate_penrose_projection(3.0)
        @test qc isa AbstractLattice{2,Float64}
        @test num_sites(qc) > 0
        @test qc.generation_method isa ProjectionMethod
        @test haskey(qc.parameters, :radius)
        @test qc.parameters[:radius] == 3.0

        # Every site is within the requested radius.
        for i in 1:num_sites(qc)
            @test norm(position(qc, i)) <= 3.0 + 1e-10
        end

        # Larger radius gives more points.
        qc_small = generate_penrose_projection(2.0)
        qc_large = generate_penrose_projection(5.0)
        @test num_sites(qc_large) > num_sites(qc_small)

        # Positions are SVector{2, Float64}.
        @test position(qc, 1) isa SVector{2,Float64}
    end

    @testset "substitution method (placeholder inflation)" begin
        qc = generate_penrose_substitution(4)
        @test num_sites(qc) > 0
        @test qc.generation_method isa SubstitutionMethod
        @test qc.parameters[:generations] == 4

        # For strict Robinson deflation, boundary half-tiles are dropped.
        # So generations 1 or 2 actually lose full tiles. 
        # By generation 4 and 5, bulk area dominates and tile count increases.
        qc_4 = generate_penrose_substitution(4)
        qc_5 = generate_penrose_substitution(5)
        @test qc_5.parameters[:n_tiles] > qc_4.parameters[:n_tiles]
        @test num_sites(qc_5) > num_sites(qc_4)
    end

    @testset "LatticeCore traits" begin
        qc = generate_penrose_projection(4.0)
        @test periodicity(qc) isa Aperiodic
        @test reciprocal_support(qc) isa HasFourierModule
        @test is_finite(qc) == true
        @test site_layout(qc) isa UniformLayout
        @test site_type(qc, 1) isa IsingSite
    end

    @testset "distance-based bond builder" begin
        qc = generate_penrose_projection(3.0)
        @test isempty(bonds(qc))

        build_nearest_neighbor_bonds!(qc; cutoff=1.2)
        # At least some bonds should be found in a 3.0-radius sample.
        @test !isempty(bonds(qc))

        for b in bonds(qc)
            @test b isa Bond{2,Float64}
            @test norm(b.vector) < 1.2
            @test 1 <= b.i <= num_sites(qc)
            @test 1 <= b.j <= num_sites(qc)
        end
    end
end
