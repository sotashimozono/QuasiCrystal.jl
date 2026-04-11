@testset "Fibonacci lattice" begin
    @testset "projection method" begin
        qc = generate_fibonacci_projection(10)

        # LatticeCore contract
        @test qc isa AbstractLattice{1,Float64}
        @test num_sites(qc) <= 10
        @test num_sites(qc) > 0
        @test qc.generation_method isa ProjectionMethod
        @test haskey(qc.parameters, :n_points)

        # Positions are SVector{1, Float64}, sorted, non-negative
        ps = [position(qc, i)[1] for i in 1:num_sites(qc)]
        @test issorted(ps)
        @test all(ps .>= 0)

        # Larger requested size gives strictly more points
        qc_large = generate_fibonacci_projection(50)
        @test num_sites(qc_large) > num_sites(qc)
        @test num_sites(qc_large) <= 50
    end

    @testset "substitution method" begin
        qc = generate_fibonacci_substitution(5)
        @test num_sites(qc) > 0
        @test qc.generation_method isa SubstitutionMethod
        @test haskey(qc.parameters, :generations)
        @test qc.parameters[:generations] == 5

        # Sequence length follows the Fibonacci numbers (off by one
        # because we prepend position 0).
        for n in 1:7
            qc_n = generate_fibonacci_substitution(n)
            expected_length = fibonacci_sequence_length(n)
            @test num_sites(qc_n) == expected_length + 1
        end

        # Ratio of long to short spacings approaches ϕ.
        qc_deep = generate_fibonacci_substitution(10)
        positions_1d = [position(qc_deep, i)[1] for i in 1:num_sites(qc_deep)]
        spacings = diff(positions_1d)
        long_spacings = filter(x -> x > 1.5, spacings)
        short_spacings = filter(x -> x <= 1.5, spacings)

        if !isempty(long_spacings) && !isempty(short_spacings)
            avg_long = sum(long_spacings) / length(long_spacings)
            avg_short = sum(short_spacings) / length(short_spacings)
            @test abs(avg_long / avg_short - ϕ) < 0.1
        end
    end

    @testset "fibonacci_sequence_length helper" begin
        @test fibonacci_sequence_length(0) == 1
        @test fibonacci_sequence_length(1) == 2
        @test fibonacci_sequence_length(2) == 3
        @test fibonacci_sequence_length(3) == 5
        @test fibonacci_sequence_length(4) == 8
        @test fibonacci_sequence_length(5) == 13
        @test fibonacci_sequence_length(6) == 21
    end

    @testset "method comparison" begin
        proj = generate_fibonacci_projection(20)
        subst = generate_fibonacci_substitution(6)   # 21 + 1 positions
        @test abs(num_sites(proj) - num_sites(subst)) < 5
    end

    @testset "distance-based bond builder" begin
        qc = generate_fibonacci_substitution(5)
        @test isempty(bonds(qc))

        build_nearest_neighbor_bonds!(qc; cutoff=ϕ + 0.5)
        @test !isempty(bonds(qc))

        # Every bond should have the expected displacement magnitude.
        for b in bonds(qc)
            @test b isa Bond{1,Float64}
            @test norm(b.vector) < ϕ + 0.5
            @test norm(b.vector) > 0
        end

        # Neighbour lists are symmetric.
        for i in 1:num_sites(qc)
            for j in neighbors(qc, i)
                @test i in neighbors(qc, j)
            end
        end
    end
end
