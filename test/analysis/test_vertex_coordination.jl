@testset "vertex coordination & vertex types" begin
    @testset "Penrose substitution: coordination + vertex_type_counts" begin
        for gen in 2:3
            qc = generate_penrose_substitution(gen)
            n = num_sites(qc)
            @test n > 0

            # Build bonds with cutoff slightly above the rhombus edge
            # length (=1 for the substitution generator).
            build_nearest_neighbor_bonds!(qc; cutoff=1.1)

            # ---- coordination ----
            coords_all = coordination(qc)
            @test length(coords_all) == n
            @test all(c -> c isa Int, coords_all)
            @test minimum(coords_all) ≥ 1
            # The 8 standard P3 vertex configurations have valences in
            # {3, 4, 5, 6, 7}. Boundary vertices on a finite patch may
            # have lower coordination, so we only assert the upper
            # bound here.
            @test maximum(coords_all) ≤ 7

            # Per-site call agrees with the bulk vector.
            for i in 1:n
                @test coordination(qc, i) == coords_all[i]
            end

            # ---- vertex_type & vertex_type_counts ----
            counts = vertex_type_counts(qc)
            @test sum(values(counts)) == n
            allowed = Set([:Sun, :Star, :Ace, :Deuce, :Jack, :Queen, :King, :Other])
            for k in keys(counts)
                @test k in allowed
            end

            # Every per-site vertex_type call returns a symbol from
            # the allowed set.
            for i in 1:n
                @test vertex_type(qc, i) in allowed
            end
        end
    end

    @testset "coordination on Fibonacci 1D" begin
        qc = generate_fibonacci_substitution(5)
        build_nearest_neighbor_bonds!(qc; cutoff=1.1)
        n = num_sites(qc)
        coords_all = coordination(qc)
        @test length(coords_all) == n
        # 1D chain: interior sites have coordination 2, endpoints 1
        # (with cutoff just above the long Fibonacci interval, no
        # second-nearest pickups).
        @test all(c -> 0 ≤ c ≤ 2, coords_all)
    end

    @testset "coordination on Ammann–Beenker substitution" begin
        qc = generate_ammann_beenker_substitution(2)
        build_nearest_neighbor_bonds!(qc; cutoff=1.1)
        coords_all = coordination(qc)
        @test length(coords_all) == num_sites(qc)
        @test minimum(coords_all) ≥ 1
        # AB tilings: the 8 standard vertex configurations have
        # coordination ≤ 8 (the all-square 8-fold "Star" centre is
        # the maximum).
        @test maximum(coords_all) ≤ 8
    end

    @testset "bounds checking" begin
        qc = generate_penrose_substitution(2)
        n = num_sites(qc)
        @test_throws BoundsError coordination(qc, 0)
        @test_throws BoundsError coordination(qc, n + 1)
        @test_throws BoundsError vertex_type(qc, 0)
        @test_throws BoundsError vertex_type(qc, n + 1)
    end
end
