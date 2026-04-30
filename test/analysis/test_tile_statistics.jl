@testset "tile statistics" begin
    @testset "tile_area constants" begin
        # Penrose
        @test tile_area(FatRhombus()) ≈ sin(2π / 5)
        @test tile_area(ThinRhombus()) ≈ sin(π / 5)
        # Ammann–Beenker
        @test tile_area(Square()) ≈ 1.0
        @test tile_area(RhombusAB()) ≈ sin(π / 4)
        # Sanity: fat > thin
        @test tile_area(FatRhombus()) > tile_area(ThinRhombus())
    end

    @testset "tile_counts / tile_density (Penrose)" begin
        qc = generate_penrose_substitution(4)
        counts = tile_counts(qc)
        # Total tile count is consistent with data.tiles and num_plaquettes.
        @test sum(values(counts)) == length(qc.tiles)
        @test sum(values(counts)) == num_plaquettes(qc)

        # Both Penrose tile types are present in a non-trivial tiling.
        @test haskey(counts, FatRhombus())
        @test haskey(counts, ThinRhombus())
        @test counts[FatRhombus()] > 0
        @test counts[ThinRhombus()] > 0

        densities = tile_density(qc)
        @test sum(values(densities)) ≈ 1.0
        for (t, c) in counts
            @test densities[t] ≈ c / length(qc.tiles)
        end
    end

    @testset "tile_counts on Fibonacci (no tiles)" begin
        # 1D Fibonacci does not populate `tiles`; the API still works.
        qc = generate_fibonacci_substitution(4)
        @test isempty(qc.tiles)
        @test isempty(tile_counts(qc))
        @test isempty(tile_density(qc))
    end

    @testset "golden_ratio_check API surface" begin
        # `golden_ratio_check` is exercised on (i) synthetic tilings
        # with engineered tile counts (covering both ok=true and
        # ok=false branches), and (ii) the actual substitution
        # generator (issue #60: now uses Robinson-triangle inflation,
        # so the asymptotic `#fat / #thin → ϕ` is recovered up to
        # boundary effects).
        function _synthetic(n_fat::Int, n_thin::Int)
            tiles = Tile{2,Float64}[]
            v0 = SVector{2,Float64}(0.0, 0.0)
            for k in 1:n_fat
                push!(
                    tiles,
                    Tile{2,Float64}(
                        [v0, v0, v0, v0], FatRhombus(), SVector{2,Float64}(k, 0.0)
                    ),
                )
            end
            for k in 1:n_thin
                push!(
                    tiles,
                    Tile{2,Float64}(
                        [v0, v0, v0, v0], ThinRhombus(), SVector{2,Float64}(k, 1.0)
                    ),
                )
            end
            positions = SVector{2,Float64}[v0]
            return QuasicrystalData{2,Float64}(
                PenroseP3(), positions, tiles, SubstitutionMethod(), Dict{Symbol,Any}()
            )
        end

        # 89 / 55 = 1.6181818... — Fibonacci ratio, ≈ ϕ to 4 decimals.
        good = _synthetic(89, 55)
        r_good = golden_ratio_check(good)
        @test r_good.expected ≈ ϕ
        @test r_good.observed ≈ 89 / 55
        @test r_good.ok == true

        # 1 / 1 ratio is far from ϕ; with default tol=0.05 → ok=false.
        bad = _synthetic(10, 10)
        r_bad = golden_ratio_check(bad)
        @test r_bad.observed ≈ 1.0
        @test r_bad.ok == false

        # Larger tol can flip it back to true.
        r_bad_loose = golden_ratio_check(bad; tol=1.0)
        @test r_bad_loose.ok == true

        # No-thin-rhombi → ArgumentError.
        only_fat = _synthetic(5, 0)
        @test_throws ArgumentError golden_ratio_check(only_fat)

        # Real generator (Robinson-triangle inflation, issue #60).
        # Generation 4 already lands inside the default tol=0.05
        # window of `ϕ`; the bulk recurrence is exact, so only
        # patch-boundary half-tiles drag the ratio.
        qc = generate_penrose_substitution(4)
        r_gen = golden_ratio_check(qc)
        @test propertynames(r_gen) == (:observed, :expected, :ok)
        @test r_gen.expected ≈ ϕ
        @test r_gen.observed > 0
        @test r_gen.ok == true
    end

    @testset "Ammann–Beenker tile mix" begin
        qc = generate_ammann_beenker_substitution(3)
        counts = tile_counts(qc)
        @test sum(values(counts)) == length(qc.tiles)
        @test haskey(counts, Square())
        @test haskey(counts, RhombusAB())

        # The AB tile mix is *not* the golden ratio — Ammann–Beenker has
        # its own characteristic ratio (asymptotically related to √2,
        # not ϕ). We therefore expect the golden-ratio relative deviation
        # to be clearly non-zero.
        ratio_rs = counts[RhombusAB()] / counts[Square()]
        @test abs(ratio_rs - ϕ) / ϕ > 0.05  # not golden
    end

    @testset "tile_perimeter" begin
        # Empty tiling
        qc_fib = generate_fibonacci_substitution(3)
        @test tile_perimeter(qc_fib) == 0.0

        # Penrose: every tile has 4 unit edges, so the perimeter is
        # bounded above by 4 * num_tiles (strict inequality whenever
        # any interior edge is shared).
        qc = generate_penrose_substitution(3)
        peri = tile_perimeter(qc)
        @test peri > 0
        @test peri < 4 * length(qc.tiles) + 1e-8

        # Lower bound: in a connected tiling the perimeter is at least
        # the number of distinct edges, which is at least num_tiles
        # (each tile contributes at least one edge to the union).
        @test peri ≥ length(qc.tiles) - 1e-8

        # Refining the tiling grows the perimeter monotonically.
        qc2 = generate_penrose_substitution(4)
        @test tile_perimeter(qc2) > peri
    end
end
