using Test
using QuasiCrystal
using Plots

@testset "plot_tiles (Plots ext)" begin
    @testset "Penrose tiles" begin
        qc = generate_penrose_substitution(2)
        @test !isempty(qc.tiles)

        # Default palette path.
        p = plot_tiles(qc)
        @test p isa Plots.Plot

        # Each preset must build without error.
        for preset in (:default, :pastel, :bw)
            @test plot_tiles(qc; palette=preset) isa Plots.Plot
        end

        # Custom Dict palette (partial dict -- unspecified types fall
        # back to the :default palette).
        custom = Dict{TileType,Any}(FatRhombus() => :crimson)
        @test plot_tiles(qc; palette=custom) isa Plots.Plot

        # Boundary toggle and legend toggle.
        @test plot_tiles(qc; show_boundary=false) isa Plots.Plot
        @test plot_tiles(qc; legend=false) isa Plots.Plot

        # Custom title forwarded.
        @test plot_tiles(qc; title="My Penrose") isa Plots.Plot
    end

    @testset "Ammann-Beenker tiles" begin
        ab = generate_ammann_beenker_substitution(2)
        @test !isempty(ab.tiles)

        p = plot_tiles(ab)
        @test p isa Plots.Plot

        # Boundary toggle and explicit Square override.
        custom = Dict{TileType,Any}(Square() => :gold)
        @test plot_tiles(ab; palette=custom, show_boundary=false) isa Plots.Plot
    end

    @testset "Error paths" begin
        # 1D quasicrystals do not carry tiles -- recipe must reject them.
        fib = generate_fibonacci_substitution(4)
        @test_throws ArgumentError plot_tiles(fib)

        # Projection-method 2D generators do not populate tiles.
        empty_qc = generate_penrose_projection(1.5)
        @test isempty(empty_qc.tiles)
        @test_throws ArgumentError plot_tiles(empty_qc)

        # Unknown palette symbol.
        qc = generate_penrose_substitution(1)
        @test_throws ArgumentError plot_tiles(qc; palette=:nope)

        # Wrong palette type (not Symbol or Dict).
        @test_throws ArgumentError plot_tiles(qc; palette=42)
    end
end
