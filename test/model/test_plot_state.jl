using QuasiCrystal, Test
using Plots

@testset "plot_state (Plots ext)" begin
    @testset "real-valued state on 2D Penrose" begin
        qc = generate_penrose_substitution(2)
        n = num_sites(qc)
        state = [Float64(i) / n for i in 1:n]
        p = plot_state(qc, state)
        @test p isa Plots.Plot

        # Forwarded kwargs must not error.
        p2 = plot_state(qc, state; colormap=:plasma, marker_size=5, title="custom")
        @test p2 isa Plots.Plot
    end

    @testset "real-valued state on 1D Fibonacci" begin
        fib = generate_fibonacci_substitution(5)
        n = num_sites(fib)
        state = collect(range(-1.0, 1.0; length=n))
        p = plot_state(fib, state)
        @test p isa Plots.Plot
    end

    @testset "discrete (Bool / Int) palette path" begin
        qc = generate_ammann_beenker_projection(3.0)
        n = num_sites(qc)
        bool_state = [iseven(i) for i in 1:n]
        p_bool = plot_state(qc, bool_state)
        @test p_bool isa Plots.Plot

        int_state = [mod(i, 3) for i in 1:n]
        p_int = plot_state(qc, int_state; marker_size=3)
        @test p_int isa Plots.Plot
    end

    @testset "complex state with mode projections" begin
        qc = generate_penrose_substitution(2)
        n = num_sites(qc)
        psi = [cis(2pi * i / n) / sqrt(n) for i in 1:n]
        for mode in (:abs2, :abs, :real, :imag, :phase)
            p = plot_state(qc, psi; mode=mode)
            @test p isa Plots.Plot
        end
        @test_throws ArgumentError plot_state(qc, psi; mode=:bogus)
    end

    @testset "input validation" begin
        qc = generate_penrose_substitution(2)
        bad = zeros(num_sites(qc) + 1)
        @test_throws DimensionMismatch plot_state(qc, bad)
    end
end
