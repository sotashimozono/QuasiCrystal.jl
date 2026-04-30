using QuasiCrystal, Test, StaticArrays
using Plots  # triggers QuasiCrystalPlotsExt

@testset "plot_acceptance_window" begin
    @testset "Fibonacci :interval_1d" begin
        qc = generate_fibonacci_projection(40)
        @test window_shape(qc) === :interval_1d

        plt = plot_acceptance_window(qc; show_hyper_points=true, n_hyper=4)
        @test plt isa Plots.Plot
        # window line + inside scatter + outside scatter + endpoint vlines
        @test length(plt.series_list) >= 2

        plt2 = plot_acceptance_window(qc; show_hyper_points=false)
        @test plt2 isa Plots.Plot
    end

    @testset "Ammann–Beenker :box_2d" begin
        qc = generate_ammann_beenker_projection(2.0)
        @test window_shape(qc) === :box_2d

        plt = plot_acceptance_window(qc; show_hyper_points=true, n_hyper=3)
        @test plt isa Plots.Plot
        # rectangle (shape) + at least one scatter overlay
        @test length(plt.series_list) >= 2

        plt2 = plot_acceptance_window(qc; show_hyper_points=false)
        @test plt2 isa Plots.Plot
    end

    @testset "Penrose P3 :box_3d" begin
        qc = generate_penrose_projection(2.0)
        @test window_shape(qc) === :box_3d

        plt = plot_acceptance_window(qc; show_hyper_points=true, n_hyper=2)
        @test plt isa Plots.Plot
        # Three sub-panels (y₁y₂, y₁y₃, y₂y₃)
        @test length(plt.subplots) >= 3

        plt2 = plot_acceptance_window(qc; show_hyper_points=false, n_hyper=2)
        @test plt2 isa Plots.Plot
    end

    @testset "substitution :none raises" begin
        fib_s = generate_fibonacci_substitution(4)
        @test window_shape(fib_s) === :none
        @test_throws ArgumentError plot_acceptance_window(fib_s)

        ab_s = generate_ammann_beenker_substitution(2)
        @test window_shape(ab_s) === :none
        @test_throws ArgumentError plot_acceptance_window(ab_s)

        pen_s = generate_penrose_substitution(2)
        @test window_shape(pen_s) === :none
        @test_throws ArgumentError plot_acceptance_window(pen_s)
    end

    @testset "unknown window_shape raises" begin
        qc = generate_fibonacci_projection(20)
        delete!(qc.parameters, :window_shape)
        @test window_shape(qc) === :unknown
        @test_throws ArgumentError plot_acceptance_window(qc)
    end
end
