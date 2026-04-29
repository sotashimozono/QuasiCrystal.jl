using QuasiCrystal, Test, StaticArrays

@testset "numerics tolerance constants" begin
    @test QuasiCrystal.VERTEX_MERGE_TOL == 1e-10
    @test QuasiCrystal.POSITION_TOLERANCE === QuasiCrystal.VERTEX_MERGE_TOL
    @test QuasiCrystal.SNAP_GRID_EPS == 1e-5
    @test QuasiCrystal.STAR_DIRECTION_TOL == 1e-4

    a = SVector(0.0, 0.0)
    b = SVector(1e-12, 0.0)
    c = SVector(0.5, 0.0)
    @test positions_equal(a, b)
    @test !positions_equal(a, c)
    @test positions_equal(a, c; tol=1.0)
end

@testset "window_shape accessor" begin
    fib_p = generate_fibonacci_projection(50)
    @test window_shape(fib_p) === :interval_1d

    fib_s = generate_fibonacci_substitution(4)
    @test window_shape(fib_s) === :none

    ab_p = generate_ammann_beenker_projection(2.0)
    @test window_shape(ab_p) === :box_2d

    ab_s = generate_ammann_beenker_substitution(2)
    @test window_shape(ab_s) === :none

    pen_p = generate_penrose_projection(2.0)
    @test window_shape(pen_p) === :box_3d

    pen_s = generate_penrose_substitution(2)
    @test window_shape(pen_s) === :none

    # Absent key falls back to :unknown.
    pen_p.parameters[:window_shape] = :unknown
    delete!(pen_p.parameters, :window_shape)
    @test window_shape(pen_p) === :unknown
end
