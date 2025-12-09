@testset "Fibonacci Lattice Tests" begin
  using LinearAlgebra

  @testset "Projection Method" begin
    # Test basic generation
    qc_data = generate_fibonacci_projection(10)
    @test length(qc_data.positions) <= 10
    @test length(qc_data.positions) > 0
    @test qc_data.generation_method isa ProjectionMethod
    @test haskey(qc_data.parameters, :n_points)

    # Test positions are sorted
    positions_1d = [p[1] for p in qc_data.positions]
    @test issorted(positions_1d)

    # Test all positions are positive or zero
    @test all(positions_1d .>= 0)

    # Test larger lattice
    qc_data_large = generate_fibonacci_projection(50)
    @test length(qc_data_large.positions) <= 50
    @test length(qc_data_large.positions) > length(qc_data.positions)
  end

  @testset "Substitution Method" begin
    # Test basic generation
    qc_data = generate_fibonacci_substitution(5)
    @test length(qc_data.positions) > 0
    @test qc_data.generation_method isa SubstitutionMethod
    @test haskey(qc_data.parameters, :generations)
    @test qc_data.parameters[:generations] == 5

    # Test sequence length follows Fibonacci numbers
    for n in 1:7
      qc_data_n = generate_fibonacci_substitution(n)
      expected_length = fibonacci_sequence_length(n)
      # +1 because we include position 0
      @test length(qc_data_n.positions) == expected_length + 1
    end

    # Test spacing ratio approaches golden ratio
    qc_data = generate_fibonacci_substitution(10)
    positions_1d = [p[1] for p in qc_data.positions]
    spacings = diff(positions_1d)
    long_spacings = filter(x -> x > 1.5, spacings)
    short_spacings = filter(x -> x <= 1.5, spacings)

    if !isempty(long_spacings) && !isempty(short_spacings)
      avg_long = sum(long_spacings) / length(long_spacings)
      avg_short = sum(short_spacings) / length(short_spacings)
      ratio = avg_long / avg_short
      # Should be close to golden ratio
      @test abs(ratio - ϕ) < 0.1
    end
  end

  @testset "Fibonacci Sequence Length" begin
    @test fibonacci_sequence_length(0) == 1
    @test fibonacci_sequence_length(1) == 2
    @test fibonacci_sequence_length(2) == 3
    @test fibonacci_sequence_length(3) == 5
    @test fibonacci_sequence_length(4) == 8
    @test fibonacci_sequence_length(5) == 13
    @test fibonacci_sequence_length(6) == 21
  end

  @testset "Method Comparison" begin
    # Both methods should produce similar number of points
    proj_data = generate_fibonacci_projection(20)
    subst_data = generate_fibonacci_substitution(6)  # Fibonacci(6) = 21 points

    # Should have similar number of points
    @test abs(length(proj_data.positions) - length(subst_data.positions)) < 5
  end
end
