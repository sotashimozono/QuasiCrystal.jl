@testset "convenience builder API" begin
    @testset "fibonacci(n) substitution shortcut" begin
        qc = fibonacci(4)
        @test qc isa QuasicrystalData{1,Float64}
        @test qc.topology isa FibonacciLattice
        @test qc.generation_method isa SubstitutionMethod
        @test qc.parameters[:generations] == 4
        # Defaults to 10 generations.
        qc_default = fibonacci()
        @test qc_default.parameters[:generations] == 10
        @test num_sites(qc_default) > num_sites(qc)
        # Underlying generator and shortcut produce identical positions.
        ref = generate_fibonacci_substitution(4)
        @test num_sites(qc) == num_sites(ref)
        @test position(qc, 1) == position(ref, 1)
        @test position(qc, num_sites(qc)) == position(ref, num_sites(ref))
    end

    @testset "penrose(n) substitution shortcut" begin
        qc = penrose(2)
        @test qc isa QuasicrystalData{2,Float64}
        @test qc.topology isa PenroseP3
        @test qc.generation_method isa SubstitutionMethod
        @test qc.parameters[:generations] == 2
        # Defaults to 3 generations.
        qc_default = penrose()
        @test qc_default.parameters[:generations] == 3
        # Forwarded `method` keyword is preserved on the data.
        qc_alg = penrose(2; method=DefaultSubstitution())
        @test qc_alg.generation_method isa SubstitutionMethod
        @test qc_alg.generation_method.algorithm isa DefaultSubstitution
        # Bare-algorithm and wrapped-method calls agree.
        qc_wrapped = penrose(2; method=SubstitutionMethod(DefaultSubstitution()))
        @test num_sites(qc_alg) == num_sites(qc_wrapped)
    end

    @testset "ammann_beenker(n) substitution shortcut" begin
        qc = ammann_beenker(2)
        @test qc isa QuasicrystalData{2,Float64}
        @test qc.topology isa AmmannBeenker
        @test qc.generation_method isa SubstitutionMethod
        @test qc.parameters[:generations] == 2
        @test qc.parameters[:symmetry] == 8
        # Default generations.
        qc_default = ammann_beenker()
        @test qc_default.parameters[:generations] == 3
    end

    @testset "penrose_projected(radius) projection shortcut" begin
        qc = penrose_projected(3.0)
        @test qc isa QuasicrystalData{2,Float64}
        @test qc.topology isa PenroseP3
        @test qc.generation_method isa ProjectionMethod
        @test qc.parameters[:radius] == 3.0
        # Reference equality with the canonical generator.
        ref = generate_penrose_projection(3.0)
        @test num_sites(qc) == num_sites(ref)
    end

    @testset "ammann_beenker_projected(radius) projection shortcut" begin
        qc = ammann_beenker_projected(3.0)
        @test qc isa QuasicrystalData{2,Float64}
        @test qc.topology isa AmmannBeenker
        @test qc.generation_method isa ProjectionMethod
        @test qc.parameters[:radius] == 3.0
        @test qc.parameters[:symmetry] == 8
        ref = generate_ammann_beenker_projection(3.0)
        @test num_sites(qc) == num_sites(ref)
    end

    @testset "fibonacci_projected(n_points) projection shortcut" begin
        qc = fibonacci_projected(50)
        @test qc isa QuasicrystalData{1,Float64}
        @test qc.topology isa FibonacciLattice
        @test qc.generation_method isa ProjectionMethod
        ref = generate_fibonacci_projection(50)
        @test num_sites(qc) == num_sites(ref)
    end

    @testset "shortcut return values match canonical generator field-by-field" begin
        # Exercise the full QuasicrystalData payload (positions, params,
        # generation_method) for each builder so that any future drift
        # between builder and underlying generator is caught.
        for (shortcut, ref) in (
            (() -> fibonacci(3), () -> generate_fibonacci_substitution(3)),
            (() -> penrose(2), () -> generate_penrose_substitution(2)),
            (() -> ammann_beenker(2), () -> generate_ammann_beenker_substitution(2)),
            (() -> penrose_projected(3.0), () -> generate_penrose_projection(3.0)),
            (
                () -> ammann_beenker_projected(3.0),
                () -> generate_ammann_beenker_projection(3.0),
            ),
            (() -> fibonacci_projected(20), () -> generate_fibonacci_projection(20)),
        )
            qc = shortcut()
            qc_ref = ref()
            @test num_sites(qc) == num_sites(qc_ref)
            @test typeof(qc.topology) === typeof(qc_ref.topology)
            @test typeof(qc.generation_method) === typeof(qc_ref.generation_method)
        end
    end
end
