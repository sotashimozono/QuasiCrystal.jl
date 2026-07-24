using QuasiCrystal
using Test

const _QC_FAMILIES = (FibonacciLattice(), PenroseP3(), AmmannBeenker())

@testset "InfiniteQuasicrystal — traits and undefined site API" begin
    for topo in _QC_FAMILIES
        inf = InfiniteQuasicrystal(topo)

        @test size_trait(inf) == InfiniteSize()
        @test !is_finite(inf)
        @test periodicity(inf) == Aperiodic()
        @test reciprocal_support(inf) == HasFourierModule()
        @test topology(inf) == TopologyTrait{:quasiperiodic}()

        # No sites exist until materialised.
        @test_throws DomainError num_sites(inf)
        @test_throws DomainError position(inf, 1)
        @test_throws DomainError neighbors(inf, 1)
        # The Monte Carlo entry guard must reject the infinite object.
        @test_throws ArgumentError require_finite(inf)
    end
end

@testset "materialize matches the direct projection generator" begin
    # Independent oracle: the lazy materialize path must reproduce, site
    # for site, what the underlying generator produces directly.
    fib = materialize(InfiniteQuasicrystal(FibonacciLattice()); n_points=50)
    @test fib isa QuasicrystalData
    @test is_finite(fib)
    @test num_sites(fib) == 50
    @test positions(fib) == positions(generate_fibonacci_projection(50))
    @test require_finite(fib) === nothing   # finite ⇒ no throw

    pen = materialize(InfiniteQuasicrystal(PenroseP3()); radius=4.0)
    @test pen isa QuasicrystalData
    @test num_sites(pen) == num_sites(generate_penrose_projection(4.0))
    @test positions(pen) == positions(generate_penrose_projection(4.0))

    ab = materialize(InfiniteQuasicrystal(AmmannBeenker()); radius=4.0)
    @test ab isa QuasicrystalData
    @test num_sites(ab) == num_sites(generate_ammann_beenker_projection(4.0))
    @test positions(ab) == positions(generate_ammann_beenker_projection(4.0))
end

@testset "materialize propagates the site layout" begin
    inf = InfiniteQuasicrystal(FibonacciLattice(); layout=UniformLayout(XYSite()))
    # UniformLayout ignores the index, so site_type resolves even on the
    # infinite object.
    @test site_type(inf, 1) == XYSite()

    d = materialize(inf; n_points=20)
    @test site_layout(d) === site_layout(inf)   # same layout object carried through
    @test site_type(d, 1) == XYSite()
    @test num_sites(d) == 20
end

@testset "materialize default cutoffs produce non-empty lattices" begin
    for topo in _QC_FAMILIES
        d = materialize(InfiniteQuasicrystal(topo))
        @test d isa QuasicrystalData
        @test is_finite(d)
        @test num_sites(d) > 0
    end
end
