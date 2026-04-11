@testset "QuasicrystalData AbstractLattice interface" begin
    @testset "AbstractQuasicrystal is a topology marker, not a lattice" begin
        # The post-migration convention: AbstractQuasicrystal{D} is a
        # dispatch key for `build_quasicrystal`, not a subtype of
        # `LatticeCore.AbstractLattice`.
        @test FibonacciLattice <: AbstractQuasicrystal{1}
        @test PenroseP3 <: AbstractQuasicrystal{2}
        @test AmmannBeenker <: AbstractQuasicrystal{2}
        @test !(FibonacciLattice <: AbstractLattice)
    end

    @testset "QuasicrystalData subtypes AbstractLattice{D, Float64}" begin
        fib = generate_fibonacci_projection(15)
        pen = generate_penrose_projection(3.0)
        ab = generate_ammann_beenker_projection(3.0)

        @test fib isa AbstractLattice{1,Float64}
        @test pen isa AbstractLattice{2,Float64}
        @test ab isa AbstractLattice{2,Float64}
    end

    @testset "LatticeCore traits on QuasicrystalData" begin
        qc = generate_fibonacci_substitution(4)
        @test periodicity(qc) isa Aperiodic
        @test reciprocal_support(qc) isa HasFourierModule
        @test is_finite(qc) == true
        @test size_trait(qc) isa FiniteSize{1}
        @test size_trait(qc).dims == (num_sites(qc),)
    end

    @testset "Bond and Connection are now LatticeCore types" begin
        # The old QuasiCrystal.Bond(src, dst, type::Int, Vector{Float64})
        # constructor is gone. Users construct LatticeCore.Bond{D, T}
        # directly, which is re-exported from QuasiCrystal.
        b = Bond{2,Float64}(1, 2, SVector(1.0, 0.0), :nearest)
        @test b.i == 1
        @test b.j == 2
        @test b.vector == SVector(1.0, 0.0)
        @test b.type === :nearest
    end

    @testset "build_quasicrystal dispatch" begin
        @test build_quasicrystal(FibonacciLattice; n_points=10) isa QuasicrystalData
        @test build_quasicrystal(PenroseP3; radius=3.0) isa QuasicrystalData
        @test build_quasicrystal(AmmannBeenker; radius=3.0) isa QuasicrystalData
        @test build_quasicrystal(
            FibonacciLattice; generator=:substitution, generations=5
        ) isa QuasicrystalData
    end

    @testset "QuasicrystalData re-exports LatticeCore vocabulary" begin
        @test PeriodicAxis() isa AbstractAxisBC
        @test IsingSite() isa AbstractSiteType
        @test UniformLayout(IsingSite()) isa AbstractSiteLayout
        @test RowMajor() isa AbstractIndexing
    end
end
