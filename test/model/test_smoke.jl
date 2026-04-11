@testset "smoke test: Fibonacci projection builds and satisfies AbstractLattice" begin
    qc = generate_fibonacci_projection(20)

    # QuasicrystalData is a concrete LatticeCore.AbstractLattice
    @test qc isa AbstractLattice{1,Float64}
    @test num_sites(qc) > 0
    @test num_sites(qc) <= 20

    # Positions are SVectors now, sorted and non-negative
    @test position(qc, 1) isa SVector{1,Float64}
    positions_1d = [position(qc, i)[1] for i in 1:num_sites(qc)]
    @test issorted(positions_1d)
    @test all(positions_1d .>= 0)

    # Traits
    @test periodicity(qc) isa Aperiodic
    @test reciprocal_support(qc) isa HasFourierModule
    @test is_finite(qc) == true

    # Connectivity is empty until build_nearest_neighbor_bonds! runs
    @test isempty(bonds(qc))
    @test all(isempty(neighbors(qc, i)) for i in 1:num_sites(qc))

    # Run the distance-based bond builder; every close pair should
    # now show up in the bond list.
    build_nearest_neighbor_bonds!(qc; cutoff=ϕ + 0.5)
    @test !isempty(bonds(qc))
    for b in bonds(qc)
        @test b isa Bond{1,Float64}
        @test 1 <= b.i <= num_sites(qc)
        @test 1 <= b.j <= num_sites(qc)
        @test norm(b.vector) < ϕ + 0.5
    end

    # Default layout is a uniform Ising layout
    @test site_layout(qc) isa UniformLayout
    @test site_type(qc, 1) isa IsingSite
end
