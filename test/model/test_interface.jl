@testset "AbstractLattice Interface Tests" begin
    using LinearAlgebra

    @testset "Bond and Connection Types" begin
        # Test Bond construction
        bond = Bond(1, 2, 1, [1.0, 0.0])
        @test bond.src == 1
        @test bond.dst == 2
        @test bond.type == 1
        @test bond.vector == [1.0, 0.0]
        @test bond isa AbstractLatticeConnection

        # Test Connection construction
        conn = Connection(1, 1, 1, 0, 1)
        @test conn.src_sub == 1
        @test conn.dst_sub == 1
        @test conn.dx == 1
        @test conn.dy == 0
        @test conn.type == 1
        @test conn isa AbstractLatticeConnection
    end

    @testset "AbstractQuasicrystal inherits AbstractLattice" begin
        # Test type hierarchy
        @test FibonacciLattice <: AbstractQuasicrystal{1}
        @test AbstractQuasicrystal{1} <: AbstractLattice{1}
        @test PenroseP3 <: AbstractQuasicrystal{2}
        @test AbstractQuasicrystal{2} <: AbstractLattice{2}
        @test AmmannBeenker <: AbstractQuasicrystal{2}
    end

    @testset "QuasicrystalData with bonds" begin
        # Generate a simple quasicrystal
        qc_data = generate_fibonacci_projection(10)

        # Test initial state
        @test num_sites(qc_data) == length(qc_data.positions)
        @test num_bonds(qc_data) == 0
        @test length(qc_data.nearest_neighbors) == num_sites(qc_data)
        @test all(isempty.(qc_data.nearest_neighbors))

        # Build nearest neighbor bonds
        build_nearest_neighbor_bonds!(qc_data; cutoff=2.0)

        # Test bonds were created
        @test num_bonds(qc_data) > 0
        @test length(get_bonds(qc_data)) == num_bonds(qc_data)

        # Test that bonds are valid
        for bond in qc_data.bonds
            @test 1 <= bond.src <= num_sites(qc_data)
            @test 1 <= bond.dst <= num_sites(qc_data)
            @test bond.src < bond.dst  # We only create bonds in one direction
            @test length(bond.vector) == 1  # 1D lattice

            # Check bond vector is correct
            expected_vector = qc_data.positions[bond.dst] - qc_data.positions[bond.src]
            @test bond.vector ≈ expected_vector
        end

        # Test nearest neighbors
        nn = get_nearest_neighbors(qc_data)
        @test length(nn) == num_sites(qc_data)

        # All neighbors should be valid site indices
        for (i, neighbors) in enumerate(nn)
            @test all(1 .<= neighbors .<= num_sites(qc_data))
        end
    end

    @testset "Interface Methods" begin
        # Test with Fibonacci lattice
        qc_data = generate_fibonacci_projection(20)
        build_nearest_neighbor_bonds!(qc_data; cutoff=2.0)

        # Test get_positions
        positions = get_positions(qc_data)
        @test positions isa Vector{Vector{Float64}}
        @test length(positions) == num_sites(qc_data)
        @test all(length.(positions) .== 1)  # 1D positions

        # Test get_bonds
        bonds = get_bonds(qc_data)
        @test bonds isa Vector{Bond}
        @test length(bonds) == num_bonds(qc_data)

        # Test get_nearest_neighbors
        nn = get_nearest_neighbors(qc_data)
        @test nn isa Vector{Vector{Int}}
        @test length(nn) == num_sites(qc_data)

        # Test num_sites and num_bonds
        @test num_sites(qc_data) isa Int
        @test num_bonds(qc_data) isa Int
        @test num_sites(qc_data) > 0
        @test num_bonds(qc_data) > 0
    end

    @testset "2D Quasicrystal Bonds" begin
        # Test with Penrose tiling
        qc_data = generate_penrose_projection(3.0)

        # Initially no bonds
        @test num_bonds(qc_data) == 0

        # Build bonds
        build_nearest_neighbor_bonds!(qc_data; cutoff=1.5)

        # Should have bonds now
        @test num_bonds(qc_data) > 0

        # Test bond properties
        for bond in qc_data.bonds
            @test length(bond.vector) == 2  # 2D
            @test norm(bond.vector) < 1.5   # Within cutoff
        end

        # Test symmetry: if i is neighbor of j, j should be neighbor of i
        nn = get_nearest_neighbors(qc_data)
        for i in 1:num_sites(qc_data)
            for j in nn[i]
                @test i in nn[j]
            end
        end
    end
end
