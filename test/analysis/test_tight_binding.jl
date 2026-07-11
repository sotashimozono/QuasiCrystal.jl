using QuasiCrystal
using LatticeCore
using LinearAlgebra
using SparseArrays
using Test

# Helper: a Fibonacci chain with populated nearest-neighbour bonds. Its NN
# graph is a simple path, so uniform-hopping tight binding has the exact
# particle-in-a-box spectrum -2t·cos(kπ/(N+1)).
function fib_chain(gen)
    qc = generate_fibonacci_substitution(gen)
    build_nearest_neighbor_bonds!(qc; cutoff=1.8)
    return qc
end

@testset "tight-binding & spectral analysis" begin
    @testset "Hamiltonian structural invariants" begin
        qc = fib_chain(7)
        N = num_sites(qc)
        nb = length(bonds(qc))
        @test nb == N - 1                       # path graph
        t = 1.3
        H = tight_binding_hamiltonian(qc; t=t)
        Hd = Matrix(H)
        @test issymmetric(Hd)
        @test all(iszero, diag(Hd))             # no on-site term
        @test nnz(H) == 2nb                      # each bond -> 2 off-diagonal entries
        @test tr(Hd) ≈ 0 atol = 1e-12
        # Tr(H²) = Σ_ij |H_ij|² = 2·nb·t²
        @test tr(Hd * Hd) ≈ 2nb * t^2 rtol = 1e-12
    end

    @testset "path-graph closed-form spectrum (Fibonacci chain)" begin
        for gen in (5, 7)
            qc = fib_chain(gen)
            N = num_sites(qc)
            t = 0.9
            E = spectrum(qc; t=t)
            exact = sort([-2t * cos(k * π / (N + 1)) for k in 1:N])
            @test E ≈ exact rtol = 1e-10
            # spectrum is symmetric about 0 for a bipartite (path) graph
            @test E ≈ -reverse(E) rtol = 1e-10
        end
    end

    @testset "on-site energy shifts the spectrum rigidly" begin
        qc = fib_chain(6)
        ε = 0.7
        @test spectrum(qc; t=1.0, onsite=ε) ≈ spectrum(qc; t=1.0) .+ ε rtol = 1e-10
        # per-site on-site vector of constant ε agrees with the scalar form
        N = num_sites(qc)
        @test spectrum(qc; onsite=fill(ε, N)) ≈ spectrum(qc; onsite=ε) rtol = 1e-10
    end

    @testset "per-bond hoppings" begin
        qc = fib_chain(6)
        nb = length(bonds(qc))
        # a uniform hopping vector reproduces the scalar-t Hamiltonian
        H_vec = tight_binding_hamiltonian(qc, fill(2.0, nb))
        H_sca = tight_binding_hamiltonian(qc; t=2.0)
        @test Matrix(H_vec) == Matrix(H_sca)
        # dimension mismatch is rejected
        @test_throws DimensionMismatch tight_binding_hamiltonian(qc, fill(1.0, nb + 1))
        # a bond-modulated chain (JL/JS) opens a gap: distinct hoppings split
        # the band, so the density of states at E=0 collapses. Here we just
        # check the spectrum changes and stays symmetric.
        hop = [isodd(k) ? 1.0 : 0.4 for k in 1:nb]
        E = eigvals(Symmetric(Matrix(tight_binding_hamiltonian(qc, hop))))
        @test E ≈ -reverse(E) rtol = 1e-10
        @test !(E ≈ spectrum(qc; t=1.0))
    end

    @testset "inverse participation ratio" begin
        N = 40
        # fully extended state -> IPR = 1/N
        ext = fill(1 / sqrt(N), N)
        @test inverse_participation_ratio(ext) ≈ 1 / N rtol = 1e-12
        # single-site state -> IPR = 1
        loc = zeros(N);
        loc[7] = 1.0
        @test inverse_participation_ratio(loc) == 1.0
        # a state on m equal sites -> IPR = 1/m (participation number = m)
        m = 5;
        s = zeros(N);
        s[1:m] .= 1 / sqrt(m)
        @test inverse_participation_ratio(s) ≈ 1 / m rtol = 1e-12
        @test_throws ArgumentError inverse_participation_ratio(zeros(N))

        # per-eigenstate IPRs are bounded in [1/N, 1]
        qc = fib_chain(6)
        Nq = num_sites(qc)
        _, iprs = inverse_participation_ratios(qc; t=1.0)
        @test length(iprs) == Nq
        @test all(1 / Nq - 1e-9 .≤ iprs .≤ 1 + 1e-9)
    end

    @testset "density of states normalization" begin
        qc = fib_chain(7)
        N = num_sites(qc)
        # histogram: integrates to the number of states
        centers, dos = density_of_states(qc; t=1.0, nbins=40)
        dE = centers[2] - centers[1]
        @test sum(dos) * dE ≈ N rtol = 1e-10
        @test length(centers) == length(dos) == 40
        # Gaussian-broadened DOS also integrates to ≈ N
        c2, dos2 = density_of_states(qc; t=1.0, nbins=200, broadening=0.05)
        dE2 = c2[2] - c2[1]
        @test sum(dos2) * dE2 ≈ N rtol = 5e-2
        @test_throws ArgumentError density_of_states(qc; nbins=0)
    end

    @testset "2D Penrose invariants" begin
        qc = generate_penrose_substitution(3)
        build_nearest_neighbor_bonds!(qc; cutoff=1.1)
        nb = length(bonds(qc))
        @test nb > 0
        t = 1.0
        Hd = Matrix(tight_binding_hamiltonian(qc; t=t))
        @test issymmetric(Hd)
        E = eigvals(Symmetric(Hd))
        @test sum(E) ≈ 0 atol = 1e-9                       # Tr H = 0
        @test sum(abs2, E) ≈ 2nb * t^2 rtol = 1e-10        # Tr H² = 2·nb·t²
        # Gershgorin: |E| ≤ t · max coordination
        maxdeg = maximum(length(neighbors(qc, i)) for i in 1:num_sites(qc))
        @test maximum(abs, E) ≤ t * maxdeg + 1e-9
    end
end
