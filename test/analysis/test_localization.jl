using QuasiCrystal
using LatticeCore
using LinearAlgebra
using Random
using Test

function fib_chain_loc(gen)
    qc = generate_fibonacci_substitution(gen)
    build_nearest_neighbor_bonds!(qc; cutoff=1.8)
    return qc
end

@testset "localization diagnostics" begin
    @testset "mean IPR: matrix vs lattice method, path-graph value" begin
        qc = fib_chain_loc(8)
        N = num_sites(qc)
        H = tight_binding_hamiltonian(qc; t=1.0)
        @test mean_inverse_participation_ratio(H) ==
            mean_inverse_participation_ratio(qc; t=1.0)
        # extended path-graph eigenstates have mean IPR ≈ 1.5/(N+1)
        @test mean_inverse_participation_ratio(qc) ≈ 1.5 / (N + 1) rtol = 0.15
    end

    @testset "ipr_scaling_exponent recovers a known power law" begin
        # synthetic data mean_ipr = C·N^{-0.7} must fit τ = 0.7
        Ns = [20, 40, 80, 160, 320]
        mis = [3.0 * N^(-0.7) for N in Ns]
        @test ipr_scaling_exponent(Ns, mis) ≈ 0.7 rtol = 1e-10
        # input validation
        @test_throws ArgumentError ipr_scaling_exponent([10], [0.1])
        @test_throws DimensionMismatch ipr_scaling_exponent([10, 20], [0.1])
        @test_throws ArgumentError ipr_scaling_exponent([10, 20], [0.1, -0.2])
    end

    # --- three physical regimes, distinguished by the exponent τ ---
    gens = 7:11

    @testset "extended: uniform hopping ⇒ τ ≈ 1" begin
        qcs = [fib_chain_loc(g) for g in gens]
        sizes, mis = ipr_scaling(qcs; t=1.0)
        τ_ext = ipr_scaling_exponent(sizes, mis)
        @test τ_ext > 0.9
    end

    @testset "localized: strong on-site disorder ⇒ τ ≈ 0" begin
        sizes = Int[];
        mis = Float64[]
        for g in gens
            qc = fib_chain_loc(g)
            N = num_sites(qc)
            Random.seed!(1234 + g)
            onsite = 12.0 .* (rand(N) .- 0.5)      # W = 12 box disorder
            push!(sizes, N)
            push!(mis, mean_inverse_participation_ratio(qc; t=1.0, onsite=onsite))
        end
        τ_loc = ipr_scaling_exponent(sizes, mis)
        @test τ_loc < 0.2
    end

    @testset "critical: off-diagonal Fibonacci model ⇒ 0 < τ < 1" begin
        sizes = Int[];
        mis = Float64[]
        for g in gens
            qc = fib_chain_loc(g)
            nb = length(bonds(qc))
            word = fibonacci_word(FibonacciLattice(), g)   # length == nb
            hop = [word[k] == 0 ? 1.0 : 0.4 for k in 1:nb] # tL=1.0, tS=0.4
            H = tight_binding_hamiltonian(qc, hop)
            push!(sizes, num_sites(qc))
            push!(mis, mean_inverse_participation_ratio(H))
        end
        τ_crit = ipr_scaling_exponent(sizes, mis)
        @test 0.3 < τ_crit < 0.9
    end

    @testset "regime ordering τ_localized < τ_critical < τ_extended" begin
        # recompute compactly and assert the physical ordering in one place
        τ = Dict{Symbol,Float64}()
        # extended
        let qcs = [fib_chain_loc(g) for g in gens]
            s, m = ipr_scaling(qcs; t=1.0)
            τ[:ext] = ipr_scaling_exponent(s, m)
        end
        # critical
        let s = Int[], m = Float64[]
            for g in gens
                qc = fib_chain_loc(g);
                nb = length(bonds(qc))
                word = fibonacci_word(FibonacciLattice(), g)
                H = tight_binding_hamiltonian(qc, [word[k] == 0 ? 1.0 : 0.4 for k in 1:nb])
                push!(s, num_sites(qc));
                push!(m, mean_inverse_participation_ratio(H))
            end
            τ[:crit] = ipr_scaling_exponent(s, m)
        end
        # localized
        let s = Int[], m = Float64[]
            for g in gens
                qc = fib_chain_loc(g);
                N = num_sites(qc)
                Random.seed!(999 + g)
                push!(s, N)
                push!(
                    m, mean_inverse_participation_ratio(qc; onsite=12.0 .* (rand(N) .- 0.5))
                )
            end
            τ[:loc] = ipr_scaling_exponent(s, m)
        end
        @test τ[:loc] < τ[:crit] < τ[:ext]
    end

    @testset "energy-window averaging" begin
        qc = fib_chain_loc(8)
        H = tight_binding_hamiltonian(qc; t=1.0)
        # a window around the band centre selects a subset and returns a valid IPR
        val = mean_inverse_participation_ratio(H; energy_window=(-0.5, 0.5))
        @test 0 < val ≤ 1
        # an out-of-band window is empty and must error
        @test_throws ArgumentError mean_inverse_participation_ratio(
            H; energy_window=(100.0, 200.0)
        )
    end
end
