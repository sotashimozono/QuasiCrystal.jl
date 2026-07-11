using QuasiCrystal
using LatticeCore
using LinearAlgebra
using Test

function fib_chain_dsf(gen)
    qc = generate_fibonacci_substitution(gen)
    build_nearest_neighbor_bonds!(qc; cutoff=1.8)
    return qc
end

# Independent first moment: ∫ ω A(k,ω) dω = ⟨k|H|k⟩ = -(2t/N) Σ_bonds cos(k·Δr).
function first_moment_exact(qc, k, t)
    N = num_sites(qc)
    s = 0.0
    for b in bonds(qc)
        s += cos(sum(k[d] * b.vector[d] for d in eachindex(k)))
    end
    return -2t / N * s
end

@testset "dynamic structure factor A(k,ω)" begin
    qc = fib_chain_dsf(7)
    N = num_sites(qc)
    t = 1.0
    kpoints = [[0.0], [0.4], [1.0], [Float64(π)]]

    omegas, A = dynamic_structure_factor(qc, kpoints; t=t, nomega=600, broadening=0.03)
    dω = omegas[2] - omegas[1]

    @testset "shape & non-negativity" begin
        @test size(A) == (length(kpoints), length(omegas))
        @test all(A .≥ 0)
    end

    @testset "spectral-weight sum rule: ∫A(k,ω)dω = 1" begin
        for j in eachindex(kpoints)
            @test sum(@view A[j, :]) * dω ≈ 1.0 rtol = 1e-3
        end
    end

    @testset "first moment: ∫ω A(k,ω)dω = ⟨k|H|k⟩" begin
        for (j, k) in enumerate(kpoints)
            m1 = sum(omegas[m] * A[j, m] for m in eachindex(omegas)) * dω
            @test m1 ≈ first_moment_exact(qc, k, t) atol = 2e-3
        end
        # k = 0 special case: ⟨0|H|0⟩ = -(2t/N)·num_bonds
        m1_gamma = sum(omegas[m] * A[1, m] for m in eachindex(omegas)) * dω
        @test m1_gamma ≈ -2t / N * length(bonds(qc)) atol = 2e-3
    end

    @testset "momentum-lattice method matches k-vector method" begin
        # a BraggPeakSet as k-source; both routes must agree on the same grid.
        peaks = bragg_peaks(qc; kmax=6.0, intensity_cutoff=1e-6)
        kvecs = [k_point(peaks, j) for j in 1:num_k_points(peaks)]
        o1, A1 = dynamic_structure_factor(qc, peaks; t=t, nomega=120, broadening=0.05)
        o2, A2 = dynamic_structure_factor(qc, kvecs; t=t, nomega=120, broadening=0.05)
        @test o1 == o2
        @test A1 ≈ A2
    end

    @testset "argument validation" begin
        @test_throws ArgumentError dynamic_structure_factor(qc, kpoints; nomega=0)
        @test_throws ArgumentError dynamic_structure_factor(qc, kpoints; broadening=0.0)
    end
end
