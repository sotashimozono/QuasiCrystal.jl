# Tests for QuasiCrystalNFFTExt — NUFFT-backed structure_factor on
# quasicrystals. Asserts:
#
# 1. The trait dispatch in LatticeCore routes `structure_factor(qc, state, ml)`
#    through our NUFFT extension (`QuasiCrystalNFFTExt`) once `using NFFT` has
#    been issued.
# 2. The NUFFT path agrees with the naive `O(N · M)` reference loop to
#    floating-point round-off, on small Penrose and Ammann–Beenker lattices,
#    for both real-valued spin-like states and complex-valued states.
# 3. A timing comparison on a larger Penrose runs to completion. We record
#    timings via `@elapsed` so CI logs surface the comparison; we do not
#    assert NUFFT < naive since `NFFT.NNDFTPlan` is mathematically the same
#    direct double loop and may be slower than the naive Julia path —
#    correctness is the load-bearing assertion here.

using QuasiCrystal
using LatticeCore
using NFFT
using LinearAlgebra
using StaticArrays
using Random
using Test

@testset "structure_factor NUFFT extension" begin
    @testset "trait dispatch through QuasiCrystalNFFTExt" begin
        # `using NFFT` should have triggered both LatticeCoreNFFTExt and
        # QuasiCrystalNFFTExt. We confirm by checking that the QC method
        # for `_structure_factor_fast` exists on `HasFourierModule` +
        # `QuasicrystalData`.
        ms = methods(LatticeCore._structure_factor_fast)
        qc_has_override = any(ms) do m
            sig = m.sig
            occursin("QuasicrystalData", string(sig)) &&
                occursin("HasFourierModule", string(sig))
        end
        @test qc_has_override
    end

    @testset "Penrose: naive vs NUFFT (real state)" begin
        Random.seed!(11)
        qc = generate_penrose_substitution(3)
        N = num_sites(qc)
        @test N > 100
        peaks = bragg_peaks(qc; kmax=8.0, intensity_cutoff=1e-6)
        @test num_k_points(peaks) > 10

        state = randn(N)
        S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        S_fast = structure_factor(qc, state, peaks)

        @test length(S_fast) == length(S_naive)
        # The NUFFT path goes through `NFFT.NNDFTPlan`, which is
        # mathematically the same sum as the naive loop. We expect
        # round-off-level agreement.
        absdiff = maximum(abs.(S_naive .- S_fast))
        scale = max(maximum(abs.(S_naive)), 1.0)
        @test absdiff / scale < 1e-10
    end

    @testset "Penrose: naive vs NUFFT (complex state)" begin
        Random.seed!(13)
        qc = generate_penrose_substitution(3)
        N = num_sites(qc)
        peaks = bragg_peaks(qc; kmax=8.0, intensity_cutoff=1e-6)

        # Random complex state of unit norm.
        state = randn(ComplexF64, N)
        state ./= norm(state)

        S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        S_fast = structure_factor(qc, state, peaks)

        absdiff = maximum(abs.(S_naive .- S_fast))
        scale = max(maximum(abs.(S_naive)), 1.0)
        @test absdiff / scale < 1e-10
    end

    @testset "Penrose: uniform state lights up Γ peak" begin
        # Sanity check that the structured path reproduces the basic
        # physics expectation: a uniform state |ψ⟩ = (1, 1, ..., 1)/√N
        # gives S(k = 0) = |Σ_n 1|² / N = N (after our normalisation).
        Random.seed!(0)
        qc = generate_penrose_substitution(3)
        N = num_sites(qc)
        peaks = bragg_peaks(qc; kmax=8.0, intensity_cutoff=1e-6)
        state = ones(Float64, N)

        S_fast = structure_factor(qc, state, peaks)
        S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        @test maximum(abs.(S_fast .- S_naive)) < 1e-10

        # Identify the Γ peak (k = 0). It may not be exactly the first
        # entry, depending on enumeration; pick the smallest |k|.
        k_norms = [norm(k_point(peaks, j)) for j in 1:num_k_points(peaks)]
        i_gamma = argmin(k_norms)
        @test k_norms[i_gamma] < 1e-10
        # |Σ_n 1|² / N = N
        @test S_fast[i_gamma] ≈ float(N) atol = 1e-8
    end

    @testset "Ammann–Beenker: naive vs NUFFT" begin
        Random.seed!(17)
        qc = generate_ammann_beenker_substitution(3)
        N = num_sites(qc)
        @test N > 100
        peaks = bragg_peaks(qc; kmax=8.0, intensity_cutoff=1e-6)
        @test num_k_points(peaks) > 10

        state = randn(ComplexF64, N)
        S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        S_fast = structure_factor(qc, state, peaks)

        absdiff = maximum(abs.(S_naive .- S_fast))
        scale = max(maximum(abs.(S_naive)), 1.0)
        @test absdiff / scale < 1e-10
    end

    @testset "Larger Penrose: completion + naive/NUFFT timing" begin
        # Larger problem (N > 1000) so the timings are visible. We do
        # not assert NUFFT < naive: NFFT.NNDFTPlan is the structured
        # direct sum and may be slower than the straight naive loop.
        # The assertion is correctness; the @elapsed values are logged
        # for the PR description / future profiling.
        Random.seed!(19)
        qc = generate_penrose_substitution(4)
        N = num_sites(qc)
        @test N > 1000
        peaks = bragg_peaks(qc; kmax=12.0, intensity_cutoff=1e-6)
        M = num_k_points(peaks)
        state = randn(ComplexF64, N)

        # Warm-up to skip first-call compilation cost.
        LatticeCore._structure_factor_naive(qc, state, peaks)
        structure_factor(qc, state, peaks)

        t_naive = @elapsed S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        t_nufft = @elapsed S_fast = structure_factor(qc, state, peaks)

        absdiff = maximum(abs.(S_naive .- S_fast))
        scale = max(maximum(abs.(S_naive)), 1.0)
        @test absdiff / scale < 1e-10

        @info "structure_factor benchmark" N M t_naive_ms = round(t_naive * 1000, digits=2) t_nufft_ms = round(t_nufft * 1000, digits=2) ratio = round(t_naive / max(t_nufft, eps()), digits=2)
    end
end
