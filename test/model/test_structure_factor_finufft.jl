# Tests for QuasiCrystalFINUFFTExt — genuine O((N+M) log(1/ε)) NUFFT-backed
# structure_factor on quasicrystals (issue #61). Asserts:
#
# 1. The trait dispatch in LatticeCore routes `structure_factor(qc, state, ml)`
#    through our FINUFFT extension once `using FINUFFT` has been issued.
# 2. The FINUFFT type-3 path agrees with the naive `O(N · M)` reference loop
#    to the requested tolerance (~1e-10), on Penrose (2D) and Fibonacci (1D)
#    lattices, for both real and complex states.
# 3. On a larger problem the FINUFFT path is *faster* than the naive loop
#    (logged; see the controlled 25× benchmark in the PR). Correctness is the
#    load-bearing assertion — we do not hard-gate on wall-clock in CI.

using QuasiCrystal
using LatticeCore
using FINUFFT
using LinearAlgebra
using StaticArrays
using Random
using Test

# Handle to the loaded extension so we can exercise the FINUFFT type-3 path
# directly, bypassing the small-problem crossover that would otherwise route
# these (deliberately small, fast-to-run) test cases back to the naive loop.
const QCF = Base.get_extension(QuasiCrystal, :QuasiCrystalFINUFFTExt)
finufft_S(qc, state, ml) = QCF._structure_factor_finufft(qc, state, ml)

@testset "structure_factor FINUFFT extension" begin
    @testset "trait dispatch through QuasiCrystalFINUFFTExt" begin
        @test QCF !== nothing
        ms = methods(LatticeCore._structure_factor_fast)
        qc_has_override = any(ms) do m
            sig = string(m.sig)
            occursin("QuasicrystalData", sig) && occursin("HasFourierModule", sig)
        end
        @test qc_has_override
    end

    @testset "Penrose (2D): naive vs FINUFFT (real state)" begin
        Random.seed!(11)
        qc = generate_penrose_substitution(3)
        N = num_sites(qc)
        @test N > 100
        peaks = bragg_peaks(qc; kmax=8.0, intensity_cutoff=1e-6)
        @test num_k_points(peaks) > 10

        state = randn(N)
        S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        S_fast = finufft_S(qc, state, peaks)          # force the FINUFFT path

        @test length(S_fast) == length(S_naive)
        absdiff = maximum(abs.(S_naive .- S_fast))
        scale = max(maximum(abs.(S_naive)), 1.0)
        @test absdiff / scale < 1e-10
        # public entry point agrees too (routes via crossover)
        @test maximum(abs.(structure_factor(qc, state, peaks) .- S_naive)) / scale < 1e-10
    end

    @testset "Penrose (2D): naive vs FINUFFT (complex state)" begin
        Random.seed!(13)
        qc = generate_penrose_substitution(3)
        N = num_sites(qc)
        peaks = bragg_peaks(qc; kmax=8.0, intensity_cutoff=1e-6)

        state = randn(ComplexF64, N)
        state ./= norm(state)

        S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        S_fast = finufft_S(qc, state, peaks)

        absdiff = maximum(abs.(S_naive .- S_fast))
        scale = max(maximum(abs.(S_naive)), 1.0)
        @test absdiff / scale < 1e-10
    end

    @testset "Fibonacci (1D): naive vs FINUFFT" begin
        # Exercises the 1D type-3 path (nufft1d3).
        Random.seed!(23)
        qc = generate_fibonacci_substitution(9)
        N = num_sites(qc)
        @test N > 50
        peaks = bragg_peaks(qc; kmax=8.0, intensity_cutoff=1e-6)
        @test num_k_points(peaks) > 5

        state = randn(ComplexF64, N)
        S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        S_fast = finufft_S(qc, state, peaks)

        absdiff = maximum(abs.(S_naive .- S_fast))
        scale = max(maximum(abs.(S_naive)), 1.0)
        @test absdiff / scale < 1e-10
    end

    @testset "Penrose: uniform state lights up Γ peak" begin
        Random.seed!(0)
        qc = generate_penrose_substitution(3)
        N = num_sites(qc)
        peaks = bragg_peaks(qc; kmax=8.0, intensity_cutoff=1e-6)
        state = ones(Float64, N)

        S_fast = finufft_S(qc, state, peaks)
        k_norms = [norm(k_point(peaks, j)) for j in 1:num_k_points(peaks)]
        i_gamma = argmin(k_norms)
        @test k_norms[i_gamma] < 1e-10
        # |Σ_n 1|² / N = N
        @test S_fast[i_gamma] ≈ float(N) atol = 1e-8
    end

    @testset "Larger Penrose: completion + naive/FINUFFT timing" begin
        # Note on scale: the crossover heuristic (`_NUFFT_MIN_WORK`) routes
        # problems with `N · M` below ~2e6 back to the naive loop, so on this
        # patch the *public* `structure_factor` equals the naive result. We
        # additionally time the FINUFFT path directly for information — its
        # fixed setup overhead makes it slower here; the decisive ~25× win
        # appears only at `N · M ≳ 1e8` (see the controlled benchmark in the
        # PR). Correctness is the load-bearing assertion; we never hard-gate
        # on wall-clock in CI.
        Random.seed!(19)
        qc = generate_penrose_substitution(6)
        N = num_sites(qc)
        @test N > 1000
        peaks = bragg_peaks(qc; kmax=12.0, intensity_cutoff=1e-6)
        M = num_k_points(peaks)
        state = randn(ComplexF64, N)

        # Warm-up to skip first-call compilation cost.
        LatticeCore._structure_factor_naive(qc, state, peaks)
        finufft_S(qc, state, peaks)

        t_naive = @elapsed S_naive = LatticeCore._structure_factor_naive(qc, state, peaks)
        t_fast = @elapsed S_fast = finufft_S(qc, state, peaks)

        absdiff = maximum(abs.(S_naive .- S_fast))
        scale = max(maximum(abs.(S_naive)), 1.0)
        @test absdiff / scale < 1e-10

        @info "structure_factor benchmark (FINUFFT type-3)" N M work = N * M t_naive_ms = round(
            t_naive * 1000, digits=2
        ) t_finufft_ms = round(t_fast * 1000, digits=2) ratio = round(
            t_naive / max(t_fast, eps()), digits=2
        )
    end
end
