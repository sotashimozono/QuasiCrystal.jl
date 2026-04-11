@testset "Fibonacci Fourier analysis" begin
    @testset "IntervalWindow Fourier transform" begin
        w = IntervalWindow(0.5)

        # q = 0 gives the total integral 2 * half_width.
        @test window_fourier(w, 0.0) ≈ 1.0
        @test window_fourier(w, SVector(0.0)) ≈ 1.0

        # Symmetric in q.
        @test window_fourier(w, 3.0) ≈ window_fourier(w, -3.0)

        # Matches the analytic formula 2 sin(qa) / q.
        for q in (0.3, 1.1, 2.7, 5.4)
            expected = 2 * sin(q * 0.5) / q
            @test isapprox(window_fourier(w, q), expected; atol=1e-12)
        end
    end

    @testset "hyper_reciprocal_lattice for Fibonacci" begin
        qc = generate_fibonacci_projection(20)
        hrl = hyper_reciprocal_lattice(qc)

        # Parameter signature: (DPhys=1, DHyper=2, DPerp=1, T=Float64, W=IntervalWindow)
        @test hrl isa HyperReciprocalLattice{1,2,1,Float64,<:IntervalWindow}

        # Hyper reciprocal basis of Z² is 2π * I₂.
        @test hrl.hyper_basis ≈ SMatrix{2,2}(2π, 0.0, 0.0, 2π)

        # parallel_proj and perp_proj are orthogonal unit 1×2 rows.
        par = hrl.parallel_proj
        perp = hrl.perp_proj
        @test size(par) == (1, 2)
        @test size(perp) == (1, 2)
        @test sum(abs2, par) ≈ 1.0 atol = 1e-10
        @test sum(abs2, perp) ≈ 1.0 atol = 1e-10
        @test abs(sum(par .* perp)) < 1e-10   # orthogonal
    end

    @testset "bragg_peaks returns a BraggPeakSet subtype of AbstractMomentumLattice" begin
        qc = generate_fibonacci_projection(20)
        peaks = bragg_peaks(qc; kmax=20.0)

        @test peaks isa BraggPeakSet{1,2,Float64}
        @test peaks isa AbstractMomentumLattice{1,Float64}
        @test num_k_points(peaks) > 0

        # Γ (hyper index (0, 0)) is always present — amplitude = 2 * half_width
        γ_idx = findfirst(==((0, 0)), peaks.hyper_indices)
        @test γ_idx !== nothing
        @test peaks.peaks[γ_idx] ≈ SVector(0.0)

        # Intensities are non-negative.
        @test all(I -> I >= -1e-12, peaks.intensities)

        # Γ has the largest intensity.
        γ_I = peaks.intensities[γ_idx]
        @test γ_I == maximum(peaks.intensities)
    end

    @testset "dominant peaks cluster around golden-ratio momenta" begin
        # The (1, 1) and (1, 0) hyper indices are among the brightest
        # Bragg peaks of the Fibonacci chain; at kmax = 20 they must
        # both be present and both non-trivial.
        qc = generate_fibonacci_projection(20)
        peaks = bragg_peaks(qc; kmax=20.0)

        indices = Set(peaks.hyper_indices)
        @test (1, 1) in indices
        @test (1, 0) in indices
        @test (0, 1) in indices
        @test (-1, 0) in indices
        @test (0, -1) in indices

        # All of these should have positive intensity (they are not
        # filtered out by the default cutoff).
        for hi in [(1, 1), (1, 0), (0, 1), (-1, 0), (0, -1)]
            i = findfirst(==(hi), peaks.hyper_indices)
            @test peaks.intensities[i] > 0
        end
    end

    @testset "intensity_cutoff filters small peaks" begin
        qc = generate_fibonacci_projection(20)
        all_peaks = bragg_peaks(qc; kmax=20.0, intensity_cutoff=0.0)
        bright_peaks = bragg_peaks(qc; kmax=20.0, intensity_cutoff=0.01)

        @test num_k_points(bright_peaks) <= num_k_points(all_peaks)
        @test all(I >= 0.01 for I in bright_peaks.intensities)
    end

    @testset "fourier_module dispatch via trait" begin
        qc = generate_fibonacci_projection(15)

        # The momentum_lattice helper dispatches on reciprocal_support,
        # which for a QuasicrystalData is HasFourierModule(). It
        # should route through fourier_module(qc) -> bragg_peaks(qc)
        # with the default cutoff.
        @test reciprocal_support(qc) isa HasFourierModule

        ml = momentum_lattice(qc)
        @test ml isa BraggPeakSet
        @test num_k_points(ml) > 0

        # Direct call with an explicit cutoff also works.
        ml_custom = fourier_module(qc; kmax=15.0)
        @test ml_custom isa BraggPeakSet
    end

    @testset "BraggPeakSet composes with structure_factor" begin
        qc = generate_fibonacci_projection(20)
        peaks = bragg_peaks(qc; kmax=10.0)

        # Run structure_factor on a trivial (all-1) state. The
        # ferromagnetic structure factor at k = 0 should be N.
        state = ones(Int8, num_sites(qc))
        S0 = structure_factor(qc, state, SVector(0.0))
        @test S0 ≈ Float64(num_sites(qc)) atol = 1e-10

        # And iterating structure_factor over the whole BraggPeakSet
        # returns one real value per peak.
        Sks = structure_factor(qc, state, peaks)
        @test Sks isa Vector{Float64}
        @test length(Sks) == num_k_points(peaks)
    end
end
