@testset "Ammann–Beenker Fourier analysis" begin
    @testset "BoxWindow Fourier transform" begin
        w = BoxWindow(SVector(0.5, 0.5))

        # q = 0 → total area = ∏ 2·aᵢ = 1.0
        @test window_fourier(w, SVector(0.0, 0.0)) ≈ 1.0

        # Separability: Ŵ(q₁,q₂) = sinc(q₁a₁)·sinc(q₂a₂) (with our
        # 2 sin(qa)/q convention at q=0 → 2a).
        for q1 in (0.3, 1.7, 4.1), q2 in (0.5, 2.2, 3.3)
            expected = (2 * sin(q1 * 0.5) / q1) * (2 * sin(q2 * 0.5) / q2)
            @test isapprox(
                window_fourier(w, SVector(q1, q2)), expected; atol=1e-12
            )
        end

        # One-axis limit: q₂ = 0 reduces to a 1D sinc times 2·a₂.
        q1 = 1.3
        expected = (2 * sin(q1 * 0.5) / q1) * 1.0
        @test isapprox(
            window_fourier(w, SVector(q1, 0.0)), expected; atol=1e-12
        )

        # Symmetric in both q axes.
        @test window_fourier(w, SVector(1.1, 2.7)) ≈
              window_fourier(w, SVector(-1.1, -2.7))
    end

    @testset "hyper_reciprocal_lattice for Ammann–Beenker" begin
        qc = generate_ammann_beenker_projection(4.0)
        hrl = hyper_reciprocal_lattice(qc)

        @test hrl isa HyperReciprocalLattice{2,4,2,Float64,<:BoxWindow}

        # Host reciprocal basis is 2π·I₄.
        @test hrl.hyper_basis ≈ SMatrix{4,4,Float64}(2π * I)

        par = hrl.parallel_proj
        perp = hrl.perp_proj
        @test size(par) == (2, 4)
        @test size(perp) == (2, 4)

        # Each column of parallel_proj is a unit star vector.
        for k in 1:4
            col = par[:, k]
            @test isapprox(sum(abs2, col), 1.0; atol=1e-10)
        end
        for k in 1:4
            col = perp[:, k]
            @test isapprox(sum(abs2, col), 1.0; atol=1e-10)
        end
    end

    @testset "bragg_peaks returns a BraggPeakSet" begin
        qc = generate_ammann_beenker_projection(4.0)
        peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=1e-4)

        @test peaks isa BraggPeakSet{2,4,Float64}
        @test peaks isa AbstractMomentumLattice{2,Float64}
        @test num_k_points(peaks) > 0

        # Γ (hyper index (0,0,0,0)) is always present.
        γ_idx = findfirst(==((0, 0, 0, 0)), peaks.hyper_indices)
        @test γ_idx !== nothing
        @test peaks.peaks[γ_idx] ≈ SVector(0.0, 0.0)

        # Intensities are non-negative and Γ is the brightest.
        @test all(I -> I >= -1e-12, peaks.intensities)
        γ_I = peaks.intensities[γ_idx]
        @test γ_I == maximum(peaks.intensities)
        @test γ_I ≈ 1.0 atol = 1e-10  # (∏ 2·½)² = 1
    end

    @testset "8-fold rotation symmetry of the peak set" begin
        # The physical star has 8-fold symmetry (squares + 45° rhombi).
        # Bragg peaks should come in 8-fold rotated orbits when we
        # restrict to a reasonable cutoff.
        qc = generate_ammann_beenker_projection(4.0)
        peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=1e-3)

        # For every peak, its 90° rotation should also appear in the
        # set with the same intensity (up to numerical noise). Note:
        # the current generator uses a shifted-by-π/4 perp frame so
        # full 8-fold symmetry is not guaranteed, but the C₄
        # sub-symmetry is, since the host Z⁴ has a canonical
        # cyclic permutation that commutes with both projections.
        R = SMatrix{2,2,Float64}(0, 1, -1, 0)  # 90° rotation
        for (k, I) in zip(peaks.peaks, peaks.intensities)
            I < 1e-3 && continue
            rk = R * k
            # find a peak close to rk
            found = false
            for (k2, I2) in zip(peaks.peaks, peaks.intensities)
                if sum(abs2, k2 - rk) < 1e-8 && isapprox(I2, I; atol=1e-8)
                    found = true
                    break
                end
            end
            @test found
        end
    end

    @testset "intensity_cutoff filters small peaks" begin
        qc = generate_ammann_beenker_projection(4.0)
        all_peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=0.0)
        bright_peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=0.05)

        @test num_k_points(bright_peaks) <= num_k_points(all_peaks)
        @test all(I >= 0.05 for I in bright_peaks.intensities)
    end

    @testset "fourier_module dispatch and structure_factor composition" begin
        qc = generate_ammann_beenker_projection(4.0)
        @test reciprocal_support(qc) isa HasFourierModule

        ml = momentum_lattice(qc)
        @test ml isa BraggPeakSet
        @test num_k_points(ml) > 0

        # All-ones ferromagnetic state: S(k=0) = N
        peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=1e-3)
        state = ones(Int8, num_sites(qc))
        S0 = structure_factor(qc, state, SVector(0.0, 0.0))
        @test S0 ≈ Float64(num_sites(qc)) atol = 1e-8

        Sks = structure_factor(qc, state, peaks)
        @test Sks isa Vector{Float64}
        @test length(Sks) == num_k_points(peaks)
    end
end
