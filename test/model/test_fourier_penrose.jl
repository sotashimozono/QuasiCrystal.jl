@testset "Penrose P3 Fourier analysis" begin
    @testset "BoxWindow{3} Fourier transform" begin
        w = BoxWindow(SVector(0.5, 0.5, 0.5))

        # q = 0 → volume = ∏ 2·aᵢ = 1.0
        @test window_fourier(w, SVector(0.0, 0.0, 0.0)) ≈ 1.0

        # Separability on a non-trivial point.
        q = SVector(0.7, 1.3, 2.1)
        expected = prod((2 * sin(qi * 0.5) / qi) for qi in q)
        @test isapprox(window_fourier(w, q), expected; atol=1e-12)

        # Even symmetry under q → -q.
        @test window_fourier(w, SVector(1.1, 2.7, 0.9)) ≈
            window_fourier(w, SVector(-1.1, -2.7, -0.9))
    end

    @testset "hyper_reciprocal_lattice for Penrose P3" begin
        qc = generate_penrose_projection(4.0)
        hrl = hyper_reciprocal_lattice(qc)

        @test hrl isa HyperReciprocalLattice{2,5,3,Float64,<:BoxWindow}
        @test hrl.hyper_basis ≈ SMatrix{5,5,Float64}(2π * I)

        par = hrl.parallel_proj
        perp = hrl.perp_proj
        @test size(par) == (2, 5)
        @test size(perp) == (3, 5)

        # Each parallel_proj column is a unit star vector.
        for k in 1:5
            @test isapprox(sum(abs2, par[:, k]), 1.0; atol=1e-10)
        end

        # Sum of the 5 physical star vectors is zero (5-fold star).
        s = par[:, 1] + par[:, 2] + par[:, 3] + par[:, 4] + par[:, 5]
        @test isapprox(sum(abs2, s), 0.0; atol=1e-10)
    end

    @testset "bragg_peaks returns a BraggPeakSet" begin
        qc = generate_penrose_projection(4.0)
        peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=1e-3)

        @test peaks isa BraggPeakSet{2,5,Float64}
        @test peaks isa AbstractMomentumLattice{2,Float64}
        @test num_k_points(peaks) > 0

        γ_idx = findfirst(==((0, 0, 0, 0, 0)), peaks.hyper_indices)
        @test γ_idx !== nothing
        @test peaks.peaks[γ_idx] ≈ SVector(0.0, 0.0)

        @test all(I -> I >= -1e-12, peaks.intensities)
        γ_I = peaks.intensities[γ_idx]
        @test γ_I == maximum(peaks.intensities)
        @test γ_I ≈ 1.0 atol = 1e-10  # (∏ 2·½)² = 1
    end

    @testset "5-fold rotation symmetry of the peak set" begin
        # The physical star has C₅ symmetry by construction (the
        # 5 columns of parallel_proj are a 72°-rotation orbit), so
        # the Bragg peak set must be closed under 72° rotation.
        qc = generate_penrose_projection(4.0)
        peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=1e-2)

        c, s = cos(2π / 5), sin(2π / 5)
        R = SMatrix{2,2,Float64}(c, s, -s, c)  # 72° rotation

        for (k, I) in zip(peaks.peaks, peaks.intensities)
            I < 1e-2 && continue
            rk = R * k
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
        qc = generate_penrose_projection(4.0)
        all_peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=0.0)
        bright_peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=0.05)

        @test num_k_points(bright_peaks) <= num_k_points(all_peaks)
        @test all(I >= 0.05 for I in bright_peaks.intensities)
    end

    @testset "fourier_module dispatch and structure_factor composition" begin
        qc = generate_penrose_projection(4.0)
        @test reciprocal_support(qc) isa HasFourierModule

        ml = momentum_lattice(qc)
        @test ml isa BraggPeakSet
        @test num_k_points(ml) > 0

        peaks = bragg_peaks(qc; kmax=5.0, intensity_cutoff=1e-2)
        state = ones(Int8, num_sites(qc))
        S0 = structure_factor(qc, state, SVector(0.0, 0.0))
        @test S0 ≈ Float64(num_sites(qc)) atol = 1e-8

        Sks = structure_factor(qc, state, peaks)
        @test Sks isa Vector{Float64}
        @test length(Sks) == num_k_points(peaks)
    end
end
