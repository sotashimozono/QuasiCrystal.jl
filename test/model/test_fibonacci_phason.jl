using QuasiCrystal
using Test
using StaticArrays
using LinearAlgebra: eigvals

@testset "fibonacci_phason" begin
    lat = FibonacciLattice()
    α = PHASON_INTERCEPT_FIBONACCI
    @test α ≈ inv(ϕ)
    @test 0 < α < 1

    @testset "fibonacci_word" begin
        # axiom is [L] ≡ [0]
        @test fibonacci_word(lat, 0) == [0]
        # gen 1: L → LS  ⇒  [0, 1]
        @test fibonacci_word(lat, 1) == [0, 1]
        # gen 2: LS → LS L  ⇒  [0, 1, 0]
        @test fibonacci_word(lat, 2) == [0, 1, 0]
        # gen 3: LS L → LS L LS  ⇒  [0, 1, 0, 0, 1]
        @test fibonacci_word(lat, 3) == [0, 1, 0, 0, 1]

        # length = F_{gen+2}  (F_1 = F_2 = 1)
        @test length(fibonacci_word(lat, 0)) == fibonacci_sequence_length(0)
        @test length(fibonacci_word(lat, 5)) == fibonacci_sequence_length(5)
        @test length(fibonacci_word(lat, 8)) == fibonacci_sequence_length(8)

        # alphabet is {0, 1} only
        w = fibonacci_word(lat, 7)
        @test all(s -> s == 0 || s == 1, w)

        # eventually-L ratio → α (equivalently #L / #S → ϕ)
        big = fibonacci_word(lat, 14)
        n_L = count(==(0), big)
        @test n_L / length(big) ≈ α atol = 1e-3
    end

    @testset "phason_orbit_at" begin
        @test phason_orbit_at(lat, 0) == 0.0
        @test phason_orbit_at(lat, 1) ≈ α
        # closed-form for θ0 ≠ 0
        @test phason_orbit_at(lat, 3; θ0=0.25) ≈ mod(0.25 + 3α, 1.0)
        # custom slope override
        @test phason_orbit_at(lat, 4; α=0.2) ≈ 0.8 atol = 1e-12
        # range [0, 1)
        for i in -3:7
            θ = phason_orbit_at(lat, i)
            @test 0 ≤ θ < 1
        end
        # equidistribution (Weyl) on the first ~10k iterates
        N = 10_000
        orbit = [phason_orbit_at(lat, i) for i in 0:(N - 1)]
        # Halfway threshold should split close to 50/50.
        @test abs(count(x -> x < 0.5, orbit) / N - 0.5) < 0.01
    end

    @testset "bond_couplings" begin
        # explicit gen 4 word = [0,1,0,0,1,0,1,0]; 7 bonds
        bs = bond_couplings(lat, 1.0, 0.5; gen=4)
        @test bs == [1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 0.5]

        # length = length(word) - 1
        for gen in (0, 1, 2, 6)
            w = fibonacci_word(lat, gen)
            b = bond_couplings(lat, 1.0, 2.0; gen=gen)
            @test length(b) == length(w) - 1
        end

        # mean over a long chain → α·JL + (1-α)·JS
        JL, JS = 1.0, 0.5
        bs_long = bond_couplings(lat, JL, JS; gen=14)
        @test sum(bs_long) / length(bs_long) ≈ α * JL + (1 - α) * JS atol = 1e-3
    end

    @testset "phason_grid_shift" begin
        # σ is a permutation of 1:M for every M
        for M in (4, 8, 13, 21, 50)
            σ = phason_grid_shift(lat, M)
            @test length(σ) == M
            @test sort(σ) == collect(1:M)
            # |θ_{σ(k)} - (θ_k + α)| < 1/M (mod 1)
            for k in 1:M
                θk = (k - 1) / M
                θs = (σ[k] - 1) / M
                diff = mod(θs - θk - α + 0.5, 1.0) - 0.5
                @test abs(diff) ≤ 1 / M
            end
        end

        # custom α override (rational slope ⇒ exact shift)
        @test phason_grid_shift(lat, 5; α=0.2) == [2, 3, 4, 5, 1]
    end

    @testset "inflation_matrix" begin
        M = inflation_matrix(lat)
        @test M == SMatrix{2,2,Int}(1, 1, 1, 0)
        # eigenvalues = (ϕ, -1/ϕ)
        evs = eigvals(Float64.(M))
        @test maximum(evs) ≈ ϕ atol = 1e-12
        @test minimum(evs) ≈ -inv(ϕ) atol = 1e-12
    end

    @testset "J_fourier_coeffs" begin
        JL, JS = 1.0, 0.5
        K_J = 20
        Ĵ = J_fourier_coeffs(lat, JL, JS, K_J)
        @test length(Ĵ) == 2K_J + 1

        # DC term = orbit average
        Ĵ_0 = Ĵ[K_J + 1]
        @test real(Ĵ_0) ≈ α * JL + (1 - α) * JS atol = 1e-12
        @test abs(imag(Ĵ_0)) < 1e-12

        # J(θ) is real ⇒ Ĵ_{-k} = conj(Ĵ_k)
        for k in 1:K_J
            @test Ĵ[K_J + 1 - k] ≈ conj(Ĵ[K_J + 1 + k]) atol = 1e-12
        end

        # reconstruct J(θ) on a dense grid via partial sum
        function J_truth(θ::Real, α::Real)
            θm = mod(θ, 1.0)
            if 0.0 ≤ θm < 1 - α
                return JL
            elseif 1 - α ≤ θm < 2 - 2α
                return JS
            else
                return JL
            end
        end

        function J_recon(θ::Real, Ĵ::Vector{ComplexF64}, K_J::Int)
            s = 0.0 + 0.0im
            for (idx, k) in enumerate((-K_J):K_J)
                s += Ĵ[idx] * cis(2π * k * θ)
            end
            return real(s)
        end

        # Sample θ away from step discontinuities (Gibbs-free zone).
        θs = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
        for θ in θs
            jt = J_truth(θ, α)
            jr = J_recon(θ, Ĵ, K_J)
            @test abs(jt - jr) < 0.05
        end

        # Custom α path stays consistent at DC.
        Ĵ_alt = J_fourier_coeffs(lat, 1.0, 0.0, 8; α=0.3)
        @test real(Ĵ_alt[8 + 1]) ≈ 0.3 atol = 1e-12
    end
end
