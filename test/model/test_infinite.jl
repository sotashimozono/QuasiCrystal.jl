using QuasiCrystal
using Test

const _QC_FAMILIES = (FibonacciLattice(), PenroseP3(), AmmannBeenker())

@testset "InfiniteQuasicrystal — traits and undefined site API" begin
    for topo in _QC_FAMILIES
        inf = InfiniteQuasicrystal(topo)

        @test size_trait(inf) == InfiniteSize()
        @test !is_finite(inf)
        @test periodicity(inf) == Aperiodic()
        @test reciprocal_support(inf) == HasFourierModule()
        @test topology(inf) == TopologyTrait{:quasiperiodic}()

        # No sites exist until materialised.
        @test_throws DomainError num_sites(inf)
        @test_throws DomainError position(inf, 1)
        @test_throws DomainError neighbors(inf, 1)
        # The Monte Carlo entry guard must reject the infinite object.
        @test_throws ArgumentError require_finite(inf)
    end
end

@testset "materialize matches the direct projection generator" begin
    # Independent oracle: the lazy materialize path must reproduce, site
    # for site, what the underlying generator produces directly.
    fib = materialize(InfiniteQuasicrystal(FibonacciLattice()); n_points=50)
    @test fib isa QuasicrystalData
    @test is_finite(fib)
    @test num_sites(fib) == 50
    @test positions(fib) == positions(generate_fibonacci_projection(50))
    @test require_finite(fib) === nothing   # finite ⇒ no throw

    pen = materialize(InfiniteQuasicrystal(PenroseP3()); radius=4.0)
    @test pen isa QuasicrystalData
    @test num_sites(pen) == num_sites(generate_penrose_projection(4.0))
    @test positions(pen) == positions(generate_penrose_projection(4.0))

    ab = materialize(InfiniteQuasicrystal(AmmannBeenker()); radius=4.0)
    @test ab isa QuasicrystalData
    @test num_sites(ab) == num_sites(generate_ammann_beenker_projection(4.0))
    @test positions(ab) == positions(generate_ammann_beenker_projection(4.0))
end

@testset "materialize propagates the site layout" begin
    inf = InfiniteQuasicrystal(FibonacciLattice(); layout=UniformLayout(XYSite()))
    # UniformLayout ignores the index, so site_type resolves even on the
    # infinite object.
    @test site_type(inf, 1) == XYSite()

    d = materialize(inf; n_points=20)
    @test site_layout(d) === site_layout(inf)   # same layout object carried through
    @test site_type(d, 1) == XYSite()
    @test num_sites(d) == 20
end

@testset "materialize default cutoffs produce non-empty lattices" begin
    for topo in _QC_FAMILIES
        d = materialize(InfiniteQuasicrystal(topo))
        @test d isa QuasicrystalData
        @test is_finite(d)
        @test num_sites(d) > 0
    end
end

# ---- RG inflate / deflate -------------------------------------------

const _PHI = (1 + sqrt(5)) / 2

# every point of `sub` lies (within tol) on a point of `base`?
function _is_subset(sub, base; tol=1e-8)
    return all(y -> minimum(sum(abs2, y - x) for x in base) <= tol^2, sub)
end

@testset "inflation_factor and lazy inflate/deflate bookkeeping" begin
    @test inflation_factor(FibonacciLattice()) ≈ _PHI
    @test inflation_factor(PenroseP3()) ≈ _PHI

    inf = InfiniteQuasicrystal(FibonacciLattice())
    @test inf.inflation == 0
    @test inflate(inf).inflation == 1
    @test inflate(inf, 3).inflation == 3
    @test deflate(inf).inflation == -1
    @test deflate(inflate(inf)).inflation == 0          # inverse
    # inflate/deflate carry topology and layout over.
    inf_xy = InfiniteQuasicrystal(FibonacciLattice(); layout=UniformLayout(XYSite()))
    @test site_type(inflate(inf_xy), 1) == XYSite()
end

@testset "inflate scales materialised positions by λ (exact self-map)" begin
    # Fibonacci (count cutoff): materialise(inflate) == ϕ · materialise(base).
    base = positions(materialize(InfiniteQuasicrystal(FibonacciLattice()); n_points=60))
    up = positions(
        materialize(inflate(InfiniteQuasicrystal(FibonacciLattice())); n_points=60)
    )
    @test length(up) == length(base)
    @test all(up[i] ≈ _PHI .* base[i] for i in eachindex(base))

    # deflate is the exact inverse.
    back = positions(
        materialize(deflate(inflate(InfiniteQuasicrystal(FibonacciLattice()))); n_points=60)
    )
    @test all(back[i] ≈ base[i] for i in eachindex(base))

    # Penrose (radius cutoff): same elementwise λ scaling.
    pbase = positions(materialize(InfiniteQuasicrystal(PenroseP3()); radius=6.0))
    pup = positions(materialize(inflate(InfiniteQuasicrystal(PenroseP3())); radius=6.0))
    @test length(pup) == length(pbase)
    @test all(pup[i] ≈ _PHI .* pbase[i] for i in eachindex(pbase))
end

@testset "self-similarity: inflated points are genuine points of the set" begin
    # The deep RG property (not just coordinate rescaling): the inflated
    # (sparser) set is a SUBSET of the denser base set — λ·S ⊆ S.
    # Fibonacci: base must cover the ϕ-scaled range, so generate more.
    fib = FibonacciLattice()
    infl = positions(materialize(inflate(InfiniteQuasicrystal(fib)); n_points=80))
    dense = positions(materialize(InfiniteQuasicrystal(fib); n_points=200))
    @test _is_subset(infl, dense)

    # Penrose: base radius must reach ϕ·R to contain the inflated points.
    pen = PenroseP3()
    p_infl = positions(materialize(inflate(InfiniteQuasicrystal(pen)); radius=5.0))
    p_dense = positions(materialize(InfiniteQuasicrystal(pen); radius=5.0 * _PHI + 1.0))
    @test _is_subset(p_infl, p_dense)
end

@testset "Ammann–Beenker inflate/deflate (now self-similar)" begin
    ab = InfiniteQuasicrystal(AmmannBeenker())
    λ = 1 + sqrt(2)
    @test inflation_factor(AmmannBeenker()) ≈ λ
    @test inflate(ab).inflation == 1
    @test deflate(inflate(ab)).inflation == 0

    # inflate scales the materialised positions by exactly λ.
    base = positions(materialize(ab; radius=8.0))
    up = positions(materialize(inflate(ab); radius=8.0))
    @test length(up) == length(base)
    @test all(up[i] ≈ λ .* base[i] for i in eachindex(base))

    # self-similarity: the inflated (sparser) set is a subset of a denser
    # base set — λ·S ⊆ S, the fix that made AB well-posed.
    infl = positions(materialize(inflate(ab); radius=6.0))
    dense = positions(materialize(ab; radius=6.0 * λ + 1.0))
    @test _is_subset(infl, dense)
end
