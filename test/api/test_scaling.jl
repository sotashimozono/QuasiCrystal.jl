using QuasiCrystal
using LatticeCore
using Test

# The three substitution families, with the generation counts their builders take.
const SUBST = (
    ("fibonacci", fibonacci, 6),
    ("penrose", penrose, 2),
    ("ammann_beenker", ammann_beenker, 1),
)

@testset "scaling_rule follows the generation method" begin
    for (name, f, n) in SUBST
        d = f(n)
        @test scaling_rule(d) === SubstitutionScaling(1)
    end
    # a projection-generated patch has no rule depth to advance
    proj = generate_fibonacci_projection(30)
    @test scaling_rule(proj) === NoScaling()
end

@testset "rescale advances the generation" begin
    for (name, f, n) in SUBST
        d = f(n)
        @test rescale(d, 0) === d
        # one step up matches building at n+1 directly
        @test num_sites(rescale(d)) == num_sites(f(n + 1))
        @test num_sites(rescale(d, 2)) == num_sites(f(n + 2))
        # and strictly grows
        @test num_sites(rescale(d)) > num_sites(d)
        # stepping back down recovers the smaller patch
        @test num_sites(rescale(rescale(d), -1)) == num_sites(d)
        # topology and generation method survive
        up = rescale(d)
        @test typeof(up.topology) === typeof(d.topology)
        @test typeof(up.generation_method) === typeof(d.generation_method)
        @test up.parameters[:generations] == d.parameters[:generations] + 1
    end
end

@testset "rescale refuses what it cannot do" begin
    d = fibonacci(2)
    # generation count may not go negative
    @test_throws ArgumentError rescale(d, -3)
    # a projection-generated patch has no substitution depth
    @test_throws ArgumentError rescale(generate_fibonacci_projection(30))
end

@testset "size_sequence over the substitution families" begin
    for (name, f, n) in SUBST
        d = f(n)
        seq = size_sequence(d, 2)
        @test length(seq) == 3
        counts = num_sites.(seq)
        @test counts == [num_sites(f(n)), num_sites(f(n + 1)), num_sites(f(n + 2))]
        @test issorted(counts)                      # each step is at least as large
    end
end

@testset "the Fibonacci chain grows by the Fibonacci numbers" begin
    # An independent check that rescale advances the substitution by exactly one generation, and not
    # merely that the patch gets bigger: `fibonacci_sequence_length` computes the word length from
    # the recurrence, without generating anything.
    #
    # The word length counts LETTERS, i.e. intervals; `num_sites` counts VERTICES, i.e. endpoints.
    # A chain of m intervals has m+1 endpoints, hence the offset.
    d = fibonacci(4)
    for (k, l) in enumerate(size_sequence(d, 4))
        gen = 4 + k - 1
        @test num_sites(l) == fibonacci_sequence_length(gen) + 1
    end
    # and the letter counts really are consecutive Fibonacci numbers
    lens = [fibonacci_sequence_length(g) for g in 4:8]
    @test all(lens[i + 2] == lens[i + 1] + lens[i] for i in 1:(length(lens) - 2))
end

@testset "the layout is carried over" begin
    d = fibonacci(5)
    layout = UniformLayout(PottsSite(3))
    relaid = QuasiCrystal._with_layout(d, layout)
    @test relaid.layout === layout
    up = rescale(relaid)
    @test up.layout === layout                      # not reset to the generator's default
    @test num_sites(up) == num_sites(fibonacci(6))
end

@testset "cell_partition is refused, with a reason" begin
    # Inflation parentage is not tracked, so this must raise rather than return something plausible.
    for (name, f, n) in SUBST
        @test_throws ArgumentError cell_partition(f(n))
        @test_throws ArgumentError cell_partition(f(n), 2)
    end
end

@testset "rescale returns a fresh patch, with bonds unbuilt" begin
    # `rescale` calls the generator, so bonds start empty exactly as after any other generator call.
    # Documented behaviour; pinned here so it cannot change silently.
    d = fibonacci(5)
    @test isempty(rescale(d).bonds)
end
