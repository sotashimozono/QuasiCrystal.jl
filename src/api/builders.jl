"""
Convenience builder API: thin one-liner wrappers around the
`generate_*_substitution` and `generate_*_projection` family that
ship the most commonly-used defaults so callers do not have to
remember per-family argument names.

These builders are pure delegators: they construct no new state,
they only forward keyword arguments to the canonical generators.
The returned `QuasicrystalData` is identical (`===`-equal in
field-by-field semantics) to what the underlying `generate_*`
function would have produced.

# Naming

| Builder                      | Underlying generator                     |
|------------------------------|------------------------------------------|
| `fibonacci(n)`               | `generate_fibonacci_substitution(n)`     |
| `penrose(n)`                 | `generate_penrose_substitution(n)`       |
| `ammann_beenker(n)`          | `generate_ammann_beenker_substitution(n)`|
| `penrose_projected(r)`       | `generate_penrose_projection(r)`         |
| `ammann_beenker_projected(r)`| `generate_ammann_beenker_projection(r)`  |

The substitution-route builders default to a small number of
generations so a bare `penrose()` returns a usable point set
without explicit arguments. The projection-route builders take a
positional `radius` because there is no universally sensible
default that scales gracefully with use case.
"""

"""
    penrose(n::Int = 3; method = SubstitutionMethod(DefaultSubstitution()))
        → QuasicrystalData{2, Float64, PenroseP3}

Build a Penrose P3 quasicrystal via `n` substitution generations.
Equivalent to `generate_penrose_substitution(n; method=method)`.

`method` accepts either an `AbstractSubstitutionAlgorithm` (which
will be wrapped in a `SubstitutionMethod`) or a fully-formed
`SubstitutionMethod`.
"""
function penrose(
    n::Int=3;
    method::Union{SubstitutionMethod,AbstractSubstitutionAlgorithm}=SubstitutionMethod(),
)
    return generate_penrose_substitution(n; method=_as_substitution_method(method))
end

"""
    ammann_beenker(n::Int = 3; method = SubstitutionMethod(DefaultSubstitution()))
        → QuasicrystalData{2, Float64, AmmannBeenker}

Build an Ammann–Beenker quasicrystal via `n` substitution
generations. Equivalent to
`generate_ammann_beenker_substitution(n; method=method)`.
"""
function ammann_beenker(
    n::Int=3;
    method::Union{SubstitutionMethod,AbstractSubstitutionAlgorithm}=SubstitutionMethod(),
)
    return generate_ammann_beenker_substitution(
        n; method=_as_substitution_method(method)
    )
end

"""
    fibonacci(n::Int = 10) → QuasicrystalData{1, Float64, FibonacciLattice}

Build a 1D Fibonacci chain via `n` substitution generations.
Equivalent to `generate_fibonacci_substitution(n)`. The default
`n = 10` produces a chain of `fibonacci_sequence_length(10) + 1`
sites — large enough for typical inspection / plotting but still
cheap to construct.
"""
function fibonacci(
    n::Int=10;
    method::Union{SubstitutionMethod,AbstractSubstitutionAlgorithm}=SubstitutionMethod(),
)
    return generate_fibonacci_substitution(n; method=_as_substitution_method(method))
end

"""
    penrose_projected(radius::Real = 5.0; method = ProjectionMethod())
        → QuasicrystalData{2, Float64, PenroseP3}

Build a Penrose P3 quasicrystal via the cut-and-project route with
physical-space radius `radius`. Equivalent to
`generate_penrose_projection(radius; method=method)`.

The default `radius = 5.0` is large enough to expose 5-fold
rotational symmetry while staying cheap to construct.
"""
function penrose_projected(
    radius::Real=5.0; method::ProjectionMethod=ProjectionMethod()
)
    return generate_penrose_projection(radius; method=method)
end

"""
    ammann_beenker_projected(radius::Real = 5.0; method = ProjectionMethod())
        → QuasicrystalData{2, Float64, AmmannBeenker}

Build an Ammann–Beenker quasicrystal via the cut-and-project route
with physical-space radius `radius`. Equivalent to
`generate_ammann_beenker_projection(radius; method=method)`.
"""
function ammann_beenker_projected(
    radius::Real=5.0; method::ProjectionMethod=ProjectionMethod()
)
    return generate_ammann_beenker_projection(radius; method=method)
end

"""
    fibonacci_projected(n_points::Int = 100; method = ProjectionMethod())
        → QuasicrystalData{1, Float64, FibonacciLattice}

Build a 1D Fibonacci chain via the cut-and-project route, taking
the first `n_points` projected sites. Equivalent to
`generate_fibonacci_projection(n_points; method=method)`.
"""
function fibonacci_projected(
    n_points::Int=100; method::ProjectionMethod=ProjectionMethod()
)
    return generate_fibonacci_projection(n_points; method=method)
end

# ---- Internal helpers -----------------------------------------------

# Accept either a fully-formed `SubstitutionMethod` or a bare
# algorithm singleton (so callers can write `method=DefaultSubstitution()`
# and have it wrapped automatically).
_as_substitution_method(m::SubstitutionMethod) = m
_as_substitution_method(alg::AbstractSubstitutionAlgorithm) = SubstitutionMethod(alg)
