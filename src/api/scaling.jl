# Scale changes for quasicrystals: the LatticeCore scaling interface.
#
# A substitution-generated patch has a natural next size — one more generation of the rule — so it
# takes part in `size_sequence` even though it has no side length to double. A projection-generated
# patch does not: its extent is set by a window radius, not by a rule depth, so it reports
# `NoScaling()` and stays out.

"""
    scaling_rule(d::QuasicrystalData) -> AbstractScalingRule

[`SubstitutionScaling`](@ref)`(1)` for a patch built by a substitution rule — one scale step is one
more generation. [`NoScaling`](@ref) for a projection-generated patch, whose extent comes from a
window radius rather than a rule depth and therefore has no canonical "next size".
"""
function LatticeCore.scaling_rule(d::QuasicrystalData)
    return d.generation_method isa SubstitutionMethod ? SubstitutionScaling(1) : NoScaling()
end

_regenerate(::FibonacciLattice, g::Int, m::SubstitutionMethod) =
    generate_fibonacci_substitution(g; method=m)
_regenerate(::PenroseP3, g::Int, m::SubstitutionMethod) =
    generate_penrose_substitution(g; method=m)
_regenerate(::AmmannBeenker, g::Int, m::SubstitutionMethod) =
    generate_ammann_beenker_substitution(g; method=m)
_regenerate(topo, ::Int, ::SubstitutionMethod) = throw(ArgumentError(
    "rescale: no substitution generator is registered for $(typeof(topo))"
))

# The generators always build with the default layout, so a non-default one has to be reattached.
# A quasicrystal layout is a property of the site type, not of the patch size, so this is sound.
function _with_layout(d::QuasicrystalData{D,T,Topo,TT}, layout) where {D,T,Topo,TT}
    layout === d.layout && return d
    return QuasicrystalData{D,T,Topo,TT,typeof(layout)}(
        d.topology, d.positions, d.tiles, d.generation_method, d.parameters,
        d.bonds, d.nearest_neighbors, layout,
    )
end

"""
    rescale(d::QuasicrystalData, k::Integer = 1) -> QuasicrystalData

The same construction `k` substitution generations further on, i.e. regenerated at
`parameters[:generations] + k`. Topology, generation method and site layout are carried over.

Only defined for substitution-generated patches; a projection-generated one raises, as does a
negative `k` that would take the generation count below zero.

!!! note "The result is a fresh patch"
    `rescale` calls the generator, so — exactly as after any other generator call — `bonds` and
    `nearest_neighbors` start empty and must be rebuilt with `build_nearest_neighbor_bonds!` before
    anything walks the graph.

!!! note "Successive generations need not be nested"
    For the Fibonacci chain the substitution extends the word by a prefix, so generation `n+1`
    contains generation `n`. Robinson inflation of a Penrose or Ammann–Beenker patch may re-centre
    and rescale, so containment is not guaranteed in general. `size_sequence` gives *the same
    construction at successive generations*, not a nested family.

```jldoctest
julia> using QuasiCrystal, LatticeCore

julia> d = fibonacci(6);

julia> scaling_rule(d)
SubstitutionScaling(1)

julia> num_sites.(size_sequence(d, 3))
4-element Vector{Int64}:
 14
 22
 35
 56
```
"""
function LatticeCore.rescale(d::QuasicrystalData, k::Integer=1)
    k == 0 && return d
    method = d.generation_method
    method isa SubstitutionMethod || throw(ArgumentError(
        "rescale is defined only for substitution-generated quasicrystals; this patch was built " *
        "by $(typeof(method)), whose extent is set by a window radius rather than a rule depth " *
        "(scaling_rule = $(LatticeCore.scaling_rule(d)))."
    ))
    haskey(d.parameters, :generations) || throw(ArgumentError(
        "rescale: this patch records no `:generations`, so its substitution depth is unknown."
    ))
    g = d.parameters[:generations]::Int + k * LatticeCore.scaling_rule(d).depth_step
    g >= 0 || throw(ArgumentError(
        "cannot step down $(-k) generation(s) from $(d.parameters[:generations]): " *
        "the generation count would be $g."
    ))
    return _with_layout(_regenerate(d.topology, g, method), d.layout)
end

"""
    cell_partition(d::QuasicrystalData, k::Integer = 1)

Not available for quasicrystals: grouping the vertices of a patch by the inflated tile they belong
to needs the parentage of each tile through the substitution, which the generators do not currently
record. Always raises.
"""
function LatticeCore.cell_partition(d::QuasicrystalData, k::Integer=1)
    return throw(ArgumentError(
        "cell_partition is not available for QuasicrystalData: inflation parentage (which vertex " *
        "of generation n sits in which tile of generation n−$k) is not tracked by the generators."
    ))
end
