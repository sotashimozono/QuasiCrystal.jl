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

function _regenerate(::FibonacciLattice, g::Int, m::SubstitutionMethod)
    return generate_fibonacci_substitution(g; method=m)
end
function _regenerate(::PenroseP3, g::Int, m::SubstitutionMethod)
    return generate_penrose_substitution(g; method=m)
end
function _regenerate(::AmmannBeenker, g::Int, m::SubstitutionMethod)
    return generate_ammann_beenker_substitution(g; method=m)
end
function _regenerate(topo, ::Int, ::SubstitutionMethod)
    return throw(
        ArgumentError(
            "rescale: no substitution generator is registered for $(typeof(topo))"
        ),
    )
end

# The generators always build with the default layout, so a non-default one has to be reattached.
# A quasicrystal layout is a property of the site type, not of the patch size, so this is sound.
function _with_layout(d::QuasicrystalData{D,T,Topo,TT}, layout) where {D,T,Topo,TT}
    layout === d.layout && return d
    return QuasicrystalData{D,T,Topo,TT,typeof(layout)}(
        d.topology,
        d.positions,
        d.tiles,
        d.generation_method,
        d.parameters,
        d.bonds,
        d.nearest_neighbors,
        layout,
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
    method isa SubstitutionMethod || throw(
        ArgumentError(
            "rescale is defined only for substitution-generated quasicrystals; this patch was built " *
            "by $(typeof(method)), whose extent is set by a window radius rather than a rule depth " *
            "(scaling_rule = $(LatticeCore.scaling_rule(d))).",
        ),
    )
    haskey(d.parameters, :generations) || throw(
        ArgumentError(
            "rescale: this patch records no `:generations`, so its substitution depth is unknown.",
        ),
    )
    g = d.parameters[:generations]::Int + k * LatticeCore.scaling_rule(d).depth_step
    g >= 0 || throw(
        ArgumentError(
            "cannot step down $(-k) generation(s) from $(d.parameters[:generations]): " *
            "the generation count would be $g.",
        ),
    )
    return _with_layout(_regenerate(d.topology, g, method), d.layout)
end

# The Fibonacci substitution word at generation `g` (axiom `L`, rule
# `L → LS`, `S → L`) — the same sequence `generate_fibonacci_substitution`
# builds, so vertex indices line up: a `g`-generation chain has this many
# letters (intervals) and one more vertex.
function _fibonacci_word(g::Int)
    seq = ['L']
    for _ in 1:g
        new = Char[]
        for s in seq
            s == 'L' ? append!(new, ('L', 'S')) : push!(new, 'L')
        end
        seq = new
    end
    return seq
end

# Apply `k` substitutions to the generation-`(g-k)` word while carrying,
# for every resulting interval, the index of the coarse (generation-`g-k`)
# letter it descends from. Returns that parentage label per fine interval.
function _fibonacci_parentage(g::Int, k::Int)
    word = _fibonacci_word(g - k)
    labels = collect(1:length(word))            # one coarse cell per coarse letter
    for _ in 1:k
        nw, nl = Char[], Int[]
        for (ch, lab) in zip(word, labels)
            if ch == 'L'
                append!(nw, ('L', 'S'))
                append!(nl, (lab, lab))
            else
                push!(nw, 'L')
                push!(nl, lab)
            end
        end
        word, labels = nw, nl
    end
    return labels                               # length = #fine intervals
end

"""
    cell_partition(d::QuasicrystalData, k::Integer = 1) -> Vector{Vector{Int}}

Group the vertices of `d` by the coarse cell — the inflated tile of
`rescale(d, -k)` — they descend from, i.e. the inflation parentage.

Implemented for **Fibonacci substitution** chains, where the parentage
is exact: the generation-`g` word is obtained from the generation-`(g-k)`
word by `k` applications of `L → LS`, `S → L`, so each coarse tile owns a
contiguous run of fine intervals. A fine vertex is assigned to the coarse
cell of the interval it starts; the final vertex joins the last cell.
There are as many cells as `rescale(d, -k)` has tiles.

Raises for projection-generated patches and the other substitution
families (Penrose / Ammann–Beenker), whose inflation may re-centre so a
containment-based parentage is not well defined.
"""
function LatticeCore.cell_partition(d::QuasicrystalData, k::Integer=1)
    if !(d.topology isa FibonacciLattice && d.generation_method isa SubstitutionMethod)
        throw(
            ArgumentError(
                "cell_partition is implemented only for Fibonacci substitution chains; " *
                "this patch is $(typeof(d.topology)) / $(typeof(d.generation_method)), " *
                "whose inflation parentage is not tracked (scaling_rule = " *
                "$(LatticeCore.scaling_rule(d))).",
            ),
        )
    end
    k >= 1 || throw(ArgumentError("cell_partition needs k >= 1, got k = $k"))
    haskey(d.parameters, :generations) ||
        throw(ArgumentError("cell_partition: this patch records no `:generations`."))
    g = d.parameters[:generations]::Int
    g - k >= 0 || throw(
        ArgumentError(
            "cannot deflate k=$k from generation $g: the coarse generation would be $(g - k).",
        ),
    )

    labels = _fibonacci_parentage(g, k)
    ncells = maximum(labels)
    cells = [Int[] for _ in 1:ncells]
    for (interval, lab) in enumerate(labels)
        push!(cells[lab], interval)             # left-endpoint vertex of the interval
    end
    push!(cells[labels[end]], length(labels) + 1)   # rightmost vertex → last cell
    return cells
end

# ---- 2D inflation parentage (Penrose + Ammann–Beenker) --------------
#
# In 2D a *vertex* cell_partition is not canonical: vertices are shared
# corners of several tiles. The clean, exact object is a *tile*
# parentage — which coarse cell each fine tile descends from — defined by
# fine-tile-centre containment in the k-inflated tiling. Neither family
# is a stone inflation (fine tiles straddle coarse *tile* boundaries),
# but a covering of coarse cells whose interiors partition the plane
# fixes each centre uniquely: coarse Robinson triangles for Penrose,
# and the coarse tiles themselves for Ammann–Beenker (its square/rhombus
# tiling already covers the plane).

# Point-in-triangle by consistent edge orientation (closed test).
@inline function _in_triangle(p, A, B, C)
    d1 = (p[1] - B[1]) * (A[2] - B[2]) - (A[1] - B[1]) * (p[2] - B[2])
    d2 = (p[1] - C[1]) * (B[2] - C[2]) - (B[1] - C[1]) * (p[2] - C[2])
    d3 = (p[1] - A[1]) * (C[2] - A[2]) - (C[1] - A[1]) * (p[2] - A[2])
    has_neg = (d1 < -1e-9) | (d2 < -1e-9) | (d3 < -1e-9)
    has_pos = (d1 > 1e-9) | (d2 > 1e-9) | (d3 > 1e-9)
    return !(has_neg & has_pos)
end

# Point-in-convex-polygon (vertices in order): `p` is inside if it lies
# on the same side of every directed edge.
@inline function _in_convex(p, verts)
    n = length(verts)
    sgn = 0
    for i in 1:n
        a = verts[i]
        b = verts[mod1(i + 1, n)]
        cr = (b[1] - a[1]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[1] - a[1])
        s = cr > 1e-9 ? 1 : (cr < -1e-9 ? -1 : 0)
        s == 0 && continue
        if sgn == 0
            sgn = s
        elseif s != sgn
            return false
        end
    end
    return true
end

# The coarse Robinson triangles at generation `g - k`, expressed in the
# same edge-1 frame as `generate_penrose_substitution(g)` (same seed and
# deflation, scaled by the generation-`g` rescale factor), so they nest
# the fine tiling exactly.
function _penrose_coarse_triangles(g::Int, k::Int)
    fine = sun_seed_blue()
    for _ in 1:g
        fine = deflate_robinson(fine)
    end
    s = isempty(fine) ? 1.0 : 1.0 / norm(fine[1].A - fine[1].B)
    coarse = sun_seed_blue()
    for _ in 1:(g - k)
        coarse = deflate_robinson(coarse)
    end
    return [(t.A * s, t.B * s, t.C * s) for t in coarse]
end

"""
    tile_parentage(d::QuasicrystalData, k::Integer = 1) -> Vector{Vector{Int}}

Group the tile indices of `d` by the coarse cell — the inflated
Robinson triangle of the generation-`(k)`-deflated tiling — each fine
tile descends from. Entry `c` lists the indices (into `d.tiles`) of the
fine tiles whose centre lies in coarse cell `c`; empty coarse cells (no
fine tile in this patch) are dropped.

The **tile** parentage is the well-posed 2D analogue of the 1D
[`cell_partition`](@ref): a *vertex* partition is not canonical because
vertices are shared corners, whereas each fine-tile centre lies in
exactly one coarse Robinson triangle. Assignment is by centre
containment, which is exact because the Robinson triangles tile the
plane (a fine rhombus may straddle a coarse *rhombus* boundary, but its
centre still falls in one coarse *triangle*).

Implemented for **Penrose substitution** patches. Raises for projection
patches, the other substitution families (whose generators are not
exactly self-similar), and out-of-range `k`.
"""
function tile_parentage(d::QuasicrystalData, k::Integer=1)
    k >= 1 || throw(ArgumentError("tile_parentage needs k >= 1, got k = $k"))
    if d.topology isa PenroseP3 && d.generation_method isa SubstitutionMethod
        return _penrose_tile_parentage(d, k)
    elseif d.topology isa AmmannBeenker && !isempty(d.tiles)
        return _ab_tile_parentage(d, k)
    else
        throw(
            ArgumentError(
                "tile_parentage is implemented for Penrose substitution patches and " *
                "Ammann–Beenker patches carrying a tiling; this patch is " *
                "$(typeof(d.topology)) / $(typeof(d.generation_method))" *
                (isempty(d.tiles) ? " with no tiles" : "") *
                ".",
            ),
        )
    end
end

function _penrose_tile_parentage(d::QuasicrystalData, k::Int)
    haskey(d.parameters, :generations) ||
        throw(ArgumentError("tile_parentage: this patch records no `:generations`."))
    g = d.parameters[:generations]::Int
    g - k >= 0 || throw(
        ArgumentError(
            "cannot deflate k=$k from generation $g: the coarse generation would be $(g - k).",
        ),
    )

    tris = _penrose_coarse_triangles(g, k)
    groups = [Int[] for _ in 1:length(tris)]
    for (ti, tile) in enumerate(d.tiles)
        c = 0
        for (ci, (A, B, C)) in enumerate(tris)
            if _in_triangle(tile.center, A, B, C)
                c = ci
                break
            end
        end
        c == 0 && error(
            "tile_parentage: fine tile $ti centre fell outside every coarse Robinson " *
            "triangle — the coarse frame does not match (unexpected).",
        )
        push!(groups[c], ti)
    end
    return filter(!isempty, groups)
end

# Ammann–Beenker is not a stone inflation — fine tiles straddle coarse
# *tile* boundaries — but every fine-tile centre lies in exactly one
# coarse tile (the coarse tiling covers the plane, unlike Penrose's
# coarse rhombi), so centre containment gives an exact tile parentage.
# The coarse tiling is the same AB tiling `k` inflations up: regenerate
# it at radius `R/λᵏ` and scale by `λᵏ` into the fine frame — exact
# because the AB point set is self-similar (`λ·S ⊆ S`, see
# `inflation_factor`).
function _ab_tile_parentage(d::QuasicrystalData, k::Int)
    haskey(d.parameters, :radius) || throw(
        ArgumentError(
            "tile_parentage: only projection-generated Ammann–Beenker tilings are " *
            "supported (this patch records no `:radius`; the AB substitution generator " *
            "is a placeholder).",
        ),
    )
    R = Float64(d.parameters[:radius])
    scale = (1 + sqrt(2))^k
    # Over-cover `d`: the coarse patch, scaled by `scale`, must reach past
    # `d`'s outermost fine tiles, so pad the coarse radius by a few coarse
    # tiles. Coarse tiles outside `d` simply carry no fine tiles.
    coarse = generate_ammann_beenker_projection(R / scale + 3.0)
    isempty(coarse.tiles) && throw(
        ArgumentError("cannot deflate k=$k: the coarse Ammann–Beenker tiling is empty.")
    )
    coarse_verts = [[scale .* v for v in t.vertices] for t in coarse.tiles]

    groups = [Int[] for _ in 1:length(coarse_verts)]
    for (ti, tile) in enumerate(d.tiles)
        c = 0
        for (ci, verts) in enumerate(coarse_verts)
            if _in_convex(tile.center, verts)
                c = ci
                break
            end
        end
        c == 0 && error(
            "tile_parentage: fine tile $ti centre fell outside every coarse " *
            "Ammann–Beenker tile (unexpected).",
        )
        push!(groups[c], ti)
    end
    return filter(!isempty, groups)
end
