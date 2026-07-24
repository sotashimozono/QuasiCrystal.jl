"""
Ammann–Beenker tiling: 8-fold rotationally symmetric cut-and-project
quasicrystal built from squares and 45° rhombi, from a 4D hypercubic
host lattice.
"""

"""
    AmmannBeenker <: AbstractQuasicrystal{2}

Topology marker for the Ammann–Beenker tiling.
"""
struct AmmannBeenker <: AbstractQuasicrystal{2} end

const SQRT2 = sqrt(2.0)

"""
    generate_ammann_beenker_projection(radius::Real;
                                       method::ProjectionMethod = ProjectionMethod())
        → QuasicrystalData{2, Float64}

Generate an Ammann–Beenker patch by projecting `Z^4` onto the 2D
physical subspace of the cut-and-project construction, accepting lattice
points whose Galois-conjugate perpendicular projection falls inside the
regular-octagon window. The point set is exactly self-similar under
`λ = 1 + √2`.

The square/rhombus tiling is reconstructed from the vertices (see
[`build_tile_bonds!`](@ref)-compatible `tiles`), so the returned
`QuasicrystalData` carries a full tiling — unlike the other projection
generators, whose `tiles` are empty.
"""
# Reconstruct the Ammann–Beenker tiling (squares + 45° rhombi) from a
# vertex set. Both prototiles are parallelograms with unit edges, and in
# the AB tiling every unit-distance vertex pair is a tile edge (the
# rhombus/square diagonals are 2sin(π/8), √2, 2cos(π/8) — never 1), so
# the faces are recovered exactly: at each vertex, consecutive
# unit-neighbours (ordered by angle) bound a tile whose fourth corner is
# the parallelogram completion `a + b - v`. A 90° corner is a square,
# 45°/135° a rhombus. Tiles are deduplicated by their vertex set.
#
# O(N²) in the edge pass, matching `build_nearest_neighbor_bonds!`.
function _reconstruct_ab_tiles(positions::Vector{SVector{2,Float64}})
    n = length(positions)
    n == 0 && return Tile{2,Float64}[]
    edgetol = 1e-6

    nbr = [Int[] for _ in 1:n]
    for i in 1:n, j in (i + 1):n
        if abs(norm(positions[i] - positions[j]) - 1.0) < edgetol
            push!(nbr[i], j)
            push!(nbr[j], i)
        end
    end

    idx = Dict{NTuple{2,Int},Int}()
    for (i, p) in enumerate(positions)
        idx[snap_to_grid(p, SNAP_GRID_EPS)] = i
    end

    tiles = Tile{2,Float64}[]
    seen = Set{NTuple{4,Int}}()
    for v in 1:n
        ns = nbr[v]
        isempty(ns) && continue
        ord = sortperm([
            atan(positions[a][2] - positions[v][2], positions[a][1] - positions[v][1]) for
            a in ns
        ])
        ns = ns[ord]
        m = length(ns)
        for k in 1:m
            a = ns[k]
            b = ns[mod1(k + 1, m)]
            va = positions[a] - positions[v]
            vb = positions[b] - positions[v]
            ang = acosd(clamp(dot(va, vb) / (norm(va) * norm(vb)), -1.0, 1.0))
            (abs(ang - 45) < 1 || abs(ang - 90) < 1 || abs(ang - 135) < 1) || continue
            w = positions[a] + positions[b] - positions[v]
            wi = get(idx, snap_to_grid(w, SNAP_GRID_EPS), 0)
            wi == 0 && continue
            key = Tuple(sort(SVector(v, a, wi, b)))
            key in seen && continue
            push!(seen, key)
            type = abs(ang - 90) < 1 ? Square() : RhombusAB()
            center = (positions[v] + positions[wi]) / 2
            push!(
                tiles,
                Tile{2,Float64}(
                    [positions[v], positions[a], positions[wi], positions[b]], type, center
                ),
            )
        end
    end
    return tiles
end

function generate_ammann_beenker_projection(
    radius::Real; method::ProjectionMethod=ProjectionMethod()
)
    theta = π / 4

    E_par = zeros(4, 2)
    for i in 1:4
        E_par[i, 1] = cos((i - 1) * theta)
        E_par[i, 2] = sin((i - 1) * theta)
    end

    # Perpendicular space is the *Galois-conjugate* eigenplane: the star
    # vector `e_k` maps to angle `3·k·45°`, not `k·45° + 45°`. The 45°
    # rotation used previously is not the conjugate, so the resulting set
    # was not self-similar. With the conjugate projection the inflation
    # `M` on `Z⁴` acts as `×(1+√2)` in physical space and `×(1−√2)` in
    # perp space, so `(1+√2)·S ⊆ S` holds exactly.
    E_perp = zeros(4, 2)
    for i in 1:4
        E_perp[i, 1] = cos((i - 1) * 3 * theta)
        E_perp[i, 2] = sin((i - 1) * 3 * theta)
    end

    # Acceptance window: the regular octagon that is the perp projection
    # of the unit 4-cube `[-1/2, 1/2]⁴` — the intersection of an
    # axis-aligned and a 45°-rotated box, with apothem `(1+√2)/2`. A box
    # (the old `window_size = 0.5` square) is the wrong shape and breaks
    # self-similarity.
    apothem = (1 + sqrt(2)) / 2
    diag = apothem * sqrt(2)
    n_max = ceil(Int, radius * 1.5)

    # Pre-typed SVector buffer (avoids the legacy `Vector{Float64}[]`
    # round-trip through plain Julia vectors).
    positions = SVector{2,Float64}[]
    P = SMatrix{2,4,Float64}(E_par')
    Q = SMatrix{2,4,Float64}(E_perp')

    tol = 1e-9
    for n1 in (-n_max):n_max,
        n2 in (-n_max):n_max, n3 in (-n_max):n_max,
        n4 in (-n_max):n_max

        lp = SVector{4,Float64}(float(n1), float(n2), float(n3), float(n4))
        pos_par = P * lp
        norm(pos_par) > radius && continue
        w = Q * lp
        if abs(w[1]) <= apothem + tol &&
            abs(w[2]) <= apothem + tol &&
            abs(w[1] + w[2]) <= diag + tol &&
            abs(w[1] - w[2]) <= diag + tol
            push!(positions, pos_par)
        end
    end
    tiles = _reconstruct_ab_tiles(positions)

    params = Dict{Symbol,Any}(
        :radius => radius,
        :n_max => n_max,
        :window_apothem => apothem,
        :window_shape => :octagon_2d,
        :n_vertices => length(positions),
        :n_tiles => length(tiles),
        :symmetry => 8,
    )
    return QuasicrystalData{2,Float64}(AmmannBeenker(), positions, tiles, method, params)
end

"""
    generate_ammann_beenker_substitution(generations::Int;
                                         method::SubstitutionMethod = SubstitutionMethod())
        → QuasicrystalData{2, Float64}

Generate an Ammann–Beenker point set by iterating the substitution
(inflation) rules with an inflation factor of `1 + √2`.

The current inflation routine is a **placeholder** that scales
tiles rather than subdividing them properly into smaller squares
and rhombi. A full implementation is a future task.
"""
function generate_ammann_beenker_substitution(
    generations::Int; method::SubstitutionMethod=SubstitutionMethod()
)
    # Start with an octagonal wreath of 8 rhombi
    star = [SVector(cos(i * π / 4), sin(i * π / 4)) for i in 0:7]

    current_rhombi = []
    for i in 0:7
        # Tile spanned by (e_i, e_{i+1})
        push!(current_rhombi, (i, mod(i+1, 8), SVector(0.0, 0.0)))
    end

    λ = 1 + sqrt(2)

    for _ in 1:generations
        new_rhombi = []
        for (i, j, w) in current_rhombi
            # Grid substitution: e_i -> e_{i-1} + e_i + e_{i+1}
            # The tile (e_i, e_j) is replaced by 9 tiles (e_a, e_b)
            # a in {i-1, i, i+1}, b in {j-1, j, j+1}

            U = [star[mod(i - 1, 8) + 1], star[i + 1], star[mod(i + 1, 8) + 1]]
            V = [star[mod(j - 1, 8) + 1], star[j + 1], star[mod(j + 1, 8) + 1]]

            w_scaled = w * λ

            for (ai, u) in enumerate(U), (bi, v) in enumerate(V)
                pos = w_scaled
                for ak in 1:(ai - 1)
                    pos += U[ak]
                end
                for bk in 1:(bi - 1)
                    pos += V[bk]
                end

                idx_a = mod(i + (ai-2), 8)
                idx_b = mod(j + (bi-2), 8)

                if idx_a != idx_b && mod(abs(idx_a - idx_b), 8) != 4
                    push!(new_rhombi, (idx_a, idx_b, pos))
                end
            end
        end
        current_rhombi = new_rhombi
    end

    # Convert to Tiles and deduplicate
    tile_dict = Dict{NTuple{2,Int},Tile{2,Float64}}()
    for (i, j, w) in current_rhombi
        v1, v2, v3, v4 = w, w + star[i + 1], w + star[i + 1] + star[j + 1], w + star[j + 1]
        c = (v1 + v3) / 2
        key = snap_to_grid(c, SNAP_GRID_EPS)

        # Determine type: Square if |i-j| == 2 or 6, Rhombus otherwise.
        diff = mod(abs(i - j), 8)
        type = (diff == 2 || diff == 6) ? Square() : RhombusAB()

        if !haskey(tile_dict, key)
            tile_dict[key] = Tile{2,Float64}([v1, v2, v3, v4], type, c)
        end
    end

    tiles = collect(values(tile_dict))

    # Collect unique vertices via stable grid snap
    pos_dict = Dict{NTuple{2,Int},SVector{2,Float64}}()
    for tile in tiles
        for v in tile.vertices
            k = snap_to_grid(v, SNAP_GRID_EPS)
            get!(pos_dict, k, v)
        end
    end
    positions = collect(values(pos_dict))

    params = Dict{Symbol,Any}(
        :generations => generations,
        :n_tiles => length(tiles),
        :n_vertices => length(positions),
        :symmetry => 8,
        :window_shape => :none,
    )
    return QuasicrystalData{2,Float64}(AmmannBeenker(), positions, tiles, method, params)
end

"""
    inflate_ammann_beenker_tiles(tiles::Vector{Tile{2, Float64}}, alg=DefaultSubstitution())

Inflation rule for AB: λe₁ = e₁₋₁ + e₁ + e₁₊₁ where λ = 1+√2.

`DefaultSubstitution` and [`AmmannBeenkerInflation`](@ref) are
fully implemented; other algorithms raise `error("not implemented")`.
"""
function inflate_ammann_beenker_tiles(tiles::Vector{Tile{2,Float64}})
    # Robust rule for AB: λe₁ = e₁₋₁ + e₁ + e₁₊₁ where λ = 1+√2
    star = [SVector(cos(i * π / 4), sin(i * π / 4)) for i in 0:7]
    λ = 1 + sqrt(2)

    current_rhombi = []
    for tile in tiles
        v = tile.vertices
        v1 = v[1]
        e1, e2 = v[2] - v1, v[4] - v1
        i = j = -1
        for k in 0:7
            if norm(e1 - star[k + 1]) < STAR_DIRECTION_TOL
                i = k
            end
            if norm(e2 - star[k + 1]) < STAR_DIRECTION_TOL
                j = k
            end
        end
        if i != -1 && j != -1
            push!(current_rhombi, (i, j, v1))
        end
    end

    new_rhombi = []
    for (i, j, w) in current_rhombi
        U = [star[mod(i - 1, 8) + 1], star[i + 1], star[mod(i + 1, 8) + 1]]
        V = [star[mod(j - 1, 8) + 1], star[j + 1], star[mod(j + 1, 8) + 1]]
        w_scaled = w * λ
        for (ai, u) in enumerate(U), (bi, v) in enumerate(V)
            pos = w_scaled
            for ak in 1:(ai - 1)
                pos += U[ak]
            end
            for bk in 1:(bi - 1)
                pos += V[bk]
            end
            idx_a, idx_b = mod(i + (ai-2), 8), mod(j + (bi-2), 8)
            if idx_a != idx_b && mod(abs(idx_a - idx_b), 8) != 4
                push!(new_rhombi, (idx_a, idx_b, pos))
            end
        end
    end

    tile_dict = Dict{Tuple{Int,Int},Tile{2,Float64}}()
    for (i, j, w) in new_rhombi
        v1, v2, v3, v4 = w, w + star[i + 1], w + star[i + 1] + star[j + 1], w + star[j + 1]
        c = (v1 + v3) / 2
        key = (round(Int, c[1]*1e5), round(Int, c[2]*1e5))
        diff = mod(abs(i - j), 8)
        type = (diff == 2 || diff == 6) ? Square() : RhombusAB()
        if !haskey(tile_dict, key)
            tile_dict[key] = Tile{2,Float64}([v1, v2, v3, v4], type, c)
        end
    end
    return collect(values(tile_dict))
end

# ---- Algorithm-dispatched inflation entry points --------------------

function inflate_ammann_beenker_tiles(
    tiles::Vector{Tile{2,Float64}}, ::Union{DefaultSubstitution,AmmannBeenkerInflation}
)
    return inflate_ammann_beenker_tiles(tiles)
end

function inflate_ammann_beenker_tiles(
    ::Vector{Tile{2,Float64}}, alg::AbstractSubstitutionAlgorithm
)
    return error(
        "$(typeof(alg)) is not yet implemented for AmmannBeenker. " *
        "Use DefaultSubstitution() or AmmannBeenkerInflation() instead.",
    )
end

# Single-dispatch on the algorithm: `AmmannBeenkerInflation` is
# AB-specific so this overload is unambiguous.
function inflate_tiles(tiles::Vector{Tile{2,Float64}}, alg::AmmannBeenkerInflation)
    return inflate_ammann_beenker_tiles(tiles, alg)
end
