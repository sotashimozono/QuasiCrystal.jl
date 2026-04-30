module QuasiCrystalPlotsExt

using QuasiCrystal
using LatticeCore
using StaticArrays
using LinearAlgebra
using Plots

"""
    QuasiCrystalPlotsExt

Optional Plots-backed visualisation extension. Loaded automatically
once the user issues `using Plots` together with `using QuasiCrystal`.

This extension provides the implementations for the public stubs
declared in `src/utils/visualization.jl`:

- `visualize_quasicrystal_positions(qc; kwargs...)` — physical-space
  scatter of a 1D or 2D quasicrystal point set.
- `plot_acceptance_window(data; show_hyper_points, kwargs...)` —
  cut-and-project acceptance window with optional overlay of the
  hyper-lattice points projected to perpendicular space.
- `plot_tiles(data; palette, ...)` — semantic tile-colouring recipe
  for 2D quasicrystals; each tile is filled with a colour selected
  from the active palette by its `TileType` tag.

Keeping these in the extension keeps `Plots` out of the main module's
load graph; users opt in by `using Plots`.
"""
QuasiCrystalPlotsExt

# ====================================================================
# visualize_quasicrystal_positions
# ====================================================================

function QuasiCrystal.visualize_quasicrystal_positions(
    qc::QuasiCrystal.QuasicrystalData{1,T}; kwargs...
) where {T}
    n = LatticeCore.num_sites(qc)
    positions_1d = [LatticeCore.position(qc, i)[1] for i in 1:n]
    y_vals = zeros(n)

    return Plots.scatter(
        positions_1d,
        y_vals;
        marker=:circle,
        markersize=6,
        label="Sites",
        ylims=(-0.5, 0.5),
        yticks=[],
        xlabel="Position",
        title="1D Quasicrystal",
        legend=:topright,
        kwargs...,
    )
end

function QuasiCrystal.visualize_quasicrystal_positions(
    qc::QuasiCrystal.QuasicrystalData{2,T}; kwargs...
) where {T}
    n = LatticeCore.num_sites(qc)
    x_vals = [LatticeCore.position(qc, i)[1] for i in 1:n]
    y_vals = [LatticeCore.position(qc, i)[2] for i in 1:n]

    return Plots.scatter(
        x_vals,
        y_vals;
        marker=:circle,
        markersize=3,
        label="Vertices",
        aspect_ratio=:equal,
        xlabel="x",
        ylabel="y",
        title="2D Quasicrystal",
        legend=:topright,
        kwargs...,
    )
end

# ====================================================================
# plot_acceptance_window
# ====================================================================

function QuasiCrystal.plot_acceptance_window(
    data::QuasiCrystal.QuasicrystalData;
    show_hyper_points::Bool=true,
    n_hyper::Union{Nothing,Int}=nothing,
    kwargs...,
)
    shape = QuasiCrystal.window_shape(data)
    if shape === :none
        throw(
            ArgumentError(
                "plot_acceptance_window requires a cut-and-project quasicrystal " *
                "with an explicit acceptance window; got window_shape == :none. " *
                "Substitution-generated quasicrystals do not carry a perp-space window.",
            ),
        )
    elseif shape === :unknown
        throw(
            ArgumentError(
                "plot_acceptance_window cannot dispatch on window_shape == :unknown. " *
                "Make sure `data.parameters[:window_shape]` is populated by the generator.",
            ),
        )
    elseif shape === :interval_1d
        return _plot_window_interval(data; show_hyper_points, n_hyper, kwargs...)
    elseif shape === :box_2d
        return _plot_window_box2d(data; show_hyper_points, n_hyper, kwargs...)
    elseif shape === :box_3d
        return _plot_window_box3d(data; show_hyper_points, n_hyper, kwargs...)
    else
        throw(
            ArgumentError(
                "plot_acceptance_window does not yet handle window_shape == :$(shape).",
            ),
        )
    end
end

# ---- helpers --------------------------------------------------------

function _perp_geometry(data::QuasiCrystal.QuasicrystalData)
    hrl = QuasiCrystal.hyper_reciprocal_lattice(data)
    hw = _direct_half_widths(data, hrl)
    return (perp_proj=hrl.perp_proj, half_widths=hw)
end

function _direct_half_widths(
    data::QuasiCrystal.QuasicrystalData, hrl::LatticeCore.HyperReciprocalLattice
)
    DPerp = size(hrl.perp_proj, 1)
    if haskey(data.parameters, :window_size)
        ws = Float64(data.parameters[:window_size])
        return fill(ws, DPerp)
    end
    return fill(0.5, DPerp)
end

_default_n_hyper(::QuasiCrystal.QuasicrystalData{2,T,QuasiCrystal.PenroseP3}) where {T} = 3
_default_n_hyper(::QuasiCrystal.QuasicrystalData) = 4

function _projected_hyper_points(
    perp_proj::SMatrix{DPerp,DHyper,T}, n::Int
) where {DPerp,DHyper,T}
    pts = SVector{DPerp,Float64}[]
    sizehint!(pts, (2n + 1)^DHyper)
    ranges = ntuple(_ -> (-n):n, DHyper)
    for idx in Iterators.product(ranges...)
        v = SVector{DHyper,Float64}(ntuple(k -> Float64(idx[k]), DHyper))
        push!(pts, perp_proj * v)
    end
    return pts
end

function _split_points(
    pts::Vector{SVector{DPerp,Float64}}, half_widths::Vector{Float64}
) where {DPerp}
    inside = SVector{DPerp,Float64}[]
    outside = SVector{DPerp,Float64}[]
    @assert length(half_widths) == DPerp
    for p in pts
        if all(abs(p[d]) <= half_widths[d] for d in 1:DPerp)
            push!(inside, p)
        else
            push!(outside, p)
        end
    end
    return inside, outside
end

# ---- :interval_1d (Fibonacci) --------------------------------------

function _plot_window_interval(
    data::QuasiCrystal.QuasicrystalData;
    show_hyper_points::Bool,
    n_hyper::Union{Nothing,Int},
    kwargs...,
)
    geom = _perp_geometry(data)
    half = geom.half_widths[1]
    n = n_hyper === nothing ? _default_n_hyper(data) : n_hyper

    plt = Plots.plot(;
        xlabel="perpendicular coordinate",
        yticks=[],
        ylims=(-1.0, 1.0),
        title="Acceptance window (interval, Fibonacci)",
        legend=:topright,
        kwargs...,
    )

    Plots.plot!(
        plt, [-half, half], [0.0, 0.0]; linewidth=4, color=:steelblue, label="window"
    )
    Plots.scatter!(
        plt,
        [-half, half],
        [0.0, 0.0];
        marker=:vline,
        markersize=12,
        color=:steelblue,
        label="",
    )

    if show_hyper_points
        pts = _projected_hyper_points(geom.perp_proj, n)
        inside, outside = _split_points(pts, [half])
        if !isempty(inside)
            xs = [p[1] for p in inside]
            ys = zeros(length(xs))
            Plots.scatter!(
                plt, xs, ys; marker=:circle, markersize=4, color=:darkorange, label="inside"
            )
        end
        if !isempty(outside)
            xs = [p[1] for p in outside]
            ys = zeros(length(xs))
            Plots.scatter!(
                plt,
                xs,
                ys;
                marker=:circle,
                markersize=3,
                color=:gray,
                label="outside",
                alpha=0.5,
            )
        end
    end

    return plt
end

# ---- :box_2d (Ammann–Beenker) --------------------------------------

function _plot_window_box2d(
    data::QuasiCrystal.QuasicrystalData;
    show_hyper_points::Bool,
    n_hyper::Union{Nothing,Int},
    kwargs...,
)
    geom = _perp_geometry(data)
    hw = geom.half_widths
    a, b = hw[1], hw[2]
    n = n_hyper === nothing ? _default_n_hyper(data) : n_hyper

    plt = Plots.plot(;
        xlabel="y₁ (perp)",
        ylabel="y₂ (perp)",
        aspect_ratio=:equal,
        title="Acceptance window (square, Ammann–Beenker)",
        legend=:topright,
        kwargs...,
    )

    rect_x = [-a, a, a, -a, -a]
    rect_y = [-b, -b, b, b, -b]
    Plots.plot!(
        plt,
        rect_x,
        rect_y;
        seriestype=:shape,
        fillcolor=:steelblue,
        fillalpha=0.2,
        linecolor=:steelblue,
        linewidth=2,
        label="window",
    )

    if show_hyper_points
        pts = _projected_hyper_points(geom.perp_proj, n)
        inside, outside = _split_points(pts, hw)
        if !isempty(inside)
            xs = [p[1] for p in inside]
            ys = [p[2] for p in inside]
            Plots.scatter!(
                plt, xs, ys; marker=:circle, markersize=4, color=:darkorange, label="inside"
            )
        end
        if !isempty(outside)
            xs = [p[1] for p in outside]
            ys = [p[2] for p in outside]
            Plots.scatter!(
                plt,
                xs,
                ys;
                marker=:circle,
                markersize=3,
                color=:gray,
                label="outside",
                alpha=0.5,
            )
        end
    end

    return plt
end

# ---- :box_3d (Penrose P3) ------------------------------------------

function _plot_window_box3d(
    data::QuasiCrystal.QuasicrystalData;
    show_hyper_points::Bool,
    n_hyper::Union{Nothing,Int},
    kwargs...,
)
    geom = _perp_geometry(data)
    hw = geom.half_widths
    @assert length(hw) == 3
    n = n_hyper === nothing ? _default_n_hyper(data) : n_hyper

    pts = if show_hyper_points
        _projected_hyper_points(geom.perp_proj, n)
    else
        SVector{3,Float64}[]
    end
    inside, outside = _split_points(pts, hw)

    panels = Plots.Plot[]
    for (i, j, label) in ((1, 2, "y₁ vs y₂"), (1, 3, "y₁ vs y₃"), (2, 3, "y₂ vs y₃"))
        a, b = hw[i], hw[j]
        sub = Plots.plot(;
            xlabel="y$(i)", ylabel="y$(j)", aspect_ratio=:equal, title=label, legend=false
        )
        rect_x = [-a, a, a, -a, -a]
        rect_y = [-b, -b, b, b, -b]
        Plots.plot!(
            sub,
            rect_x,
            rect_y;
            seriestype=:shape,
            fillcolor=:steelblue,
            fillalpha=0.2,
            linecolor=:steelblue,
            linewidth=2,
        )
        if !isempty(inside)
            xs = [p[i] for p in inside]
            ys = [p[j] for p in inside]
            Plots.scatter!(sub, xs, ys; marker=:circle, markersize=3, color=:darkorange)
        end
        if !isempty(outside)
            xs = [p[i] for p in outside]
            ys = [p[j] for p in outside]
            Plots.scatter!(
                sub, xs, ys; marker=:circle, markersize=2, color=:gray, alpha=0.4
            )
        end
        push!(panels, sub)
    end

    return Plots.plot(
        panels...;
        layout=(1, 3),
        plot_title="Acceptance window (cube, Penrose P3) — orthogonal projections",
        size=(1200, 400),
        kwargs...,
    )
end

# ====================================================================
# plot_tiles -- semantic tile-type colouring
# ====================================================================

const _PALETTE_DEFAULT = Dict{QuasiCrystal.TileType,Any}(
    QuasiCrystal.FatRhombus()  => RGB(0.96, 0.69, 0.36),
    QuasiCrystal.ThinRhombus() => RGB(0.39, 0.55, 0.78),
    QuasiCrystal.Square()      => RGB(0.34, 0.60, 0.74),
    QuasiCrystal.RhombusAB()   => RGB(0.90, 0.55, 0.35),
)

const _PALETTE_PASTEL = Dict{QuasiCrystal.TileType,Any}(
    QuasiCrystal.FatRhombus()  => RGB(0.99, 0.86, 0.71),
    QuasiCrystal.ThinRhombus() => RGB(0.74, 0.83, 0.93),
    QuasiCrystal.Square()      => RGB(0.74, 0.87, 0.91),
    QuasiCrystal.RhombusAB()   => RGB(0.99, 0.81, 0.71),
)

const _PALETTE_BW = Dict{QuasiCrystal.TileType,Any}(
    QuasiCrystal.FatRhombus()  => RGB(0.85, 0.85, 0.85),
    QuasiCrystal.ThinRhombus() => RGB(0.55, 0.55, 0.55),
    QuasiCrystal.Square()      => RGB(0.85, 0.85, 0.85),
    QuasiCrystal.RhombusAB()   => RGB(0.55, 0.55, 0.55),
)

const _PALETTE_PRESETS = (:default, :pastel, :bw)

function _resolve_palette(palette)
    if palette isa Symbol
        if palette === :default
            return _PALETTE_DEFAULT
        elseif palette === :pastel
            return _PALETTE_PASTEL
        elseif palette === :bw
            return _PALETTE_BW
        else
            throw(
                ArgumentError(
                    "plot_tiles: unknown palette $(repr(palette)). " *
                    "Use one of $(_PALETTE_PRESETS) or pass a " *
                    "Dict{TileType, Color}.",
                ),
            )
        end
    elseif palette isa AbstractDict
        merged = copy(_PALETTE_DEFAULT)
        for (k, v) in palette
            merged[k] = v
        end
        return merged
    else
        throw(
            ArgumentError(
                "plot_tiles: palette must be a preset Symbol " *
                "$(_PALETTE_PRESETS) or a Dict{TileType, Color}; " *
                "got $(typeof(palette)).",
            ),
        )
    end
end

_tile_label(::QuasiCrystal.FatRhombus)  = "Fat rhombus"
_tile_label(::QuasiCrystal.ThinRhombus) = "Thin rhombus"
_tile_label(::QuasiCrystal.Square)      = "Square"
_tile_label(::QuasiCrystal.RhombusAB)   = "Rhombus (AB)"
_tile_label(t::QuasiCrystal.TileType)   = string(QuasiCrystal.tile_type_symbol(t))

@inline function _tile_polygon_xy(tile::QuasiCrystal.Tile{2,T}) where {T}
    v = tile.vertices
    n = length(v)
    xs = Vector{Float64}(undef, n + 1)
    ys = Vector{Float64}(undef, n + 1)
    @inbounds for k in 1:n
        xs[k] = Float64(v[k][1])
        ys[k] = Float64(v[k][2])
    end
    xs[n + 1] = xs[1]
    ys[n + 1] = ys[1]
    return xs, ys
end

function _group_tiles_by_type(tiles)
    groups = Vector{Pair{QuasiCrystal.TileType,Vector{Int}}}()
    seen = Dict{Any,Int}()
    for (i, t) in enumerate(tiles)
        ty = t.type
        key = typeof(ty)
        idx = get(seen, key, 0)
        if idx == 0
            push!(groups, ty => Int[i])
            seen[key] = length(groups)
        else
            push!(groups[idx].second, i)
        end
    end
    return groups
end

function QuasiCrystal.plot_tiles(
    ::QuasiCrystal.QuasicrystalData{1,T}; kwargs...
) where {T}
    throw(
        ArgumentError(
            "plot_tiles: only 2D quasicrystals carry tiles; got D=1. " *
            "Use visualize_quasicrystal_positions for 1D point sets.",
        ),
    )
end

function QuasiCrystal.plot_tiles(
    data::QuasiCrystal.QuasicrystalData{2,T};
    palette=:default,
    show_boundary::Bool=true,
    boundary_color=:black,
    boundary_width::Real=0.5,
    legend::Bool=true,
    title::AbstractString="Quasicrystal tiling",
    kwargs...,
) where {T}
    tiles = data.tiles
    isempty(tiles) && throw(
        ArgumentError(
            "plot_tiles: data.tiles is empty. The projection-method " *
            "generators do not currently populate tiles; use " *
            "generate_*_substitution(...) instead.",
        ),
    )

    pal = _resolve_palette(palette)

    p = Plots.plot(;
        aspect_ratio=:equal,
        xlabel="x",
        ylabel="y",
        title=title,
        legend=legend ? :topright : false,
        kwargs...,
    )

    groups = _group_tiles_by_type(tiles)

    for (tile_type, indices) in groups
        fill_color = get(pal, tile_type, get(_PALETTE_DEFAULT, tile_type, :gray))
        label = _tile_label(tile_type)

        shapes = Vector{Plots.Shape}(undef, length(indices))
        @inbounds for (k, i) in enumerate(indices)
            xs, ys = _tile_polygon_xy(tiles[i])
            shapes[k] = Plots.Shape(xs, ys)
        end

        if show_boundary
            Plots.plot!(
                p,
                shapes;
                fillcolor=fill_color,
                linecolor=boundary_color,
                linewidth=boundary_width,
                label=label,
            )
        else
            Plots.plot!(
                p,
                shapes;
                fillcolor=fill_color,
                linecolor=fill_color,
                linewidth=0,
                label=label,
            )
        end
    end

    return p
end

# ====================================================================
# plot_state -- per-site state colour mapping
# ====================================================================

@inline _project_complex(z, ::Val{:abs2}) = abs2(z)
@inline _project_complex(z, ::Val{:abs}) = abs(z)
@inline _project_complex(z, ::Val{:real}) = real(z)
@inline _project_complex(z, ::Val{:imag}) = imag(z)
@inline _project_complex(z, ::Val{:phase}) = angle(z)

function _project_complex(state::AbstractVector{<:Complex}, mode::Symbol)
    mode in (:abs2, :abs, :real, :imag, :phase) || throw(
        ArgumentError(
            "plot_state: unsupported mode=$(mode). " *
            "Use one of :abs2, :abs, :real, :imag, :phase.",
        ),
    )
    v = Val(mode)
    return [_project_complex(z, v) for z in state]
end

_mode_label(mode::Symbol) = mode === :abs2 ? "|state|^2" :
    mode === :abs ? "|state|" :
    mode === :real ? "Re(state)" :
    mode === :imag ? "Im(state)" :
    mode === :phase ? "arg(state)" :
    "state"

function _plot_state_discrete(
    qc::QuasiCrystal.QuasicrystalData{D,T},
    state::AbstractVector;
    marker_size::Real,
    palette,
    title::AbstractString,
    kwargs...,
) where {D,T}
    n = LatticeCore.num_sites(qc)
    length(state) == n || throw(
        DimensionMismatch(
            "plot_state: length(state)=$(length(state)) != num_sites(qc)=$(n)."
        ),
    )
    levels = sort!(unique(state))
    level_to_idx = Dict(v => k for (k, v) in enumerate(levels))
    color_idx = [level_to_idx[s] for s in state]
    pal = palette === nothing ? Plots.palette(:tab10, max(length(levels), 2)) : palette

    if D == 1
        xs = [LatticeCore.position(qc, i)[1] for i in 1:n]
        ys = zeros(n)
        return Plots.scatter(
            xs,
            ys;
            marker=:circle,
            markersize=marker_size,
            zcolor=color_idx,
            color=pal,
            colorbar=true,
            colorbar_title="state level",
            ylims=(-0.5, 0.5),
            yticks=[],
            xlabel="Position",
            title=title,
            label="",
            kwargs...,
        )
    elseif D == 2
        xs = [LatticeCore.position(qc, i)[1] for i in 1:n]
        ys = [LatticeCore.position(qc, i)[2] for i in 1:n]
        return Plots.scatter(
            xs,
            ys;
            marker=:circle,
            markersize=marker_size,
            zcolor=color_idx,
            color=pal,
            colorbar=true,
            colorbar_title="state level",
            aspect_ratio=:equal,
            xlabel="x",
            ylabel="y",
            title=title,
            label="",
            kwargs...,
        )
    else
        throw(ArgumentError("plot_state: only D=1 and D=2 lattices are supported (got D=$(D))."))
    end
end

function _plot_state_continuous(
    qc::QuasiCrystal.QuasicrystalData{D,T},
    values::AbstractVector{<:Real};
    marker_size::Real,
    colormap,
    cbar_title::AbstractString,
    title::AbstractString,
    kwargs...,
) where {D,T}
    n = LatticeCore.num_sites(qc)
    length(values) == n || throw(
        DimensionMismatch(
            "plot_state: length(state)=$(length(values)) != num_sites(qc)=$(n)."
        ),
    )
    if D == 1
        xs = [LatticeCore.position(qc, i)[1] for i in 1:n]
        ys = zeros(n)
        return Plots.scatter(
            xs,
            ys;
            marker=:circle,
            markersize=marker_size,
            zcolor=values,
            color=colormap,
            colorbar=true,
            colorbar_title=cbar_title,
            ylims=(-0.5, 0.5),
            yticks=[],
            xlabel="Position",
            title=title,
            label="",
            kwargs...,
        )
    elseif D == 2
        xs = [LatticeCore.position(qc, i)[1] for i in 1:n]
        ys = [LatticeCore.position(qc, i)[2] for i in 1:n]
        return Plots.scatter(
            xs,
            ys;
            marker=:circle,
            markersize=marker_size,
            zcolor=values,
            color=colormap,
            colorbar=true,
            colorbar_title=cbar_title,
            aspect_ratio=:equal,
            xlabel="x",
            ylabel="y",
            title=title,
            label="",
            kwargs...,
        )
    else
        throw(ArgumentError("plot_state: only D=1 and D=2 lattices are supported (got D=$(D))."))
    end
end

function QuasiCrystal.plot_state(
    qc::QuasiCrystal.QuasicrystalData{D,T},
    state::AbstractVector{<:Union{Bool,Integer}};
    colormap=:viridis,
    marker_size::Real=4,
    palette=nothing,
    title::AbstractString="QuasiCrystal state (discrete)",
    kwargs...,
) where {D,T}
    return _plot_state_discrete(
        qc,
        state;
        marker_size=marker_size,
        palette=palette,
        title=title,
        kwargs...,
    )
end

function QuasiCrystal.plot_state(
    qc::QuasiCrystal.QuasicrystalData{D,T},
    state::AbstractVector{<:Real};
    colormap=:viridis,
    marker_size::Real=4,
    title::AbstractString="QuasiCrystal state",
    kwargs...,
) where {D,T}
    return _plot_state_continuous(
        qc,
        state;
        marker_size=marker_size,
        colormap=colormap,
        cbar_title="state",
        title=title,
        kwargs...,
    )
end

function QuasiCrystal.plot_state(
    qc::QuasiCrystal.QuasicrystalData{D,T},
    state::AbstractVector{<:Complex};
    mode::Symbol=:abs2,
    colormap=:viridis,
    marker_size::Real=4,
    title::AbstractString="QuasiCrystal state",
    kwargs...,
) where {D,T}
    values = _project_complex(state, mode)
    return _plot_state_continuous(
        qc,
        values;
        marker_size=marker_size,
        colormap=colormap,
        cbar_title=_mode_label(mode),
        title=title,
        kwargs...,
    )
end

end # module QuasiCrystalPlotsExt
