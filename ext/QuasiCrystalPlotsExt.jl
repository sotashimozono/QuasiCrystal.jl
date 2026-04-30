module QuasiCrystalPlotsExt

using QuasiCrystal
using LatticeCore
using StaticArrays
using LinearAlgebra
using Plots

"""
    QuasiCrystalPlotsExt

Plots-backed visualisation extension. Provides

- `visualize_quasicrystal_positions(qc; kwargs...)` — physical-space
  scatter of a 1D or 2D quasicrystal point set.
- `plot_acceptance_window(data; show_hyper_points, kwargs...)` —
  cut-and-project acceptance window with optional overlay of the
  hyper-lattice points projected to perpendicular space.

Loaded automatically once `using Plots` runs in the same session as
`using QuasiCrystal`. The function stubs themselves are declared in
`src/utils/visualization.jl` so that the public API is visible even
when `Plots` has not been loaded yet.
"""
QuasiCrystalPlotsExt

# ====================================================================
# visualize_quasicrystal_positions
# ====================================================================

function QuasiCrystal.visualize_quasicrystal_positions(
    qc::QuasiCrystal.QuasicrystalData{1,T}; kwargs...
) where {T}
    positions_1d = [LatticeCore.position(qc, i)[1] for i in 1:LatticeCore.num_sites(qc)]
    y_vals = zeros(length(positions_1d))

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
    x_vals = [LatticeCore.position(qc, i)[1] for i in 1:LatticeCore.num_sites(qc)]
    y_vals = [LatticeCore.position(qc, i)[2] for i in 1:LatticeCore.num_sites(qc)]

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

# ---- public entry point --------------------------------------------

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
                "plot_acceptance_window does not yet handle window_shape == :$(shape)."
            ),
        )
    end
end

# ---- helpers --------------------------------------------------------

# Pull the perp_proj matrix and a half-width vector for the
# perpendicular acceptance window. The perp projection comes from
# `hyper_reciprocal_lattice` (the canonical source of perp
# geometry); the half-widths come from the generator's own
# direct-space filter so that the "inside" / "outside" colouring
# in the recipe matches the integer points the generator actually
# accepts.
function _perp_geometry(data::QuasiCrystal.QuasicrystalData)
    hrl = QuasiCrystal.hyper_reciprocal_lattice(data)
    hw = _direct_half_widths(data, hrl)
    return (perp_proj=hrl.perp_proj, half_widths=hw)
end

# AB and Penrose store `:window_size` and use the same numeric
# value on every perp axis. Fibonacci hardcodes `acceptance_width
# = 1.0` (half-width 0.5) and does not stash it in `parameters`,
# so we fall back on the generator's literal value there. If a
# user supplies a custom generator that does record
# `:window_size`, we honour it.
function _direct_half_widths(
    data::QuasiCrystal.QuasicrystalData, hrl::LatticeCore.HyperReciprocalLattice
)
    DPerp = size(hrl.perp_proj, 1)
    if haskey(data.parameters, :window_size)
        ws = Float64(data.parameters[:window_size])
        return fill(ws, DPerp)
    end
    # Fibonacci-style generators with no :window_size key default
    # to the orthonormal-frame half-width 0.5, matching
    # `generate_fibonacci_projection`'s `abs(pos_perp) <= 0.5`
    # filter. If a custom generator deviates from this convention
    # it should populate `:window_size` explicitly.
    return fill(0.5, DPerp)
end

# Default integer enumeration box. Penrose has DHyper = 5, so we keep
# the radius small to avoid 11^5 ≈ 1.6e5 points dominating the figure.
_default_n_hyper(::QuasiCrystal.QuasicrystalData{2,T,QuasiCrystal.PenroseP3}) where {T} = 3
_default_n_hyper(::QuasiCrystal.QuasicrystalData) = 4

# Enumerate every integer point in [-n, n]^DHyper exactly once,
# returning their perpendicular-space projections. We strip the
# `2π` factor that lives on `hyper_basis` because we want the
# *direct-space* integer host points (Z^DHyper), not the reciprocal
# ones — the educational picture is "which n ∈ Z^DHyper survive the
# perp-window filter".
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

# Return (inside, outside) splits — Vector{SVector{DPerp,Float64}}
# pairs — using axis-aligned half-widths.
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

    # Window: thick segment on y = 0 with end caps.
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

    # Filled rectangle.
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
#
# We cannot draw a true 3D acceptance cube cleanly inside the
# default Plots GR backend without surface meshes, so the recipe
# falls back on the classic "three orthogonal projections" view —
# the y₁y₂, y₁y₃, and y₂y₃ planes side by side — each carrying the
# corresponding 2D rectangular window cross-section.

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

end # module QuasiCrystalPlotsExt
