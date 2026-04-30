"""
Plots-based visualisation entry points for `QuasicrystalData`.

This file declares the public visualisation API as bare
`function ... end` stubs. Real implementations live in the
`QuasiCrystalPlotsExt` package extension
(`ext/QuasiCrystalPlotsExt.jl`) and load automatically once the user
issues `using Plots` together with `using QuasiCrystal`.

Keeping the stubs here lets QuasiCrystal export the symbols without
forcing `Plots` (a comparatively heavy dependency) into every
downstream load. Calling any of them before `Plots` is in scope
produces the standard "no method matching" error from Julia, since
no methods are defined until the extension activates.

# API

| Function                              | Purpose                                    |
|---------------------------------------|--------------------------------------------|
| [`visualize_quasicrystal_positions`]  | Scatter the site positions (1D / 2D).      |
| [`plot_acceptance_window`]            | Cut-and-project window + hyper points.     |
| [`plot_tiles`]                        | Polygon-fill 2D tiles, colour by type.     |
| [`plot_state`]                        | Colour-mapped per-site state vector.       |
"""

"""
    visualize_quasicrystal_positions(qc::QuasicrystalData{D, T}; kwargs...)

Scatter the vertex positions of `qc`. Implemented in
`QuasiCrystalPlotsExt`; requires `using Plots` at the call site.

Supported dimensions are `D in {1, 2}`. Extra `kwargs` are forwarded
to `Plots.scatter`.
"""
function visualize_quasicrystal_positions end

"""
    plot_acceptance_window(data::QuasicrystalData; show_hyper_points=true, kwargs...)

Visualise the **acceptance window** of a cut-and-project
quasicrystal `data`, optionally overlaying the integer
hyper-lattice points projected onto perpendicular space.

The recipe dispatches on [`window_shape`](@ref):

| `window_shape(data)` | Drawn primitive                                       |
|----------------------|-------------------------------------------------------|
| `:interval_1d`       | Horizontal segment with the window endpoints          |
| `:box_2d`            | Filled rectangle in the 2D perpendicular plane        |
| `:box_3d`            | Three orthogonal projections of the perpendicular cube |
| `:none` / `:unknown` | Raises `ArgumentError`                                |

When `show_hyper_points = true` the recipe enumerates a bounded
neighbourhood of the host lattice `Z^DHyper`, projects every point
through `perp_proj`, and scatters them on top of the window:
points inside the window (the ones that survive the cut-and-project
filter) are drawn in one colour, points outside in a contrasting
colour.

`plot_acceptance_window` is shipped by the `QuasiCrystalPlotsExt`
package extension; the caller must `using Plots` to trigger the
extension before invoking the function.

# Keyword arguments

- `show_hyper_points::Bool = true` — overlay the projected
  `Z^DHyper` lattice. Set to `false` to draw the window alone.
- `n_hyper::Int` — half-width of the integer enumeration box per
  axis. Defaults to a small value (3 for high-D Penrose, 4
  otherwise) to keep the call cheap.
- Any other keyword is forwarded to the underlying `Plots`
  primitive (`title`, `size`, `markersize`, ...).
"""
function plot_acceptance_window end

"""
    plot_tiles(data::QuasicrystalData{2, T};
               palette = :default,
               show_boundary::Bool = true,
               boundary_color = :black,
               boundary_width::Real = 0.5,
               legend::Bool = true,
               kwargs...) -> Plots.Plot

Render the 2D quasicrystalline tiling carried by `data` as a list of
filled polygons, colouring each tile according to its semantic
[`TileType`](@ref). Implemented in `QuasiCrystalPlotsExt`; requires
`using Plots` at the call site.

# Tile-type colouring

Each tile is filled with the colour assigned to its `TileType`
singleton. The `palette` keyword controls the assignment:

- `palette = :default` -- a balanced two-tone scheme
  (warm fat / cool thin for Penrose,
  cool square / warm rhombus for Ammann-Beenker).
- `palette = :pastel` -- softer pastel variants of the defaults.
- `palette = :bw` -- light grey vs. dark grey for printer-friendly
  black-and-white figures.
- `palette::AbstractDict{<:TileType, <:Any}` -- explicit mapping;
  any tile types not present in the dict fall back to the `:default`
  palette, so a *partial* dict only overrides the listed types.

# Boundary

When `show_boundary` is `true` (default), the tile edges are stroked
with `boundary_color` (default `:black`) and `boundary_width`
(default `0.5`). Set `show_boundary = false` for a flat tile-fill
render with no edges.

# Legend

When `legend = true` (default) one legend entry per distinct
`TileType` is emitted (`"Fat rhombus"`, `"Thin rhombus"`, `"Square"`,
`"Rhombus (AB)"`).

Returns the `Plots.Plot` object so callers can compose with
`plot!`, `savefig`, etc.

# Examples

```julia
using Plots, QuasiCrystal
qc = generate_penrose_substitution(3)
plot_tiles(qc)                                # default palette
plot_tiles(qc; palette = :pastel)             # pastel palette
plot_tiles(qc; palette = :bw, legend = false) # B/W, no legend

# Custom palette: only override the fat rhombus colour
plot_tiles(qc; palette = Dict(FatRhombus() => :crimson))

# Ammann-Beenker tiling, no edge strokes
ab = generate_ammann_beenker_substitution(3)
plot_tiles(ab; show_boundary = false)
```

# Errors

- `ArgumentError` if `D != 2` (no tiles in 1D quasicrystals).
- `ArgumentError` if `data.tiles` is empty (the projection-method
  generators do not currently populate tiles; use a substitution
  generator instead).
- `ArgumentError` for an unknown `palette` symbol.
"""
function plot_tiles end

"""
    plot_state(data::QuasicrystalData, state;
               colormap = :viridis,
               marker_size::Real = 4,
               palette = nothing,
               mode::Symbol = :abs2,
               title::AbstractString = "QuasiCrystal state",
               kwargs...) -> Plots.Plot

Visualise a per-site state vector on top of the quasicrystal scatter.
`state[i]` is mapped onto the marker colour of vertex `i`. Implemented
in `QuasiCrystalPlotsExt`; requires `using Plots` at the call site.

Three input forms are supported, dispatched on `eltype(state)`:

* Real-valued state -- `state::AbstractVector{<:Real}`. Continuous
  `colormap` (default `:viridis`).
* Discrete state -- `state::AbstractVector{<:Union{Bool, Integer}}`.
  Categorical palette (`:tab10` default; override via `palette`).
* Complex state -- `state::AbstractVector{<:Complex}`. Projected to
  a real scalar via `mode in (:abs2, :abs, :real, :imag, :phase)`
  (default `:abs2`) before colour-mapping.

# Examples

```julia
using QuasiCrystal, Plots

qc = generate_penrose_substitution(3)
state = [sin(p[1]) * cos(p[2]) for p in get_positions(qc)]
plot_state(qc, state; colormap = :plasma, marker_size = 5)
```

```julia
fib = generate_fibonacci_substitution(6)
occ = [iseven(i) for i in 1:num_sites(fib)]
plot_state(fib, occ; marker_size = 6)
```
"""
function plot_state end
