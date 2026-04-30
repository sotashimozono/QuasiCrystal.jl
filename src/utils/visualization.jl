"""
Plotting API stubs for quasicrystal visualisation.

The actual `Plots`-backed methods live in
`QuasiCrystalPlotsExt` and are loaded automatically once the
caller issues `using Plots`. This file only declares the public
function names so that

- `using QuasiCrystal` exports a stable symbol regardless of whether
  `Plots` has been loaded yet, and
- `?visualize_quasicrystal_positions` / `?plot_acceptance_window`
  return useful documentation at the REPL even before the extension
  has been triggered.

Calling either function before `using Plots` raises a `MethodError`
just like every other Julia generic — that is intentional and
matches the behaviour of `LatticeCore.fourier_module` /
`structure_factor` when their NUFFT extension is not loaded.
"""

"""
    visualize_quasicrystal_positions(qc::QuasicrystalData{1,T}; kwargs...)
    visualize_quasicrystal_positions(qc::QuasicrystalData{2,T}; kwargs...)

Scatter the physical positions of a quasicrystal point set.
1D lattices render as a horizontal line of dots; 2D lattices render
as a square-aspect scatter plot.

This is a pedagogical / smoke-test helper. The real implementation
ships in the `QuasiCrystalPlotsExt` package extension, so the call
site must `using Plots` first. Without `Plots` loaded, dispatch
falls through and Julia raises a `MethodError`.

For the richer cut-and-project visualisation see
[`plot_acceptance_window`](@ref).
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

This is meant as an **educational** visual: it gives a direct,
geometric picture of the construction behind every `*_projection`
generator. The number of hyper-lattice points enumerated is
deliberately bounded so calls remain interactive in a notebook.

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
