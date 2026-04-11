"""
Legacy visualisation helper for quasicrystal point sets.

The long-term plan (see the migration PR discussion) is to replace
this file with a [`LatticeCore`](@ref) package extension that
provides a generic `plot_lattice(::AbstractLattice{D, T})` function
shared between `Lattice2D` and `QuasiCrystal`. Until that lands,
this module keeps a minimal Plots-based scatter routine so
downstream notebooks that call `visualize_quasicrystal_positions`
continue to work.
"""

"""
    visualize_quasicrystal_positions(qc::QuasicrystalData{1, T}; kwargs...)

Scatter a 1D quasicrystal's positions along a line. Requires
`Plots` to be loaded at the call site.
"""
function visualize_quasicrystal_positions(qc::QuasicrystalData{1,T}; kwargs...) where {T}
    positions_1d = [position(qc, i)[1] for i in 1:num_sites(qc)]
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

"""
    visualize_quasicrystal_positions(qc::QuasicrystalData{2, T}; kwargs...)

Scatter a 2D quasicrystal's positions. Requires `Plots`.
"""
function visualize_quasicrystal_positions(qc::QuasicrystalData{2,T}; kwargs...) where {T}
    x_vals = [position(qc, i)[1] for i in 1:num_sites(qc)]
    y_vals = [position(qc, i)[2] for i in 1:num_sites(qc)]

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
