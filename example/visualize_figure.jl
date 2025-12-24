using QuasiCrystal
using Plots
using LinearAlgebra

# Build quasicrystal data from the library functions
# Build quasicrystal data from the library functions
function build_quasicrystal_data(;
    model::Symbol=:penrose,
    generator::Symbol=:projection,
    radius::Real=3.0,
    generations::Int=4,
    n_points::Int=200,
)
    if model == :penrose
        return if generator == :projection
            generate_penrose_projection(radius)
        else
            generate_penrose_substitution(generations)
        end
    elseif model == :ammann_beenker
        return if generator == :projection
            generate_ammann_beenker_projection(radius)
        else
            generate_ammann_beenker_substitution(generations)
        end
    elseif model == :fibonacci
        return if generator == :projection
            generate_fibonacci_projection(n_points)
        else
            generate_fibonacci_substitution(generations)
        end
    else
        error(
            "Unsupported model: $(model). Choose :penrose, :ammann_beenker, or :fibonacci."
        )
    end
end

# Roughly estimate a reasonable cutoff from the closest pair distance
function estimate_cutoff(qc; sample::Int=300)
    n = min(length(qc.positions), sample)
    min_dist = Inf
    for i in 1:n, j in (i + 1):n
        d = norm(qc.positions[i] - qc.positions[j])
        if d > 1e-4 && d < min_dist
            min_dist = d
        end
    end
    return isfinite(min_dist) ? 1.2 * min_dist : nothing
end

# Ensure bonds exist (build if missing)
function ensure_bonds!(qc; cutoff::Union{Nothing,Real}=nothing)
    local_cutoff = cutoff === nothing ? estimate_cutoff(qc) : cutoff
    local_cutoff === nothing &&
        error("Could not estimate cutoff; pass cutoff=... manually.")
    build_nearest_neighbor_bonds!(qc; cutoff=local_cutoff)
    return qc
end

# Draw positions (and tiles if desired) without bonds
function draw_quasicrystal_picture(;
    model::Symbol=:penrose,
    generator::Symbol=:projection,
    radius::Real=3.0,
    generations::Int=4,
    n_points::Int=200,
    show_tiles::Bool=true,
    tile_kwargs::NamedTuple=NamedTuple(),
    scatter_kwargs::NamedTuple=NamedTuple(),
)
    qc = build_quasicrystal_data(;
        model=model,
        generator=generator,
        radius=radius,
        generations=generations,
        n_points=n_points,
    )
    isempty(qc.positions) && error("No positions generated; increase radius or n_points.")
    dim = length(qc.positions[1])
    if dim == 2
        if show_tiles && !isempty(qc.tiles)
            return visualize_quasicrystal_tiles(qc; tile_kwargs...)
        else
            return visualize_quasicrystal_positions(qc; scatter_kwargs...)
        end
    else
        return visualize_quasicrystal_positions(qc; scatter_kwargs...)
    end
end

# Draw with bonds (builds bonds automatically if missing)
function draw_quasicrystal_with_bonds(;
    model::Symbol=:penrose,
    generator::Symbol=:projection,
    radius::Real=3.0,
    generations::Int=4,
    n_points::Int=200,
    cutoff::Union{Nothing,Real}=nothing,
    scatter_kwargs::NamedTuple=NamedTuple(),
    bond_kwargs::NamedTuple=NamedTuple(),
)
    qc = build_quasicrystal_data(;
        model=model,
        generator=generator,
        radius=radius,
        generations=generations,
        n_points=n_points,
    )
    ensure_bonds!(qc; cutoff=cutoff)
    isempty(qc.bonds) && error("No bonds were built; try increasing cutoff.")
    dim = length(qc.positions[1])
    if dim == 1
        xs = [p[1] for p in qc.positions]
        ys = zeros(length(xs))
        p = Plots.plot(;
            xlabel="Position", yticks=[], legend=false, title="1D quasicrystal bonds"
        )
        for b in qc.bonds
            x1, x2 = xs[b.src], xs[b.dst]
            Plots.plot!(
                p, [x1, x2], [0.0, 0.0]; color=:gray, alpha=0.6, label="", bond_kwargs...
            )
        end
        Plots.scatter!(
            p,
            xs,
            ys;
            marker=:circle,
            markersize=5,
            color=:steelblue,
            label="Sites",
            scatter_kwargs...,
        )
        return p
    else
        xs = [p[1] for p in qc.positions]
        ys = [p[2] for p in qc.positions]
        p = Plots.plot(;
            aspect_ratio=:equal, xlabel="x", ylabel="y", title="2D quasicrystal bonds"
        )
        for b in qc.bonds
            p1 = qc.positions[b.src]
            p2 = qc.positions[b.dst]
            Plots.plot!(
                p,
                [p1[1], p2[1]],
                [p1[2], p2[2]];
                color=:gray,
                alpha=0.6,
                label="",
                bond_kwargs...,
            )
        end
        Plots.scatter!(
            p,
            xs,
            ys;
            marker=:circle,
            markersize=4,
            color=:steelblue,
            label="Sites",
            scatter_kwargs...,
        )
        return p
    end
end
#=
# Example renders (run individually to keep plots clear)
penrose_pic = draw_quasicrystal_picture(;
    model=:penrose,
    generator=:projection,
    radius=3.0,
    show_tiles=false,
    scatter_kwargs=(markersize=4, markercolor=:steelblue),
)
plot!(penrose_pic, title="Penrose Quasicrystal (Projection Method)")
savefig(penrose_pic, "figure1")

ammann_pic = draw_quasicrystal_picture(;
    model=:ammann_beenker,
    generator=:projection,
    radius=3.0,
    show_tiles=false,
    scatter_kwargs=(markersize=3, markercolor=:darkorange),
)
plot!(ammann_pic, title="Ammann-Beenker Quasicrystal (Projection Method)")
savefig(ammann_pic, "figure2")

fibonacci_pic = draw_quasicrystal_picture(;
    model=:fibonacci,
    generator=:substitution,
    generations=8,
    n_points=250,
    scatter_kwargs=(markersize=6, markercolor=:forestgreen),
)
plot!(fibonacci_pic, title="Fibonacci Chain (Substitution Method)")
savefig(fibonacci_pic, "figure3")
=#
# Example renders with bonds (run individually)
penrose_bonds = draw_quasicrystal_with_bonds(;
    model=:penrose,
    generator=:projection,
    radius=10,
    scatter_kwargs=(markersize=4, markercolor=:steelblue),
)

plot!(penrose_pic; title="Penrose Quasicrystal (Projection Method)")
savefig(penrose_pic, "figure1")

ammann_bonds = draw_quasicrystal_with_bonds(;
    model=:ammann_beenker,
    generator=:projection,
    radius=3.5,
    scatter_kwargs=(markersize=3, markercolor=:darkorange),
)
plot!(ammann_pic; title="Ammann-Beenker Quasicrystal (Projection Method)")
savefig(ammann_pic, "figure2")

fibonacci_bonds = draw_quasicrystal_with_bonds(;
    model=:fibonacci,
    generator=:substitution,
    generations=8,
    cutoff=1.7,
    scatter_kwargs=(markersize=6, markercolor=:forestgreen),
)
plot!(fibonacci_pic; title="Fibonacci Chain (Substitution Method)")
savefig(fibonacci_pic, "figure3")
