"""
    generate_figures.jl

Build the visual-verification figures that are embedded in the
QuasiCrystal documentation. Called from `docs/make.jl` before
`makedocs`, so every doc build re-generates the images from the
current source tree. If the Fourier pipeline ever regresses, the
resulting PNGs will look visibly wrong (or crash) and the build
will fail — which is exactly the "verify visually in the docs"
loop we want.

All figures land in `docs/src/assets/figures/`.
"""

ENV["GKSwstype"] = "100"   # headless gr backend

using Plots
using LatticeCore
using QuasiCrystal

const FIG_DIR = joinpath(@__DIR__, "src", "assets", "figures")
mkpath(FIG_DIR)

function savepng(p, name)
    path = joinpath(FIG_DIR, name)
    savefig(p, path)
    @info "wrote $(relpath(path, @__DIR__))"
    return path
end

# ---- 1. Fibonacci chain --------------------------------------------

@info "Fibonacci"
qc_fib = generate_fibonacci_projection(34)
savepng(
    plot_lattice(qc_fib; title="Fibonacci chain ($(num_sites(qc_fib)) sites)"),
    "fibonacci_real.png",
)

fib_peaks = bragg_peaks(qc_fib; kmax=20.0, intensity_cutoff=1e-4)
savepng(
    diffraction_pattern(
        fib_peaks;
        title="Fibonacci diffraction (log₁₀ I / I_max)",
        log_intensity=true,
    ),
    "fibonacci_diffraction.png",
)

# ---- 2. Ammann–Beenker ---------------------------------------------

@info "Ammann-Beenker"
qc_ab = generate_ammann_beenker_projection(8.0)
savepng(
    plot_lattice(qc_ab; title="Ammann–Beenker ($(num_sites(qc_ab)) sites)"),
    "ammann_beenker_real.png",
)

ab_peaks = bragg_peaks(qc_ab; kmax=8.0, intensity_cutoff=1e-4)
savepng(
    diffraction_pattern(
        ab_peaks;
        title="Ammann–Beenker diffraction (log₁₀ I / I_max)",
        log_intensity=true,
        marker_scale=14.0,
    ),
    "ammann_beenker_diffraction.png",
)

# ---- 3. Penrose P3 -------------------------------------------------

@info "Penrose P3"
qc_pen = generate_penrose_projection(8.0)
savepng(
    plot_lattice(qc_pen; title="Penrose P3 ($(num_sites(qc_pen)) sites)"),
    "penrose_real.png",
)

pen_peaks = bragg_peaks(qc_pen; kmax=8.0, intensity_cutoff=1e-4)
savepng(
    diffraction_pattern(
        pen_peaks;
        title="Penrose P3 diffraction (log₁₀ I / I_max)",
        log_intensity=true,
        marker_scale=14.0,
    ),
    "penrose_diffraction.png",
)

@info "All figures regenerated." sites_fib = num_sites(qc_fib) sites_ab = num_sites(
    qc_ab
) sites_pen = num_sites(qc_pen) peaks_fib = num_k_points(fib_peaks) peaks_ab = num_k_points(
    ab_peaks
) peaks_pen = num_k_points(pen_peaks)
