# QuasiCrystal.jl

`QuasiCrystal.jl` provides cut-and-project quasicrystal point sets
(Fibonacci chain, Ammann–Beenker tiling, Penrose P3 tiling),
together with a unified Fourier-analysis pipeline that returns a
`BraggPeakSet` directly usable with `LatticeCore.structure_factor`
and the `momentum_lattice` / `fourier_module` trait machinery.

## Quick start

```julia
using Plots, LatticeCore, QuasiCrystal

qc    = generate_penrose_projection(8.0)
peaks = bragg_peaks(qc; kmax = 8.0, intensity_cutoff = 1e-4)

plot_lattice(qc)                              # real space
diffraction_pattern(peaks; log_intensity=true) # reciprocal space
```

## Where to look next

- The [Cut-and-project construction](guide/cut_and_project.md)
  page collects the mathematical background — enough to read
  the rest of the docs without pulling out a textbook.
- The [Fourier analysis pipeline](guide/fourier_analysis.md)
  page walks through how `bragg_peaks` is built from
  `hyper_reciprocal_lattice` and an acceptance window.
- The **Gallery** pages show real-space and reciprocal-space
  plots for each shipped quasicrystal. These figures are
  regenerated from the live source tree on every doc build, so
  they double as an end-to-end visual regression test of the
  Fourier pipeline.
- The [API](api.md) page lists every exported function.
