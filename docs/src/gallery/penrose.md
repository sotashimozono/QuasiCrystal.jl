# Penrose P3 tiling

Penrose P3 is the 5-fold quasicrystal built from the 5D host
lattice $\mathbb{Z}^5$. Its physical star consists of the five
vectors $v_k = (\cos 2\pi k / 5, \sin 2\pi k / 5)$. The
perpendicular space is 3D, and the currently shipped generator
uses a cubic acceptance window
(`BoxWindow{3}(SVector(0.5, 0.5, 0.5))`).

## Real space

The point set is dense but aperiodic, and the distribution of
nearest-neighbour distances shows the hallmark 5-fold structure
even without drawing the rhombi.

![Penrose P3 real space](../assets/figures/penrose_real.png)

```julia
using Plots, LatticeCore, QuasiCrystal
qc = generate_penrose_projection(8.0)
plot_lattice(qc; title="Penrose P3 ($(num_sites(qc)) sites)")
```

## Reciprocal space — Bragg peaks

The diffraction pattern is the clearest visual check of the
whole Fourier stack: the five physical star vectors are a
$C_5$ orbit, and the Bragg peak set must be closed under 72°
rotation about Γ. The plot below shows the ten-fold symmetric
pattern (C₅ combined with the even parity
$I(-k) = I(k)$) that results.

![Penrose P3 diffraction pattern](../assets/figures/penrose_diffraction.png)

```julia
peaks = bragg_peaks(qc; kmax = 8.0, intensity_cutoff = 1e-4)
diffraction_pattern(peaks; log_intensity = true, marker_scale = 14.0,
                    title = "Penrose P3 diffraction (log₁₀ I / I_max)")
```

## What to check visually

- Γ is the single brightest dot at the centre.
- The next ring of five (or ten, counting the $k \to -k$
  partner) bright peaks forms a regular decagon around Γ.
- Sub-rings of fainter peaks appear between the main decagons,
  closed under the same 72° rotation.
- The canonical Penrose pattern has a sharper pentagonal
  window in a 2D Galois-conjugate perpendicular plane — that
  variant is future work, see the follow-up note in
  `src/core/fourier/fourier.jl`.
