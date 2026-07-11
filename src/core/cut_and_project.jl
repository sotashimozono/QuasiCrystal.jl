"""
Unified cut-and-project framework for the quasicrystal topology markers.

A cut-and-project quasicrystal lives in a higher-dimensional host
lattice `Z^{D_hyper}`, split into a physical subspace `E_∥`
(`D_par`-dimensional) and a perpendicular / internal space `E_⊥`
(`D_perp = D_hyper - D_par`-dimensional). Points of `Z^{D_hyper}`
whose perpendicular projection lands inside an acceptance window
`W ⊂ E_⊥` are kept; their parallel projection forms the
quasicrystal in `E_∥`.

This file exposes the dimensions in a topology-agnostic way so that
downstream code (phason orbits, substitution rules, Fourier analysis
of bond patterns) can be written generically and dispatched on the
lattice marker.

The 1D Fibonacci chain (`FibonacciLattice`) and the 2D
Ammann-Beenker / Penrose tilings share the same primitive API; the
formal correspondence between them is documented in
`notes/cut_and_project_unified.md`.
"""

# ---- Dimensions --------------------------------------------------------------

"""
    CutAndProjectDimensions(D_par, D_hyper, D_perp)

Three integers describing the cut-and-project lift:

- `D_par`   — dimension of physical space `E_∥` (where the
              quasicrystal sits and Bragg peaks live).
- `D_hyper` — dimension of the host integer lattice `Z^{D_hyper}`
              (the lift).
- `D_perp`  — dimension of the perpendicular / internal space
              `E_⊥`, where the acceptance window lives. By
              construction `D_perp = D_hyper - D_par`.
"""
struct CutAndProjectDimensions
    D_par::Int
    D_hyper::Int
    D_perp::Int
end

function CutAndProjectDimensions(; D_par, D_hyper)
    return CutAndProjectDimensions(D_par, D_hyper, D_hyper - D_par)
end

"""
    cut_and_project_dimensions(::AbstractQuasicrystal) → CutAndProjectDimensions

Return the `(D_par, D_hyper, D_perp)` triple for the lattice topology:

| Lattice            | `D_par` | `D_hyper` | `D_perp` |
|--------------------|---------|-----------|----------|
| `FibonacciLattice` | 1       | 2         | 1        |
| `AmmannBeenker`    | 2       | 4         | 2        |
| `PenroseP3`        | 2       | 5         | 3        |

The Penrose `D_perp = 3` is the layout of `generate_penrose_projection`
and `hyper_reciprocal_lattice(::PenroseP3)` in this package today —
not the canonical 2D Galois-conjugate embedding. See
`notes/cut_and_project_unified.md` for the formalism comparison.
"""
function cut_and_project_dimensions end

function cut_and_project_dimensions(::FibonacciLattice)
    return CutAndProjectDimensions(; D_par=1, D_hyper=2)
end
cut_and_project_dimensions(::AmmannBeenker) = CutAndProjectDimensions(; D_par=2, D_hyper=4)
cut_and_project_dimensions(::PenroseP3) = CutAndProjectDimensions(; D_par=2, D_hyper=5)
