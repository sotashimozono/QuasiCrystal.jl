# Cut-and-project framework: 1D ↔ 2D correspondence

QuasiCrystal.jl ships three concrete topology markers
(`FibonacciLattice`, `AmmannBeenker`, `PenroseP3`) that all share the
same cut-and-project mathematical structure. This note pins down the
formal correspondence between the 1D Fibonacci chain and the 2D
Ammann-Beenker / Penrose tilings so downstream code can dispatch on
the lattice marker and treat dimensions / inflation / phason orbits
uniformly.

## 1. The cut-and-project setup

A cut-and-project quasicrystal is built from:

* a **host lattice** `Z^{D_hyper}` (the lift),
* a splitting `R^{D_hyper} = E_∥ ⊕ E_⊥` into physical space
  (`D_par`-dimensional) and internal / perpendicular space
  (`D_perp = D_hyper - D_par`-dimensional),
* an **acceptance window** `W ⊂ E_⊥`.

The quasicrystal is

```
Λ = { π_∥(n) : n ∈ Z^{D_hyper},  π_⊥(n) ∈ W }
```

i.e. take every host point whose perpendicular shadow lands inside
`W`, and project it to `E_∥`.

The three shipped topologies populate this triple:

| Lattice            | `D_par` | `D_hyper` | `D_perp` | Window in `E_⊥`              |
|--------------------|---------|-----------|----------|------------------------------|
| `FibonacciLattice` | 1       | 2         | 1        | interval `[-½, ½] / √(1+φ⁻²)` |
| `AmmannBeenker`    | 2       | 4         | 2        | regular octagon              |
| `PenroseP3`        | 2       | 5         | 3        | box (current code; see §5)   |

The Fibonacci interval and the AB octagon are documented in the
`IntervalWindow` / `BoxWindow` definitions in
`src/core/fourier/window.jl`. The Penrose `D_perp = 3` reflects the
current Z⁵ embedding (`generate_penrose_projection`), not the canonical
five-fold 2D Galois-conjugate window — see §5 for the discrepancy.

## 2. Inflation symmetry — the unifying linear map

Each quasicrystal carries an integer linear map `M : Z^{D_hyper} → Z^{D_hyper}`
that respects the `E_∥ / E_⊥` splitting. Restricted to `E_∥` it acts as
the **Perron–Frobenius eigenvalue** `λ_PF` of the substitution
(physical inflation); restricted to `E_⊥` it acts as its **Galois
conjugate** `λ_GC` (contracting in `E_⊥`, since `|λ_GC| < 1`).

| Lattice            | `λ_PF`             | `λ_GC`              | Algebraic class               |
|--------------------|--------------------|---------------------|-------------------------------|
| `FibonacciLattice` | `φ ≈ 1.618`        | `-1/φ ≈ -0.618`     | Golden mean, deg 2 PV root    |
| `AmmannBeenker`    | `σ = 1+√2 ≈ 2.414` | `-1/σ ≈ -0.414`     | Silver mean, deg 2 PV root    |
| `PenroseP3`        | `φ² ≈ 2.618`       | `-1/φ² ≈ -0.382`    | Squared golden, deg 2 PV root |

Each `λ_PF` is a PV (Pisot–Vijayaraghavan) number — algebraic integer
> 1 whose Galois conjugates all lie strictly inside the unit disk.
The PV property is what makes the diffraction pattern pure-point
(Meyer's theorem; see Baake & Grimm *Aperiodic Order* Vol. 1, Ch. 9).

**Fibonacci `M`** (verified in `core/model/fibonacci_phason.jl`):

```
M = ((1, 1),
     (1, 0))           eigvals (φ, -1/φ)   on Z²
```

`M · (1, 1)ᵀ = (2, 1)`, generating the Fibonacci numbers `(F_{n+1}, F_n)`.

**Ammann-Beenker `M`** (TODO): the canonical `Z⁴` silver-mean
inflation. Construction outline:

* Beenker (1982) showed the AB substitution `σ_AB` on the 8-fold star
  basis `e_k = (cos((k-1)π/4), sin((k-1)π/4))` for `k = 1,…,4` (the
  first half of the 8 directions) acts as a 4×4 integer matrix with
  characteristic polynomial `(x² - 2x - 1)² = x⁴ - 4x³ + 2x² + 4x + 1`,
  giving eigvals `(σ, σ, -1/σ, -1/σ)`.
* Trace = 4, determinant = 1; both `E_∥` and `E_⊥` 2-dimensional, so
  eigvals are 2-fold degenerate on each subspace.
* The explicit matrix in the canonical basis is documented in
  Baake & Grimm *Aperiodic Order* Vol. 1, §6.3, and in Fuchs, Mosseri,
  Vidal (2018) "Hyperuniformity and Fourier-modular structure of
  inflation tilings" (cited as [15] in arXiv 2310.20517).
* The naïve "neighbor-sum" matrix `M = I + R + R^{-1}` where `R` is
  the 4-cycle has eigvals `(3, 1, 1, -1)` — these are correct for a
  *different* (rational-Perron) substitution but NOT the silver-mean
  inflation; do not confuse the two.

**Penrose P3 `M`** (TODO): the canonical `Z⁵` golden-squared inflation
acting on the five-fold star `e_k = (cos((k-1)·2π/5), sin((k-1)·2π/5))`
for `k = 1,…,5`. Characteristic polynomial
`(x² - 3x + 1)² (x - 1)` — see Senechal *Quasicrystals and Geometry*
(1995) Ch. 8 or de Bruijn (1981) for the explicit form. Eigvals
`(φ², φ², -1/φ², -1/φ², 1)`; the unit eigenvalue corresponds to the
`(1,1,1,1,1)` diagonal direction of Z⁵ that sums to a constant under
five-fold cyclic permutation (de Bruijn's "trivial axis").

These two `M` definitions land in a follow-up commit; the Phase α
slice exposes only `cut_and_project_dimensions` so that downstream
code can already dispatch on the lattice marker for dimension-aware
operations (e.g. constructing the right-sized `SVector` for a phason
orbit point).

## 3. Phason orbits — vector-valued for `D_perp ≥ 2`

The 1D phason orbit on the Fibonacci chain is

```
θ_i = mod(θ_0 + i · α, 1) ∈ R/Z,    α = 1/φ
```

(see `phason_orbit_at(::FibonacciLattice, i; θ0, α)` in
`core/model/fibonacci_phason.jl`).

For 2D Ammann-Beenker the analogue is

```
θ_n ∈ R²/Z²,   θ_{n+1} = θ_n + α  (mod Z²)
```

where `α ∈ R²` is irrational (lives in the silver-mean rank-2 module
`Z + Zσ` per axis). For Penrose, `θ_n ∈ R³/Z³` (or `R²/Z²` in the
canonical 2D Galois-conjugate embedding).

The shipped API will be

```julia
phason_orbit_at(::AmmannBeenker, n::SVector{2,Int};
                θ0::SVector{2,Float64}, α::SMatrix{2,2,Float64})
    → SVector{2,Float64}
```

with the inflation matrix's perpendicular block providing the natural
`α` (the Galois conjugate of the parallel projection). This lands in
the same follow-up commit as the inflation matrices.

## 4. 1D ↔ 2D as restrictions of the same framework

Both the 1D Fibonacci chain and the 2D Ammann-Beenker tiling sit in
the same cut-and-project framework, differing only in the dimensions
`(D_par, D_hyper, D_perp)` and the choice of inflation algebraic
integer (golden vs silver mean).

Concretely:

* Fibonacci is a 1D quasicrystal whose phason torus is `R/Z`. Its
  inflation matrix is the smallest nontrivial PV matrix: 2×2 with
  golden eigvals.
* Ammann-Beenker is a 2D quasicrystal whose phason torus is `R²/Z²`.
  Its inflation matrix is the canonical 4×4 silver-mean matrix.

The `cut_and_project_dimensions(::AbstractQuasicrystal)` accessor is
the type-level switch that downstream code uses to allocate the right
`SVector`/`SMatrix` sizes, choose between scalar and vector phason
arithmetic, and dispatch the right rank of `inflation_matrix`.

## 5. Penrose `D_perp = 3` caveat

The current `generate_penrose_projection` and
`hyper_reciprocal_lattice(::PenroseP3)` use a `Z⁵`-host embedding with
a 3D perpendicular space. The canonical Penrose construction
(de Bruijn 1981) uses a 2D Galois-conjugate perpendicular space
because the trivial `(1,1,1,1,1)` direction of `Z⁵` is constant on
the orbit and can be quotiented out. The 3D embedding is a valid lift,
but it duplicates information in the trivial direction.

A future revision will add a `PenroseP3Canonical <: AbstractQuasicrystal{2}`
marker that quotients out the trivial axis and exposes the 2D `D_perp`
canonical form. The current `PenroseP3` marker is preserved as the
`D_perp = 3` non-canonical embedding for compatibility with the
existing Bragg-peak enumeration code in `core/fourier/fourier.jl`.

## References

* M. Baake, U. Grimm. *Aperiodic Order Vol. 1: A Mathematical Invitation.*
  Cambridge University Press, 2013. DOI:[10.1017/CBO9781139025256][bg].
  Chapters 4 (Fibonacci), 6 (Ammann-Beenker), 7 (Penrose).
* L. Boyle, P. J. Steinhardt. "Self-similar one-dimensional
  quasilattices." [arXiv:1608.08220][bs], 2016.
* C. L. Fuchs, J. C. Mosseri, J. Vidal. "Hyperuniformity and
  Fourier-modular structure of inflation tilings."
  (Cited as [15] in arXiv 2310.20517; lookup via doiget.)
* D. Roca. "Hyperuniformity and number rigidity of inflation tilings."
  [arXiv:2310.20517][roca], 2023. References Baake-Grimm-Mañibo for AB.
* N. G. de Bruijn. "Algebraic theory of Penrose's non-periodic
  tilings of the plane." Nederl. Akad. Wetensch. Proc. Ser. A 84
  (1981) 39-52, 53-66. (Two-paper foundation of the canonical Z⁵
  Penrose embedding.)
* M. Senechal. *Quasicrystals and Geometry.* Cambridge, 1995.
  Chapter 8 (Penrose) and Appendix A (cut-and-project formalism).

[bg]: https://doi.org/10.1017/CBO9781139025256
[bs]: https://arxiv.org/abs/1608.08220
[roca]: https://arxiv.org/abs/2310.20517
