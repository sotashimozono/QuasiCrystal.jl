"""
Acceptance windows and their Fourier transforms.

A cut-and-project quasicrystal is defined by a higher-dimensional
host lattice, a projection onto physical (`E_∥`) and perpendicular
(`E_⊥`) subspaces, and an **acceptance window** `W ⊂ E_⊥`. The
brightness of a Bragg peak at the projected reciprocal lattice
vector `g_∥ = π_∥(g)` is

    I(g) ∝ |Ŵ(g_⊥)|²

where `g_⊥ = π_⊥(g)` and `Ŵ` is the Fourier transform of the
window evaluated at the perpendicular-space momentum.

This file defines the concrete `AcceptanceWindow` subtypes that
`QuasiCrystal` ships, together with their analytic Fourier
transforms. `AcceptanceWindow` itself is an abstract type living
in `LatticeCore` (exported via QuasiCrystal's `using LatticeCore`);
individual quasicrystals attach the window they need to their
`HyperReciprocalLattice`.
"""

# ---- IntervalWindow -------------------------------------------------

"""
    IntervalWindow{T}(half_width::T)

1D interval acceptance window `[-half_width, +half_width]`. Used
by the Fibonacci chain (host lattice `Z²`, perpendicular space `R`).
Its Fourier transform is a scaled sinc:

```math
\\hat{W}(q) \\;=\\; \\int_{-a}^{a} e^{-i q x}\\, dx
         \\;=\\; \\frac{2 \\sin(q a)}{q},
\\qquad a \\equiv \\text{half\\_width}.
```

At `q = 0` we evaluate the limit `2a` directly to avoid a division
by zero.
"""
struct IntervalWindow{T<:AbstractFloat} <: AcceptanceWindow
    half_width::T
end

"""
    window_fourier(w::AcceptanceWindow, q) → T

Evaluate `Ŵ(q)` — the Fourier transform of the acceptance window
at perpendicular-space momentum `q`. Returns a real number (the
scaled window integrals used here are all real because the windows
are centred on the origin and symmetric).
"""
function window_fourier end

function window_fourier(w::IntervalWindow{T}, q::Real) where {T}
    a = w.half_width
    qa = T(q) * a
    if isapprox(qa, zero(T); atol=1e-12)
        return 2 * a
    end
    return 2 * sin(qa) / T(q)
end

# Allow passing an `SVector{1}` — this makes the Fibonacci generator
# code simpler because the perpendicular projection is naturally a
# one-component SVector even though the window is scalar.
function window_fourier(w::IntervalWindow{T}, q::SVector{1,<:Real}) where {T}
    return window_fourier(w, q[1])
end

# ---- BoxWindow ------------------------------------------------------

"""
    BoxWindow{DPerp, T}(half_widths::SVector{DPerp, T})

Axis-aligned hyper-rectangular acceptance window
`∏ᵢ [-half_widthsᵢ, +half_widthsᵢ]` in `DPerp`-dimensional
perpendicular space. Used by the Ammann–Beenker point-set generator
(host `Z⁴`, perpendicular space `R²`, square window).

Its Fourier transform is a product of 1D sincs:

```math
\\hat{W}(q) \\;=\\; \\prod_{i=1}^{D_\\perp} \\frac{2 \\sin(q_i a_i)}{q_i},
\\qquad a_i \\equiv \\text{half\\_widths}_i.
```
"""
struct BoxWindow{DPerp,T<:AbstractFloat} <: AcceptanceWindow
    half_widths::SVector{DPerp,T}
end

function window_fourier(w::BoxWindow{DPerp,T}, q::SVector{DPerp,<:Real}) where {DPerp,T}
    result = one(T)
    for i in 1:DPerp
        a = w.half_widths[i]
        qi = T(q[i])
        qa = qi * a
        if isapprox(qa, zero(T); atol=1e-12)
            result *= 2 * a
        else
            result *= 2 * sin(qa) / qi
        end
    end
    return result
end
