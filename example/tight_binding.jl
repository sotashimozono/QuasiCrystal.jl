"""
Tight-binding model implementation for quasicrystals.

This module implements quantum mechanics calculations on quasicrystal structures,
focusing on the tight-binding approximation commonly used in condensed matter physics.
"""

using QuasiCrystal
using LinearAlgebra
using SparseArrays

"""
    TightBindingModel{T}
Represents a tight-binding Hamiltonian on a quasicrystal.

# Fields
- `H::SparseMatrixCSC{T,Int}`: Hamiltonian matrix
- `positions::Vector{Vector{Float64}}`: Site positions
- `n_sites::Int`: Number of sites
- `is_hermitian::Bool`: Whether the Hamiltonian is Hermitian
"""
struct TightBindingModel{T}
  H::SparseMatrixCSC{T,Int}
  positions::Vector{Vector{Float64}}
  n_sites::Int
  is_hermitian::Bool

  function TightBindingModel(
    H::SparseMatrixCSC{T,Int}, positions::Vector{Vector{Float64}}
  ) where {T}
    n = size(H, 1)
    @assert n == length(positions) "Hamiltonian size must match number of positions"
    is_herm = ishermitian(H)
    return new{T}(H, positions, n, is_herm)
  end
end

"""
    build_fibonacci_tight_binding(qc_data::QuasicrystalData{1,T}; 
                                   t::Real=1.0, 
                                   onsite_modulation::Real=0.0) where T
Build tight-binding model for 1D Fibonacci lattice.

# Arguments
- `qc_data`: Fibonacci quasicrystal data
- `t`: Hopping parameter (default: 1.0)
- `onsite_modulation`: On-site energy modulation based on spacing type (default: 0.0)

# Returns
- `TightBindingModel`: Tight-binding model with nearest-neighbor hopping
"""
function build_fibonacci_tight_binding(
  qc_data::QuasicrystalData{1,T}; t::Real=1.0, onsite_modulation::Real=0.0
) where {T}
  n = length(qc_data.positions)

  # Extract 1D positions and sort
  positions_1d = [p[1] for p in qc_data.positions]
  perm = sortperm(positions_1d)
  sorted_positions = [qc_data.positions[i] for i in perm]
  sorted_1d = positions_1d[perm]

  # Build Hamiltonian with nearest-neighbor hopping
  I_idx = Int[]
  J_idx = Int[]
  V_vals = ComplexF64[]

  # On-site energies (can vary based on local environment)
  if onsite_modulation != 0.0
    # Determine spacing type for each site based on neighbors
    for i in 1:n
      onsite = 0.0
      if i > 1
        spacing = sorted_1d[i] - sorted_1d[i - 1]
        # Long spacing gets positive energy, short gets negative
        onsite += onsite_modulation * (spacing > 1.4 ? 1.0 : -1.0)
      end
      if i < n
        spacing = sorted_1d[i + 1] - sorted_1d[i]
        onsite += onsite_modulation * (spacing > 1.4 ? 1.0 : -1.0)
      end
      onsite /= (i == 1 || i == n) ? 1 : 2  # Average for interior sites

      push!(I_idx, i)
      push!(J_idx, i)
      push!(V_vals, onsite)
    end
  end

  # Hopping terms
  for i in 1:(n - 1)
    # Nearest-neighbor hopping
    push!(I_idx, i)
    push!(J_idx, i + 1)
    push!(V_vals, -t)

    push!(I_idx, i + 1)
    push!(J_idx, i)
    push!(V_vals, -t)
  end

  H = sparse(I_idx, J_idx, V_vals, n, n)

  return TightBindingModel(H, sorted_positions)
end

"""
    build_penrose_tight_binding(qc_data::QuasicrystalData{2,T}; 
                                 t::Real=1.0,
                                 cutoff::Real=1.5) where T
Build tight-binding model for 2D Penrose tiling.

# Arguments
- `qc_data`: Penrose quasicrystal data
- `t`: Hopping parameter (default: 1.0)
- `cutoff`: Distance cutoff for neighbors (default: 1.5)

# Returns
- `TightBindingModel`: Tight-binding model with nearest-neighbor hopping
"""
function build_penrose_tight_binding(
  qc_data::QuasicrystalData{2,T}; t::Real=1.0, cutoff::Real=1.5
) where {T}
  n = length(qc_data.positions)

  # Build neighbor list based on distance
  I_idx = Int[]
  J_idx = Int[]
  V_vals = ComplexF64[]

  for i in 1:n
    for j in i:n
      pos_i = qc_data.positions[i]
      pos_j = qc_data.positions[j]
      dist = norm(pos_i - pos_j)

      if i == j
        # On-site energy (zero for now)
        push!(I_idx, i)
        push!(J_idx, i)
        push!(V_vals, 0.0)
      elseif dist < cutoff && dist > 0.1
        # Nearest-neighbor hopping
        push!(I_idx, i)
        push!(J_idx, j)
        push!(V_vals, -t)

        push!(I_idx, j)
        push!(J_idx, i)
        push!(V_vals, -t)
      end
    end
  end

  H = sparse(I_idx, J_idx, V_vals, n, n)

  return TightBindingModel(H, qc_data.positions)
end

"""
    solve_eigenspectrum(model::TightBindingModel; k::Int=min(model.n_sites, 100))
Solve for eigenvalues and eigenvectors of the tight-binding Hamiltonian.

# Arguments
- `model`: Tight-binding model
- `k`: Number of eigenvalues to compute (default: min(n_sites, 100))

# Returns
- `eigenvalues`: Vector of eigenvalues
- `eigenvectors`: Matrix of eigenvectors (columns)
"""
function solve_eigenspectrum(model::TightBindingModel; k::Int=min(model.n_sites, 100))
  if k >= model.n_sites - 1
    # Use dense solver for small systems
    E, V = eigen(Matrix(model.H))
    return real(E), V
  else
    # Use sparse solver for large systems
    E, V = eigs(model.H; nev=k, which=:SM)
    perm = sortperm(real(E))
    return real(E[perm]), V[:, perm]
  end
end

"""
    compute_dos(eigenvalues::Vector{<:Real}; 
                n_bins::Int=100, 
                sigma::Real=0.1)
Compute density of states using Gaussian broadening.

# Arguments
- `eigenvalues`: Energy eigenvalues
- `n_bins`: Number of energy bins (default: 100)
- `sigma`: Gaussian broadening width (default: 0.1)

# Returns
- `energies`: Energy grid points
- `dos`: Density of states at each energy
"""
function compute_dos(eigenvalues::Vector{<:Real}; n_bins::Int=100, sigma::Real=0.1)
  E_min, E_max = extrema(eigenvalues)
  E_range = E_max - E_min
  E_min -= 0.1 * E_range
  E_max += 0.1 * E_range

  energies = range(E_min, E_max; length=n_bins)
  dos = zeros(n_bins)

  # Gaussian broadening
  for E in eigenvalues
    for (i, E_grid) in enumerate(energies)
      dos[i] += exp(-((E - E_grid) / sigma)^2) / (sigma * sqrt(π))
    end
  end

  # Normalize
  dos ./= length(eigenvalues)

  return collect(energies), dos
end

"""
    compute_ipr(eigenvector::Vector{<:Number})
Compute Inverse Participation Ratio (IPR) to measure localization.

IPR = Σᵢ |ψᵢ|⁴ / (Σᵢ |ψᵢ|²)²

For extended states: IPR ~ 1/N (small)
For localized states: IPR ~ 1 (large)

# Arguments
- `eigenvector`: Eigenstate wavefunction

# Returns
- `ipr`: Inverse Participation Ratio
"""
function compute_ipr(eigenvector::Vector{<:Number})
  psi2 = abs2.(eigenvector)
  ipr = sum(psi2 .^ 2) / sum(psi2)^2
  return ipr
end

"""
    compute_all_iprs(eigenvectors::Matrix{<:Number})
Compute IPR for all eigenvectors.

# Arguments
- `eigenvectors`: Matrix of eigenvectors (columns are eigenstates)

# Returns
- `iprs`: Vector of IPR values
"""
function compute_all_iprs(eigenvectors::Matrix{<:Number})
  n_states = size(eigenvectors, 2)
  iprs = zeros(n_states)

  for i in 1:n_states
    iprs[i] = compute_ipr(eigenvectors[:, i])
  end

  return iprs
end

export TightBindingModel
export build_fibonacci_tight_binding, build_penrose_tight_binding
export solve_eigenspectrum, compute_dos
export compute_ipr, compute_all_iprs
