module QuasiCrystal

using LinearAlgebra
using SparseArrays
using StaticArrays
using Plots

# ---- LatticeCore imports --------------------------------------------
#
# Everything that was previously duplicated inside QuasiCrystal
# (`AbstractLattice{D}`, `Bond`, `Connection`, `UnitCell`, etc.) now
# lives in LatticeCore. We import the symbols we specialise and
# re-export the rest so `using QuasiCrystal` continues to work for
# downstream code.

using LatticeCore

import LatticeCore:
    AbstractLattice,
    num_sites,
    position,
    positions,
    neighbors,
    boundary,
    bonds,
    neighbor_bonds,
    site_layout,
    site_type,
    sublattice,
    num_sublattices,
    size_trait,
    is_bipartite,
    periodicity,
    reciprocal_support,
    topology,
    basis_vectors,
    reciprocal_lattice,
    fourier_module,
    momentum_lattice,
    to_real,
    to_lattice,
    to_hyper,
    Bond,
    bond_center,
    LatticeBoundary,
    AbstractAxisBC,
    PeriodicAxis,
    OpenAxis,
    TwistedAxis,
    NoModifier,
    SSD,
    apply_axis_bc,
    axis_phase,
    bond_weight,
    AbstractCoordinate,
    RealSpace,
    LatticeCoord,
    HigherDimCoord,
    AbstractIndexing,
    RowMajor,
    ColMajor,
    Snake,
    site_index,
    lattice_coord,
    AbstractSiteType,
    IsingSite,
    PottsSite,
    XYSite,
    HeisenbergSite,
    EmptySite,
    element_type,
    state_type,
    random_state,
    zero_state,
    domain,
    AbstractSiteLayout,
    UniformLayout,
    SublatticeLayout,
    ExplicitLayout,
    TopologyTrait,
    Periodic,
    Aperiodic,
    AbstractReciprocalSupport,
    HasReciprocal,
    HasFourierModule,
    NoReciprocal,
    AbstractSizeTrait,
    FiniteSize,
    InfiniteSize,
    QuasiInfiniteSize,
    is_finite,
    materialize,
    require_finite,
    AbstractMomentumLattice,
    PeriodicMomentumLattice,
    num_k_points,
    k_point,
    reciprocal_basis,
    monkhorst_pack,
    gamma_centered,
    structure_factor,
    AcceptanceWindow,
    HyperReciprocalLattice,
    BraggPeakSet,
    AbstractLatticeElement,
    VertexCenter,
    BondCenter,
    PlaquetteCenter,
    CellCenter,
    AbstractBoundaryCondition,
    AbstractBoundaryModifier

# ---- QuasiCrystal source files ---------------------------------------

include("core/abstractquasicrystals.jl")
include("core/interface.jl")
include("core/model/fibonacci.jl")
include("core/model/penrose.jl")
include("core/model/ammann_beenker.jl")
include("core/fourier/window.jl")
include("core/fourier/fourier.jl")
include("utils/visualization.jl")

# ---- Exports ---------------------------------------------------------

# QuasiCrystal-local types
export AbstractQuasicrystal, AbstractGenerationMethod
export ProjectionMethod, SubstitutionMethod
export QuasicrystalData, Tile
export GOLDEN_RATIO, ϕ
export FibonacciLattice, PenroseP3, AmmannBeenker
export generate_fibonacci_projection, generate_fibonacci_substitution
export generate_penrose_projection, generate_penrose_substitution
export generate_ammann_beenker_projection, generate_ammann_beenker_substitution
export fibonacci_sequence_length
export build_quasicrystal
export get_positions, get_bonds, get_nearest_neighbors, num_bonds
export build_nearest_neighbor_bonds!
export visualize_quasicrystal_positions
# Fourier analysis
export IntervalWindow, window_fourier
export hyper_reciprocal_lattice, bragg_peaks

# Re-exports from LatticeCore
export AbstractLattice
export Bond, bond_center
export LatticeBoundary, AbstractAxisBC, PeriodicAxis, OpenAxis, TwistedAxis
export AbstractBoundaryCondition
export AbstractBoundaryModifier, NoModifier, SSD
export AbstractIndexing, RowMajor, ColMajor, Snake
export AbstractSiteLayout, UniformLayout, SublatticeLayout, ExplicitLayout
export AbstractSiteType, IsingSite, PottsSite, XYSite, HeisenbergSite, EmptySite
export AbstractCoordinate, RealSpace, LatticeCoord, HigherDimCoord
export AbstractLatticeElement, VertexCenter, BondCenter, PlaquetteCenter, CellCenter
export num_sites, position, positions, neighbors, boundary, bonds, neighbor_bonds
export site_layout, site_type, sublattice, num_sublattices
export site_index, lattice_coord
export is_bipartite, topology
export TopologyTrait
export Periodic, Aperiodic, periodicity
export AbstractReciprocalSupport, HasReciprocal, HasFourierModule, NoReciprocal
export reciprocal_support
export AbstractSizeTrait, FiniteSize, InfiniteSize, QuasiInfiniteSize
export size_trait, is_finite
export basis_vectors, reciprocal_lattice, fourier_module, momentum_lattice
export to_real, to_lattice, to_hyper
export monkhorst_pack, gamma_centered, structure_factor
export AcceptanceWindow, HyperReciprocalLattice, BraggPeakSet
export AbstractMomentumLattice, PeriodicMomentumLattice
export num_k_points, k_point, reciprocal_basis
export apply_axis_bc, axis_phase, bond_weight
export element_type, state_type, random_state, zero_state, domain
export materialize, require_finite

end # module QuasiCrystal
