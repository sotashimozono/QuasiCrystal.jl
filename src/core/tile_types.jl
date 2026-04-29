"""
Semantic tile-type taxonomy for quasicrystalline tilings.

Replaces the legacy `Tile.type::Int` (1/2) encoding with named
singleton subtypes of [`TileType`](@ref). Each generator family has
its own concrete tile types — Penrose P3 uses [`FatRhombus`](@ref) /
[`ThinRhombus`](@ref); Ammann–Beenker uses [`Square`](@ref) /
[`RhombusAB`](@ref).

Use [`tile_type_symbol`](@ref) to obtain the LatticeCore plaquette
tag, e.g. `:fat_rhombus`, `:thin_rhombus`, `:square`, `:rhombus`.
"""

"""
    abstract type TileType end

Root of the semantic tile-type hierarchy. Concrete subtypes are
singleton structs (e.g. [`FatRhombus`](@ref)) that replace the legacy
integer `Tile.type` field.
"""
abstract type TileType end

# --- Penrose P3 ---------------------------------------------------

"""Fat rhombus (Penrose P3, 72°/108° interior angles)."""
struct FatRhombus <: TileType end

"""Thin rhombus (Penrose P3, 36°/144° interior angles)."""
struct ThinRhombus <: TileType end

# --- Ammann–Beenker -----------------------------------------------

"""Square tile of the Ammann–Beenker tiling."""
struct Square <: TileType end

"""45°/135° rhombus tile of the Ammann–Beenker tiling."""
struct RhombusAB <: TileType end

# --- Symbol mapping ------------------------------------------------

"""
    tile_type_symbol(t::TileType) → Symbol

Return the canonical plaquette tag (`Symbol`) for a tile type, used
when promoting a [`Tile`](@ref) into a `LatticeCore.Plaquette`.
"""
tile_type_symbol(::FatRhombus) = :fat_rhombus
tile_type_symbol(::ThinRhombus) = :thin_rhombus
tile_type_symbol(::Square) = :square
tile_type_symbol(::RhombusAB) = :rhombus
