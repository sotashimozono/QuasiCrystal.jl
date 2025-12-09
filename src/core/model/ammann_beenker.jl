"""
Ammann-Beenker tiling (octagonal quasicrystal) implementation.
The Ammann-Beenker tiling has 8-fold rotational symmetry.
"""

"""
    AmmannBeenker <: AbstractQuasicrystal{2}
Ammann-Beenker tiling with 8-fold rotational symmetry.
This tiling uses squares and 45-degree rhombi.
"""
struct AmmannBeenker <: AbstractQuasicrystal{2} end

const SQRT2 = sqrt(2.0)

"""
    generate_ammann_beenker_projection(radius::Real; method::ProjectionMethod=ProjectionMethod())
Generate Ammann-Beenker tiling using projection from 4D hypercubic lattice.
# Arguments
- `radius::Real`: approximate radius of the generated pattern
- `method::ProjectionMethod`: generation method (default: ProjectionMethod)
# Returns
- `QuasicrystalData{2,Float64}`: generated Ammann-Beenker tiling data
"""
function generate_ammann_beenker_projection(
  radius::Real; method::ProjectionMethod=ProjectionMethod()
)
  # Define projection from 4D to 2D for octagonal symmetry
  theta = π / 4  # 45 degrees

  # Projection matrix to parallel space (physical 2D)
  E_par = zeros(4, 2)
  for i in 1:4
    E_par[i, 1] = cos((i - 1) * theta)
    E_par[i, 2] = sin((i - 1) * theta)
  end

  # Projection matrix to perpendicular space (window in 2D)
  E_perp = zeros(4, 2)
  for i in 1:4
    E_perp[i, 1] = cos((i - 1) * theta + π / 4)
    E_perp[i, 2] = sin((i - 1) * theta + π / 4)
  end

  # Acceptance window size
  window_size = 0.5

  # Generate 4D lattice points and project
  positions = Vector{Float64}[]
  n_max = ceil(Int, radius * 1.5)

  for n1 in (-n_max):n_max, n2 in (-n_max):n_max, n3 in (-n_max):n_max, n4 in (-n_max):n_max
    lattice_point = [float(n1), float(n2), float(n3), float(n4)]

    # Project to parallel space (physical 2D)
    pos_par = E_par' * lattice_point

    # Check if within radius
    if norm(pos_par) > radius
      continue
    end

    # Project to perpendicular space
    pos_perp = E_perp' * lattice_point

    # Check if within acceptance window (square)
    if all(abs.(pos_perp) .<= window_size)
      push!(positions, pos_par)
    end
  end

  tiles = Tile{2,Float64}[]  # Tiles to be constructed

  params = Dict{Symbol,Any}(
    :radius => radius,
    :n_max => n_max,
    :window_size => window_size,
    :n_vertices => length(positions),
    :symmetry => 8,
  )

  return QuasicrystalData{2,Float64}(positions, tiles, method, params)
end

"""
    generate_ammann_beenker_substitution(generations::Int; method::SubstitutionMethod=SubstitutionMethod())
Generate Ammann-Beenker tiling using substitution rules.
# Arguments
- `generations::Int`: number of inflation steps
- `method::SubstitutionMethod`: generation method (default: SubstitutionMethod)
# Returns
- `QuasicrystalData{2,Float64}`: generated Ammann-Beenker tiling data
"""
function generate_ammann_beenker_substitution(
  generations::Int; method::SubstitutionMethod=SubstitutionMethod()
)
  # Initial tiles: start with squares and rhombi
  initial_tiles = Tile{2,Float64}[]

  # Create initial square
  square_vertices = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
  square_center = [0.5, 0.5]
  push!(initial_tiles, Tile{2,Float64}(square_vertices, 1, square_center))  # type 1 = square

  # Create initial rhombi in octagonal pattern
  for i in 0:7
    angle = i * π / 4
    # Add 45-degree rhombus
    v1 = [cos(angle), sin(angle)]
    v2 = v1 + [cos(angle + π / 4), sin(angle + π / 4)]
    v3 = v2 + [cos(angle + π), sin(angle + π)]
    v4 = v1 + [cos(angle + π), sin(angle + π)]

    center = (v1 + v2 + v3 + v4) / 4
    push!(initial_tiles, Tile{2,Float64}([v1, v2, v3, v4], 2, center))  # type 2 = rhombus
  end

  tiles = initial_tiles

  # Apply substitution rules
  for gen in 1:generations
    tiles = inflate_ammann_beenker_tiles(tiles)
  end

  # Extract unique positions using a Set for O(n) complexity
  position_set = Set{Vector{Float64}}()
  for tile in tiles
    for v in tile.vertices
      push!(position_set, v)
    end
  end
  positions = collect(position_set)

  params = Dict{Symbol,Any}(
    :generations => generations,
    :n_tiles => length(tiles),
    :n_vertices => length(positions),
    :symmetry => 8,
  )

  return QuasicrystalData{2,Float64}(positions, tiles, method, params)
end

"""
    inflate_ammann_beenker_tiles(tiles::Vector{Tile{2,Float64}})
Apply Ammann-Beenker substitution rules.
"""
function inflate_ammann_beenker_tiles(tiles::Vector{Tile{2,Float64}})
  new_tiles = Tile{2,Float64}[]
  inflation_factor = 1 + SQRT2

  for tile in tiles
    if tile.type == 1  # Square
      # Inflate square
      append!(new_tiles, inflate_ab_square(tile, inflation_factor))
    else  # Rhombus
      # Inflate rhombus
      append!(new_tiles, inflate_ab_rhombus(tile, inflation_factor))
    end
  end

  return new_tiles
end

"""
    inflate_ab_square(tile::Tile{2,Float64}, factor::Float64)
Inflate an Ammann-Beenker square.
"""
function inflate_ab_square(tile::Tile{2,Float64}, factor::Float64)
  # Simplified: scale the square
  v = tile.vertices
  scaled_tile = Tile{2,Float64}([factor * vertex for vertex in v], 1, factor * tile.center)
  return [scaled_tile]
end

"""
    inflate_ab_rhombus(tile::Tile{2,Float64}, factor::Float64)
Inflate an Ammann-Beenker rhombus.
"""
function inflate_ab_rhombus(tile::Tile{2,Float64}, factor::Float64)
  # Simplified: scale the rhombus
  v = tile.vertices
  scaled_tile = Tile{2,Float64}([factor * vertex for vertex in v], 2, factor * tile.center)
  return [scaled_tile]
end

export AmmannBeenker,
  generate_ammann_beenker_projection, generate_ammann_beenker_substitution
