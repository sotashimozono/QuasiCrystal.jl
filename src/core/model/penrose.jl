"""
Penrose P3 (rhombus) tiling implementation.
Penrose tilings are aperiodic tilings discovered by Roger Penrose with 5-fold rotational symmetry.
"""

"""
    PenroseP3 <: AbstractQuasicrystal{2}
Penrose P3 rhombus tiling with 5-fold rotational symmetry.
This tiling uses two types of rhombi: fat (72°) and thin (36°).
"""
struct PenroseP3 <: AbstractQuasicrystal{2} end

"""
    generate_penrose_projection(radius::Real; method::ProjectionMethod=ProjectionMethod())
Generate Penrose P3 tiling using the projection method from a 5D hypercubic lattice.
# Arguments
- `radius::Real`: approximate radius of the generated pattern
- `method::ProjectionMethod`: generation method (default: ProjectionMethod)
# Returns
- `QuasicrystalData{2,Float64}`: generated Penrose tiling data
"""
function generate_penrose_projection(
  radius::Real; method::ProjectionMethod=ProjectionMethod()
)
  # Define basis vectors for projection from 5D to 2D
  # These define the acceptance window in 5D space
  theta = 2π / 5

  # Projection matrices
  # E_par: parallel space (physical 2D)
  E_par = zeros(5, 2)
  for i in 1:5
    E_par[i, 1] = cos((i - 1) * theta)
    E_par[i, 2] = sin((i - 1) * theta)
  end

  # E_perp: perpendicular space (acceptance window in 3D)
  E_perp = zeros(5, 3)
  for i in 1:5
    E_perp[i, 1] = cos(2 * (i - 1) * theta)
    E_perp[i, 2] = sin(2 * (i - 1) * theta)
    E_perp[i, 3] = cos(3 * (i - 1) * theta)
  end

  # Acceptance window size
  window_size = 0.5

  # Generate 5D lattice points and project
  positions = Vector{Float64}[]
  n_max = ceil(Int, radius * 1.5)

  for n1 in (-n_max):n_max,
    n2 in (-n_max):n_max,
    n3 in (-n_max):n_max,
    n4 in (-n_max):n_max,
    n5 in (-n_max):n_max

    lattice_point = [float(n1), float(n2), float(n3), float(n4), float(n5)]

    # Project to parallel space (physical 2D)
    pos_par = E_par' * lattice_point

    # Check if within radius
    if norm(pos_par) > radius
      continue
    end

    # Project to perpendicular space
    pos_perp = E_perp' * lattice_point

    # Check if within acceptance window (hypercube)
    if all(abs.(pos_perp) .<= window_size)
      push!(positions, pos_par)
    end
  end

  # Generate tiles from the Delaunay triangulation or direct construction
  # For now, return positions
  tiles = Tile{2,Float64}[]  # To be implemented: tile construction

  params = Dict{Symbol,Any}(
    :radius => radius,
    :n_max => n_max,
    :window_size => window_size,
    :n_vertices => length(positions),
  )

  return QuasicrystalData{2,Float64}(positions, tiles, method, params)
end

"""
    generate_penrose_substitution(generations::Int; method::SubstitutionMethod=SubstitutionMethod())
Generate Penrose P3 tiling using substitution (inflation) rules.
# Arguments
- `generations::Int`: number of inflation steps
- `method::SubstitutionMethod`: generation method (default: SubstitutionMethod)
# Returns
- `QuasicrystalData{2,Float64}`: generated Penrose tiling data
"""
function generate_penrose_substitution(
  generations::Int; method::SubstitutionMethod=SubstitutionMethod()
)
  # Initial tiles: start with a set of fat and thin rhombi
  # Fat rhombus: 72° angle
  # Thin rhombus: 36° angle

  # Define initial fat rhombus vertices
  angle_fat = deg2rad(72)
  angle_thin = deg2rad(36)

  # Start with a simple configuration
  initial_tiles = Tile{2,Float64}[]

  # Create 5 fat rhombi in a star pattern around origin
  for i in 0:4
    angle = i * 2π / 5
    v1 = [0.0, 0.0]
    v2 = [cos(angle), sin(angle)]
    v3 = v2 + [cos(angle + angle_fat), sin(angle + angle_fat)]
    v4 = [cos(angle + angle_fat), sin(angle + angle_fat)]

    center = (v1 + v2 + v3 + v4) / 4
    push!(initial_tiles, Tile{2,Float64}([v1, v2, v3, v4], 1, center))  # type 1 = fat
  end

  tiles = initial_tiles

  # Apply substitution rules 'generations' times
  for gen in 1:generations
    tiles = inflate_penrose_tiles(tiles)
  end

  # Extract all unique positions using a Set for O(n) complexity
  position_set = Set{Vector{Float64}}()
  for tile in tiles
    for v in tile.vertices
      push!(position_set, v)
    end
  end
  positions = collect(position_set)

  params = Dict{Symbol,Any}(
    :generations => generations, :n_tiles => length(tiles), :n_vertices => length(positions)
  )

  return QuasicrystalData{2,Float64}(positions, tiles, method, params)
end

"""
    inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}})
Apply Penrose substitution rules to inflate tiles.
"""
function inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}})
  new_tiles = Tile{2,Float64}[]

  for tile in tiles
    # Apply inflation based on tile type
    if tile.type == 1  # Fat rhombus
      # Fat rhombus inflates to 1 fat + 2 thin
      append!(new_tiles, inflate_fat_rhombus(tile))
    else  # Thin rhombus
      # Thin rhombus inflates to 1 fat
      append!(new_tiles, inflate_thin_rhombus(tile))
    end
  end

  return new_tiles
end

"""
    inflate_fat_rhombus(tile::Tile{2,Float64})
Inflate a fat rhombus into smaller tiles.

NOTE: This is a simplified placeholder implementation that scales the tile
rather than applying proper Penrose inflation rules. A complete implementation
would subdivide the fat rhombus into 1 fat and 2 thin rhombi according to
the Penrose matching rules.
"""
function inflate_fat_rhombus(tile::Tile{2,Float64})
  # TODO: Implement proper Penrose inflation rules
  # Proper rule: Fat rhombus → 1 fat + 2 thin rhombi
  # Current: Simplified scaling for demonstration
  v = tile.vertices

  # Scale and create new tiles
  scaled_tile = Tile{2,Float64}([ϕ * vertex for vertex in v], 1, ϕ * tile.center)

  return [scaled_tile]
end

"""
    inflate_thin_rhombus(tile::Tile{2,Float64})
Inflate a thin rhombus into smaller tiles.

NOTE: This is a simplified placeholder implementation. A complete implementation
would subdivide the thin rhombus into 1 fat rhombus according to Penrose rules.
"""
function inflate_thin_rhombus(tile::Tile{2,Float64})
  # TODO: Implement proper Penrose inflation rules
  # Proper rule: Thin rhombus → 1 fat rhombus
  # Current: Simplified scaling for demonstration
  v = tile.vertices

  scaled_tile = Tile{2,Float64}([ϕ * vertex for vertex in v], 2, ϕ * tile.center)

  return [scaled_tile]
end

export PenroseP3, generate_penrose_projection, generate_penrose_substitution
