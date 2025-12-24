"""
Fibonacci lattice (1D quasicrystal) implementation.
The Fibonacci lattice is the simplest quasicrystal, with two spacings L (long) and S (short)
arranged according to the Fibonacci sequence.
"""

"""
    FibonacciLattice <: AbstractQuasicrystal{1}
One-dimensional Fibonacci quasicrystal lattice.
Spacing ratio is the golden ratio: L/S = φ = (1+√5)/2
"""
struct FibonacciLattice <: AbstractQuasicrystal{1} end

"""
    generate_fibonacci_projection(n_points::Int; method::ProjectionMethod=ProjectionMethod())
Generate Fibonacci lattice using projection from 2D square lattice.
# Arguments
- `n_points::Int`: approximate number of lattice points
- `method::ProjectionMethod`: generation method (default: ProjectionMethod)
# Returns
- `QuasicrystalData{1,Float64}`: generated Fibonacci lattice data
"""
function generate_fibonacci_projection(
    n_points::Int; method::ProjectionMethod=ProjectionMethod()
)
    # Project from 2D square lattice onto a line with irrational slope
    # The slope is 1/φ where φ is the golden ratio

    slope = 1 / ϕ
    acceptance_width = 1.0

    positions = Float64[]
    n_max = ceil(Int, n_points * 1.5)

    for n1 in 0:n_max, n2 in 0:n_max
        # Point in 2D lattice
        point_2d = [float(n1), float(n2)]

        # Project to parallel space (the line)
        # Direction along line: [1, slope]
        direction = [1.0, slope]
        direction = direction / norm(direction)
        pos_par = dot(point_2d, direction)

        # Project to perpendicular space
        perp_direction = [-slope, 1.0]
        perp_direction = perp_direction / norm(perp_direction)
        pos_perp = dot(point_2d, perp_direction)

        # Check acceptance window
        if abs(pos_perp) <= acceptance_width / 2
            push!(positions, pos_par)
        end

        if length(positions) >= n_points
            break
        end
    end

    # Sort positions
    sort!(positions)

    # Trim to requested number
    if length(positions) > n_points
        positions = positions[1:n_points]
    end

    # Convert to vector of vectors for consistency
    positions_vec = [[p] for p in positions]

    tiles = []  # No tiles for 1D lattice

    params = Dict{Symbol,Any}(
        :n_points => length(positions), :slope => slope, :method_name => "projection"
    )

    return QuasicrystalData{1,Float64}(positions_vec, tiles, method, params)
end

"""
    generate_fibonacci_substitution(generations::Int; method::SubstitutionMethod=SubstitutionMethod())
Generate Fibonacci lattice using substitution rules.
Substitution rules: L -> LS, S -> L
# Arguments
- `generations::Int`: number of substitution generations
- `method::SubstitutionMethod`: generation method (default: SubstitutionMethod)
# Returns
- `QuasicrystalData{1,Float64}`: generated Fibonacci lattice data
"""
function generate_fibonacci_substitution(
    generations::Int; method::SubstitutionMethod=SubstitutionMethod()
)
    # Start with L
    sequence = ['L']

    # Apply substitution rules
    for gen in 1:generations
        new_sequence = Char[]
        for symbol in sequence
            if symbol == 'L'
                append!(new_sequence, ['L', 'S'])
            else  # symbol == 'S'
                push!(new_sequence, 'L')
            end
        end
        sequence = new_sequence
    end

    # Convert sequence to positions
    # Let S = 1, L = φ
    L_spacing = ϕ
    S_spacing = 1.0

    positions = [0.0]
    current_pos = 0.0

    for symbol in sequence
        if symbol == 'L'
            current_pos += L_spacing
        else
            current_pos += S_spacing
        end
        push!(positions, current_pos)
    end

    # Convert to vector of vectors
    positions_vec = [[p] for p in positions]

    tiles = []  # No tiles for 1D lattice

    params = Dict{Symbol,Any}(
        :generations => generations,
        :n_points => length(positions),
        :sequence_length => length(sequence),
        :L_spacing => L_spacing,
        :S_spacing => S_spacing,
        :method_name => "substitution",
    )

    return QuasicrystalData{1,Float64}(positions_vec, tiles, method, params)
end

"""
    fibonacci_sequence_length(n::Int)
Calculate the length of Fibonacci sequence after n generations.
F(n) follows the Fibonacci numbers: 1, 2, 3, 5, 8, 13, ...
"""
function fibonacci_sequence_length(n::Int)
    if n <= 0
        return 1
    end
    a, b = 1, 2
    for i in 2:n
        a, b = b, a + b
    end
    return b
end

export FibonacciLattice, generate_fibonacci_projection, generate_fibonacci_substitution
export fibonacci_sequence_length
