"""
Fibonacci lattice (1D quasicrystal): the simplest cut-and-project
quasicrystal, with two spacings `L` and `S` in golden ratio.
"""

"""
    FibonacciLattice <: AbstractQuasicrystal{1}

Topology marker for a 1D Fibonacci chain. Pass to
[`build_quasicrystal`](@ref) or call
[`generate_fibonacci_projection`](@ref) /
[`generate_fibonacci_substitution`](@ref) directly.
"""
struct FibonacciLattice <: AbstractQuasicrystal{1} end

"""
    generate_fibonacci_projection(n_points::Int;
                                  method::ProjectionMethod = ProjectionMethod())
        → QuasicrystalData{1, Float64}

Generate a Fibonacci chain by projecting the integer lattice `Z^2`
onto a line of slope `1/ϕ`, accepting points whose perpendicular
projection falls inside a unit window. Returns a 1D
[`QuasicrystalData`](@ref).
"""
function generate_fibonacci_projection(
    n_points::Int; method::ProjectionMethod=ProjectionMethod()
)
    slope = 1 / ϕ
    acceptance_width = 1.0

    raw_positions = Float64[]
    n_max = ceil(Int, n_points * 1.5)

    direction = SVector(1.0, slope)
    direction_unit = direction / norm(direction)
    perp_direction = SVector(-slope, 1.0)
    perp_unit = perp_direction / norm(perp_direction)

    for n1 in 0:n_max, n2 in 0:n_max
        point_2d = SVector(float(n1), float(n2))
        pos_par = dot(point_2d, direction_unit)
        pos_perp = dot(point_2d, perp_unit)

        if abs(pos_perp) <= acceptance_width / 2
            push!(raw_positions, pos_par)
        end
        length(raw_positions) >= n_points && break
    end

    sort!(raw_positions)
    if length(raw_positions) > n_points
        raw_positions = raw_positions[1:n_points]
    end

    positions = [SVector{1,Float64}(p) for p in raw_positions]
    tiles = Tile{1,Float64}[]

    params = Dict{Symbol,Any}(
        :n_points => length(positions),
        :slope => slope,
        :method_name => "projection",
        :window_shape => :interval_1d,
    )
    return QuasicrystalData{1,Float64}(FibonacciLattice(), positions, tiles, method, params)
end

"""
    generate_fibonacci_substitution(generations::Int;
                                    method::SubstitutionMethod = SubstitutionMethod())
        → QuasicrystalData{1, Float64}

Generate a Fibonacci chain by iterating the substitution rules
`L → LS`, `S → L` from the axiom `L`. The resulting positions are
a cumulative sum of `L`-spacings (`= ϕ`) and `S`-spacings (`= 1`).
"""
function generate_fibonacci_substitution(
    generations::Int; method::SubstitutionMethod=SubstitutionMethod()
)
    sequence = ['L']
    for _ in 1:generations
        new_sequence = Char[]
        for symbol in sequence
            if symbol == 'L'
                append!(new_sequence, ['L', 'S'])
            else
                push!(new_sequence, 'L')
            end
        end
        sequence = new_sequence
    end

    L_spacing = ϕ
    S_spacing = 1.0

    raw_positions = [0.0]
    current_pos = 0.0
    for symbol in sequence
        current_pos += symbol == 'L' ? L_spacing : S_spacing
        push!(raw_positions, current_pos)
    end

    positions = [SVector{1,Float64}(p) for p in raw_positions]
    tiles = Tile{1,Float64}[]

    params = Dict{Symbol,Any}(
        :generations => generations,
        :n_points => length(positions),
        :sequence_length => length(sequence),
        :L_spacing => L_spacing,
        :S_spacing => S_spacing,
        :method_name => "substitution",
        :window_shape => :none,
    )
    return QuasicrystalData{1,Float64}(FibonacciLattice(), positions, tiles, method, params)
end

"""
    fibonacci_sequence_length(n::Int) → Int

Length of the Fibonacci substitution sequence after `n` generations,
starting from the axiom `L`: `1, 2, 3, 5, 8, 13, ...`.
"""
function fibonacci_sequence_length(n::Int)
    n <= 0 && return 1
    a, b = 1, 2
    for _ in 2:n
        a, b = b, a + b
    end
    return b
end
