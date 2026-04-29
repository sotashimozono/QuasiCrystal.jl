using QuasiCrystal, Test
using StaticArrays

@testset "API gap: LatticeCore optional interface" begin
    # Build a small Penrose patch (has bonds and plaquettes). The
    # substitution generator populates `tiles`, which the
    # `plaquettes` machinery needs; the projection generator
    # currently does not.
    pen = generate_penrose_substitution(2)
    build_nearest_neighbor_bonds!(pen; cutoff=1.05)
    @test num_bonds(pen) > 0
    @test num_plaquettes(pen) > 0

    # Also exercise a 1D Fibonacci instance — it has no plaquettes
    # but should still respond to vertex/bond queries.
    fib = generate_fibonacci_substitution(5)
    build_nearest_neighbor_bonds!(fib; cutoff=GOLDEN_RATIO + 0.05)
    @test num_bonds(fib) > 0

    @testset "neighbor_bonds(lat, i) (issue #31)" begin
        # Pick a site that has at least one neighbour bond.
        i = findfirst(!isempty, get_nearest_neighbors(pen))
        @test i !== nothing
        bs = collect(neighbor_bonds(pen, i))
        @test !isempty(bs)
        # Every returned bond must touch site i and must come from
        # the actual bond list (preserving the stored `:nearest` tag
        # and unwrapped vector — not a freshly fabricated bond).
        for b in bs
            @test b.i == i || b.j == i
            @test b.type === :nearest
            @test b in get_bonds(pen)
        end
        # And the count must match the adjacency-list degree.
        @test length(bs) == length(get_nearest_neighbors(pen)[i])

        # 1D sanity check on Fibonacci.
        i1 = findfirst(!isempty, get_nearest_neighbors(fib))
        @test i1 !== nothing
        @test all(b -> b.i == i1 || b.j == i1, neighbor_bonds(fib, i1))
    end

    @testset "element_positions iterators (issue #32)" begin
        # VertexCenter: identical to `positions(data)`.
        vps = collect(element_positions(pen, VertexCenter()))
        @test vps == get_positions(pen)
        @test length(vps) == num_sites(pen)

        # BondCenter: midpoint of every stored bond.
        bps = collect(element_positions(pen, BondCenter()))
        @test length(bps) == num_bonds(pen)
        for (k, b) in enumerate(get_bonds(pen))
            @test bps[k] ≈ bond_center(pen, b)
        end

        # PlaquetteCenter: each Plaquette's `center` field.
        pps = collect(element_positions(pen, PlaquetteCenter()))
        @test length(pps) == num_plaquettes(pen)
        for (k, p) in enumerate(plaquettes(pen))
            @test pps[k] == p.center
        end
    end

    @testset "element_neighbors PlaquetteCenter (issue #33)" begin
        ps = plaquettes(pen)
        # Find a plaquette index that actually has a neighbour
        # (Penrose tilings always have some plaquette with a shared
        # edge, so this should succeed on a non-degenerate patch).
        ks = [k for k in 1:length(ps) if !isempty(element_neighbors(pen, PlaquetteCenter(), k))]
        @test !isempty(ks)
        k = first(ks)
        nbrs = element_neighbors(pen, PlaquetteCenter(), k)
        # Every reported neighbour must share at least one boundary
        # edge (vertex pair) with plaquette k.
        function _edges(p)
            vs = p.vertices
            n = length(vs)
            return Set(Tuple{Int,Int}[
                (min(vs[i], vs[mod1(i + 1, n)]), max(vs[i], vs[mod1(i + 1, n)])) for i in 1:n
            ])
        end
        e_k = _edges(ps[k])
        for j in nbrs
            @test j != k
            @test !isempty(intersect(e_k, _edges(ps[j])))
        end
    end

    @testset "positions(lat) returns the backing vector (issue #34)" begin
        # Specialisation should hand back the actual storage.
        @test positions(pen) === pen.positions
        @test positions(fib) === fib.positions
        # And it must agree with `get_positions`.
        @test positions(pen) == get_positions(pen)
    end

    @testset "topology trait :quasiperiodic (issue #35)" begin
        @test topology(pen) === TopologyTrait{:quasiperiodic}()
        @test topology(fib) === TopologyTrait{:quasiperiodic}()
        ab = generate_ammann_beenker_substitution(2)
        @test topology(ab) === TopologyTrait{:quasiperiodic}()
    end

    @testset "incident query (issue #36)" begin
        # Pick a vertex with at least one bond.
        v = findfirst(!isempty, get_nearest_neighbors(pen))
        @test v !== nothing

        # Vertex ↔ Bond round-trip.
        v_bonds = incident(pen, VertexCenter(), BondCenter(), v)
        @test !isempty(v_bonds)
        for kb in v_bonds
            ends = incident(pen, BondCenter(), VertexCenter(), kb)
            @test v in ends
            @test length(ends) == 2
        end

        # Vertex ↔ Plaquette round-trip: pick a vertex that lies on a tile.
        v_with_p = findfirst(i -> any(p -> i in p.vertices, plaquettes(pen)), 1:num_sites(pen))
        @test v_with_p !== nothing
        v_plaqs = incident(pen, VertexCenter(), PlaquetteCenter(), v_with_p)
        @test !isempty(v_plaqs)
        for kp in v_plaqs
            verts = incident(pen, PlaquetteCenter(), VertexCenter(), kp)
            @test v_with_p in verts
        end

        # Bond ↔ Plaquette: pick a plaquette and verify the boundary
        # edges resolve back to bonds (when those bonds exist).
        ps = plaquettes(pen)
        for kp in 1:length(ps)
            p_bonds = incident(pen, PlaquetteCenter(), BondCenter(), kp)
            for kb in p_bonds
                p_back = incident(pen, BondCenter(), PlaquetteCenter(), kb)
                @test kp in p_back
            end
        end
    end
end
