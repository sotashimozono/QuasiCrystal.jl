@testset "substitution algorithm dispatch (#41)" begin
    # Penrose: DefaultSubstitution and RobinsonTriangleInflation are
    # both fully wired through `inflate_tiles`; DirectTileInflation
    # is an explicit not-implemented shell.
    qc = generate_penrose_substitution(1)
    tiles = qc.tiles
    @test !isempty(tiles)

    inflated_default = QuasiCrystal.inflate_penrose_tiles(tiles, DefaultSubstitution())
    inflated_robin = QuasiCrystal.inflate_penrose_tiles(tiles, RobinsonTriangleInflation())
    @test length(inflated_default) == length(inflated_robin)
    @test length(inflated_robin) > length(tiles)

    # Public single-dispatch API.
    @test length(inflate_tiles(tiles, RobinsonTriangleInflation())) ==
        length(inflated_robin)

    @test_throws ErrorException QuasiCrystal.inflate_penrose_tiles(
        tiles, DirectTileInflation()
    )

    # Ammann–Beenker: DefaultSubstitution and AmmannBeenkerInflation
    # are fully wired; the Penrose-specific Robinson algorithm raises.
    ab = generate_ammann_beenker_substitution(1)
    ab_tiles = ab.tiles
    @test !isempty(ab_tiles)
    inflated_ab_default = QuasiCrystal.inflate_ammann_beenker_tiles(
        ab_tiles, DefaultSubstitution()
    )
    @test length(inflated_ab_default) > length(ab_tiles)
    @test length(inflate_tiles(ab_tiles, AmmannBeenkerInflation())) ==
        length(inflated_ab_default)

    @test_throws ErrorException QuasiCrystal.inflate_ammann_beenker_tiles(
        ab_tiles, RobinsonTriangleInflation()
    )
end
