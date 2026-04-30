using Aqua
using QuasiCrystal

@testset "Aqua quality checks" begin
    # Full Aqua sweep, with fine-grained skips below for known
    # cross-package method ambiguities (LatticeCore / NearestNeighbors).
    Aqua.test_all(QuasiCrystal; ambiguities=false)

    # Run ambiguity check restricted to QuasiCrystal itself, so that
    # ambiguities introduced by upstream packages do not fail CI.
    @testset "ambiguities (QuasiCrystal-only)" begin
        Aqua.test_ambiguities(QuasiCrystal)
    end
end
