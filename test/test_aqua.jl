using Aqua
using QuasiCrystal

@testset "Aqua quality checks" begin
    # Full Aqua sweep. Plots was previously ignored in `stale_deps`
    # because the dispatch only fired through the top-level `using
    # Plots`; it is now a `[weakdeps]` trigger for
    # `QuasiCrystalPlotsExt` so the default Aqua check passes.
    Aqua.test_all(QuasiCrystal; ambiguities=false)

    # Run ambiguity check restricted to QuasiCrystal itself, so that
    # ambiguities introduced by upstream packages do not fail CI.
    @testset "ambiguities (QuasiCrystal-only)" begin
        Aqua.test_ambiguities(QuasiCrystal)
    end
end
