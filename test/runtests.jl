ENV["GKSwstype"] = "100"

using QuasiCrystal, Test
using Plots, LinearAlgebra
const dirs = ["model"]

@testset "tests" begin
    test_args = copy(ARGS)
    println("Passed arguments ARGS = $(test_args) to tests.")
    @time for dir in dirs
        dirpath = joinpath(@__DIR__, dir)
        println("\nTest $(dirpath)")
        files = sort(
            filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(dirpath))
        )
        if isempty(files)
            println("  No test files found in $(dirpath).")
        else
            for f in files
                filepath = joinpath(dirpath, f)
                @time begin
                    println("  Including $(filepath)")
                    include(filepath)
                end
            end
        end
    end
end
