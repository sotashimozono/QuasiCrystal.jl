using QuasiCrystal
using Documenter
using Downloads

assets_dir = joinpath(@__DIR__, "src", "assets")
mkpath(assets_dir)
favicon_path = joinpath(assets_dir, "favicon.ico")

Downloads.download("https://github.com/sotashimozono.png", favicon_path)

# Regenerate all visual-verification figures from the live source
# tree so every doc build re-checks the Fourier pipeline end-to-end.
include(joinpath(@__DIR__, "generate_figures.jl"))

makedocs(;
    sitename="QuasiCrystal.jl",
    modules=[QuasiCrystal],
    format=Documenter.HTML(;
        canonical="https://codes.sota-shimozono.com/QuasiCrystal.jl/stable/",
        prettyurls=get(ENV, "CI", "false") == "true",
        mathengine=MathJax3(
            Dict(
                :tex => Dict(
                    :inlineMath => [["\$", "\$"], ["\\(", "\\)"]],
                    :tags => "ams",
                    :packages => ["base", "ams", "autoload", "physics"],
                ),
            ),
        ),
        assets=["assets/favicon.ico"],
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => [
            "guide/cut_and_project.md",
            "guide/fourier_analysis.md",
        ],
        "Gallery (visual verification)" => [
            "gallery/fibonacci.md",
            "gallery/ammann_beenker.md",
            "gallery/penrose.md",
        ],
        "API" => "api.md",
    ],
    warnonly=[:missing_docs, :cross_references],
)

deploydocs(; repo="github.com/sotashimozono/QuasiCrystal.jl.git", devbranch="master")
