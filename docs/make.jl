using QuasiCrystal
using Documenter

makedocs(; sitename="QuasiCrystal.jl", modules=[QuasiCrystal], pages=["Home" => "index.md"])

deploydocs(; repo="github.com/sotashimozono/QuasiCrystal.jl.git")
