module QuasiCrystal
using LinearAlgebra, Plots
include("core/abstractlattice.jl")
include("core/abstractquasicrystals.jl")
include("core/topology.jl")
include("core/lattice.jl")
include("core/interface.jl")
include("core/model/fibonacci.jl")
include("core/model/penrose.jl")
include("core/model/ammann_beenker.jl")
include("utils/visualization.jl")
end
