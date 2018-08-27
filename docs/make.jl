using Documenter, DynamicalSystems, PyPlot
PyPlot.ioff()

pkg"add DiffEqCallbacks LinearAlgebra DynamicalSystemsBase ChaosTools"

# Add TimeseriesPrediction for documentation:
using Pkg
pkg"add TimeseriesPrediction#0.7-inplacerec"
using TimeseriesPrediction


makedocs(modules=[DynamicalSystems], doctest=false)

deploydocs(
    deps   = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4",
    "python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalSystems.jl.git",
    julia  = "0.7",
    osname = "linux"
)
