using Documenter, DynamicalSystems, PyPlot
PyPlot.ioff()

using TimeseriesPrediction

makedocs(modules=[DynamicalSystems], doctest=false)

deploydocs(
    deps   = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4",
    "python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalSystems.jl.git",
    julia  = "1.0",
    osname = "linux"
)
