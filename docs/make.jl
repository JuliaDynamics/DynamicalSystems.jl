using Documenter, DynamicalSystems, PyPlot

PyPlot.ioff()

makedocs(modules=[DynamicalSystems], doctest=false)

deploydocs(
    deps   = Deps.pip("Tornado>=4.0.0,<5.0.0", "mkdocs",
    "mkdocs-material" ,"python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalSystems.jl.git",
    julia  = "0.7",
    osname = "linux"
)
