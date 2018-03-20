using Documenter, DynamicalSystems

makedocs(modules=[DynamicalSystems,TimeseriesPrediction], doctest=false)

deploydocs(
    deps   = Deps.pip("Tornado>=4.0.0,<5.0.0", "mkdocs",
    "mkdocs-material" ,"python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalSystems.jl.git",
    julia  = "0.6",
    osname = "linux"
)
