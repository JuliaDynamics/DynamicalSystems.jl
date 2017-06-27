using Documenter, DynamicalSystems

makedocs(modules=[DynamicalSystems], doctest=false)

deploydocs(
    deps   = Deps.pip("mkdocs", "mkdocs-material" ,"python-markdown-math", "pygments"),
    repo   = "github.com/Datseris/DynamicalSystems.jl.git",
    julia  = "0.6",
    osname = "linux"
)
