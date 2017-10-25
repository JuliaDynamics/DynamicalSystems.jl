using Documenter, DynamicalSystems

makedocs(modules=[DynamicalSystems], doctest=false)

deploydocs(
    deps   = Deps.pip("mkdocs==0.16.3",
    "mkdocs-material" ,"python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalSystems.jl.git",
    julia  = "0.6",
    osname = "linux"
)
