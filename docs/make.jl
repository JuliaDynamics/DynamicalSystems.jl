using DynamicalSystems, DelayEmbeddings, ChaosTools, DynamicalSystemsBase
using Documenter, PyPlot, Literate, DocumenterMarkdown

PyPlot.ioff()
cd(@__DIR__)

makedocs(modules=[DynamicalSystems, ChaosTools, DynamicalSystemsBase, DelayEmbeddings],
doctest=false, root = @__DIR__, format = :markdown)

close("all")

if !Sys.iswindows()
    deploydocs(
        deps   = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4",
        "python-markdown-math", "pygments", "pymdown-extensions"),
        repo   = "github.com/JuliaDynamics/DynamicalSystems.jl.git",
        target = "site",
        make = () -> run(`mkdocs build`)
    )
end
