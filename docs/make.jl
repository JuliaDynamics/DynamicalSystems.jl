using DynamicalSystems, DelayEmbeddings, ChaosTools, DynamicalSystemsBase
using RecurrenceAnalysis, Documenter, PyPlot, DocumenterMarkdown
using InteractiveChaos

PyPlot.ioff()
cd(@__DIR__)

makedocs(modules=[DynamicalSystems, ChaosTools, DynamicalSystemsBase,
DelayEmbeddings, RecurrenceAnalysis, InteractiveChaos],
doctest=false, root = @__DIR__, format = Markdown())

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

PyPlot.ion()
