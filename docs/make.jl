using DynamicalSystems, TimeseriesPrediction
using Documenter, PyPlot, Literate, DocumenterMarkdown

PyPlot.ioff()
cd(@__DIR__)

# Expand Spatio-temporal examples using Literate:
function replace_includes(str)

    included = ["1Dfield_temporalprediction.jl",
    "2Dfield_crossprediction.jl", "2Dfield_temporalprediction.jl"]

    path = dirname(dirname(pathof(TimeseriesPrediction)))*"/examples/"

    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    return str
end
# Literate it:
Literate.markdown("src/tsprediction/stexamples.jl", "src/tsprediction/";
                  name = "stexamples", preprocess = replace_includes)
# Literate.notebook("src/tsprediction/stexamples.jl", "src/tsprediction/";
#                   name = "stexamples", preprocess = replace_includes)


makedocs(modules=[DynamicalSystems, TimeseriesPrediction],
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
