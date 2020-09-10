cd(@__DIR__)
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using Pkg
Pkg.activate(@__DIR__)
CI && Pkg.instantiate()

using DynamicalSystems, DelayEmbeddings, ChaosTools, DynamicalSystemsBase
using RecurrenceAnalysis, Documenter, PyPlot
using DocumenterTools: Themes

# %%
# download the themes
for file in ("juliadynamics-lightdefs.scss", "juliadynamics-darkdefs.scss", "juliadynamics-style.scss")
    download("https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/$file", joinpath(@__DIR__, file))
end
# create the themes
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "juliadynamics-style.scss"), String)
    theme = read(joinpath(@__DIR__, "juliadynamics-$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "juliadynamics-$(w).scss"), header*"\n"*theme)
end
# compile the themes
Themes.compile(joinpath(@__DIR__, "juliadynamics-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "juliadynamics-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

PyPlot.ioff()

makedocs(
modules=[DynamicalSystems, ChaosTools, DynamicalSystemsBase,
         DelayEmbeddings, RecurrenceAnalysis],
doctest=false,
sitename= "DynamicalSystems.jl",
root = @__DIR__,
format = Documenter.HTML(
    prettyurls = CI,
    assets = [
        "assets/logo.ico",
        asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
    collapselevel = 1,
    ),
pages = [
    "Introduction" => "index.md",
    "Contents" => "contents.md",
    "Dynamical System Definition" => "ds/general.md",
    "Predefined Dynamical Systems" => "ds/predefined.md",
    "Numerical Data" => "embedding/dataset.md",
    "DelayEmbeddings" => [
        "Delay Coordinates Embedding" => "embedding/reconstruction.md",
        "Optimal DCE Parameters" => "embedding/estimate.md",
    ],
    "ChaosTools" => [
       "Orbit Diagrams & PSOS" => "chaos/orbitdiagram.md",
       "Lyapunov Exponents" => "chaos/lyapunovs.md",
       "Detecting & Categorizing Chaos" => "chaos/chaos_detection.md",
       "Entropies and Dimensions" => "chaos/entropies.md",
       "Nonlinear Timeseries Analysis" => "chaos/nlts.md",
       "Periodicity" => "chaos/periodicity.md",
       "Choosing a solver" => "chaos/choosing.md",
    ],
    "RecurrenceAnalysis" => [
        "Recurrence Plots" => "rqa/rplots.md",
        "Recurrence Quantification Analysis" => "rqa/quantification.md",
        "Windowed RQA" => "rqa/windowed.md",
    ],
    "Advanced Documentation" => "advanced.md",
    "Contributor Guide" => "contributors_guide.md",
],
)

close("all")

if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/DynamicalSystems.jl.git",
        target = "build",
        push_preview = true
    )
end

PyPlot.ion()
