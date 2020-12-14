cd(@__DIR__)
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using Pkg
Pkg.activate(@__DIR__)
CI && Pkg.instantiate()

using DynamicalSystems
using Entropies, RecurrenceAnalysis, DelayEmbeddings, ChaosTools, DynamicalSystemsBase
using Neighborhood
using Documenter, PyPlot
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

# %% Build docs
PyPlot.ioff()
cd(@__DIR__)
ENV["JULIA_DEBUG"] = "Documenter"


makedocs(
modules=[DynamicalSystems, ChaosTools, DynamicalSystemsBase,
         DelayEmbeddings, RecurrenceAnalysis, Entropies, Neighborhood],
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
    "Dynamical systems" => [
        "Dynamical System Definition" => "ds/general.md",
        "Predefined Dynamical Systems" => "ds/predefined.md",
        "Numerical Data" => "embedding/dataset.md",
    ],
    "Entropies" => [
        "Entropies & Probabilities" => "entropies/api.md",
        "Probabilities Estimators" => "entropies/estimators.md",
    ],
    "DelayEmbeddings" => [
        "Delay Coordinates Embedding" => "embedding/reconstruction.md",
        "Traditional Optimal Embedding" => "embedding/traditional.md",
        "Unified Optimal Embedding" => "embedding/unified.md",
    ],
    "ChaosTools" => [
       "Orbit Diagrams & PSOS" => "chaos/orbitdiagram.md",
       "Lyapunov Exponents" => "chaos/lyapunovs.md",
       "Detecting & Categorizing Chaos" => "chaos/chaos_detection.md",
       "Fractal Dimension" => "chaos/fractaldim.md",
       "Nonlinear Timeseries Analysis" => "chaos/nlts.md",
       "Periodicity & Ergodicity" => "chaos/periodicity.md",
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

if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/DynamicalSystems.jl.git",
        target = "build",
        push_preview = true
    )
end

PyPlot.close("all")
PyPlot.ion()
