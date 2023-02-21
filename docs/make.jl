cd(@__DIR__)

using DynamicalSystems

# Reexported
using ComplexityMeasures,
    RecurrenceAnalysis,
    DelayEmbeddings,
    ChaosTools,
    DynamicalSystemsBase,
    StateSpaceSets,
    Attractors,
    FractalDimensions,
    TimeseriesSurrogates

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")


# Also bring in visualizations from interactive dynamics docs:
using InteractiveDynamics
infile = joinpath(pkgdir(InteractiveDynamics), "docs", "src", "dynamicalsystems.md")
outfile = joinpath(@__DIR__, "src", "dynamicalsystems_interactive.md")
cp(infile, outfile)

pages =  [
    "Introduction" => "index.md",
    "Contents" => "contents.md",
    "Dynamical systems" => [
        "Dynamical System Definition" => "ds/general.md",
        "Predefined Dynamical Systems" => "ds/predefined.md",
        "Numerical Data" => "embedding/dataset.md",
        "Integrators" => "ds/integrators.md",
    ],
    "DelayEmbeddings" => [
        "Delay Coordinates Embedding" => "embedding/reconstruction.md",
        "Traditional Optimal Embedding" => "embedding/traditional.md",
        "Unified Optimal Embedding" => "embedding/unified.md",
        ],
    "Entropies" => [
        "Entropies & Probabilities" => "entropies/api.md",
        "Probabilities Estimators" => "entropies/estimators.md",
    ],
    "ChaosTools" => [
       "Orbit Diagrams & PSOS" => "chaos/orbitdiagram.md",
       "Lyapunov Exponents" => "chaos/lyapunovs.md",
       "Detecting & Categorizing Chaos" => "chaos/chaos_detection.md",
       "Fractal Dimension" => "chaos/fractaldim.md",
       "Nonlinear Timeseries Analysis" => "chaos/nlts.md",
       "Fixed points & Periodicity" => "chaos/periodicity.md",
       "Attractors, Basins, Tipping Points" => "chaos/basins.md",
    ],
    "RecurrenceAnalysis" => [
        "Recurrence Plots" => "rqa/rplots.md",
        "Recurrence Quantification Analysis" => "rqa/quantification.md",
        "Windowed RQA" => "rqa/windowed.md",
        "Recurrence Networks" => "rqa/networks.md",
    ],
    "Interactive GUIs" => "dynamicalsystems_interactive.md",
    "Contributor Guide" => "contributors_guide.md",
]

build_docs_with_style(pages, TransitionIndicators;
    authors = "George Datseris <datseris.george@gmail.com>",
    expandfirst = ["index.md"], #  this is the first script that loads colorscheme
)
