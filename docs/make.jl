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
cp(infile, outfile; force = true)

pages =  [
    "Introduction" => "index.md",
    "Overarching tutorial" => "tutorial.md",
    "Contents" => "contents.md",
    "Interactive GUIs" => "dynamicalsystems_interactive.md",
    "Contributor Guide" => "contributors_guide.md",
]

build_docs_with_style(pages, DynamicalSystems,
    ComplexityMeasures,
    RecurrenceAnalysis,
    DelayEmbeddings,
    ChaosTools,
    DynamicalSystemsBase,
    StateSpaceSets,
    Attractors,
    FractalDimensions,
    TimeseriesSurrogates,
    InteractiveDynamics;
    authors = "George Datseris <datseris.george@gmail.com>",
    expandfirst = ["index.md"], #  this is the first script that loads colorscheme
)
