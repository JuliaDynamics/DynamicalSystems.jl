cd(@__DIR__)

using DynamicalSystems
# TODO: cross-reference docstirngs directly into package's repos!

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

import Literate
Literate.markdown(
    joinpath(@__DIR__, "src", "tutorial.jl"), joinpath(@__DIR__, "src");
    credit = false
)

pages =  [
    "Introduction" => "index.md",
    "Overarching tutorial" => "tutorial.md",
    "Contents" => "contents.md",
    "Animations, GUIs, Visuals" => "visualizations.md",
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
    PredefinedDynamicalSystems;
    authors = "George Datseris <datseris.george@gmail.com>",
    expandfirst = ["index.md"],
    # We need to remove the cross references because we don't list here
    # the whole `DynamicalSystem` API...
    warnonly = [:doctest, :missing_docs, :cross_references],
)
