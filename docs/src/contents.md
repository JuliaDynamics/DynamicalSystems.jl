# [Contents](@id contents)

When you do `using DynamicalSystems` in your Julia session, the module
re-exports and brings into scope all submodules (Julia packages) that compose **DynamicalSystems.jl**. These are listed in this page.
Of course, you could be using these packages directly instead of adding `DynamicalSystems`.
However, doing `using DynamicalSystems` provides the environment all these packages were designed to work together in, and so we recommend to simply install `DynamicalSystems` and use that.

## Exported submodules

The submodules that compose **DynamicalSystems.jl** are the following packages, which are re-exported by `DynamicalSystems`:

**Core**
- [`StateSpaceSets`](@ref)
- [`DynamicalSystemsBase`](@ref)

**For observed/measured data and timeseries**
- [`ComplexityMeasures`](@ref)
- [`RecurrenceAnalysis`](@ref)
- [`DelayEmbeddings`](@ref)
- [`FractalDimensions`](@ref)
- [`TimeseriesSurrogates`](@ref)
- [`SignalDecomposition`](@ref)

**For dynamical system instances**
- [`PredefinedDynamicalSystems`](@ref)
- [`ChaosTools`](@ref)
- [`Attractors`](@ref)
- [`PeriodicOrbits`](@ref)

At the very end of this page, a full list of exported names is presented.


## Core

```@docs
StateSpaceSets
DynamicalSystemsBase
```

## For observed/measured data

```@docs
ComplexityMeasures
RecurrenceAnalysis
DelayEmbeddings
FractalDimensions
TimeseriesSurrogates
SignalDecomposition
```

## For dynamical system instances

```@docs
PredefinedDynamicalSystems
ChaosTools
Attractors
PeriodicOrbits
```

## All exported names

This section lists all exported names of the **DynamicalSystems.jl** library. We do not list their documentation in any way here. This list is only meant as a quantitative listing of features, as well as perhaps helping searching via the search bar. To actually learn how to use all these exported names you need to use above-linked documentation of the respective submodules!

The total exported names are:

```@example MAIN
using DynamicalSystems
all_exported_names = names(DynamicalSystems)
length(all_exported_names)
```

And they are:

```@example MAIN
using DisplayAs
DisplayAs.unlimited(all_exported_names)
```