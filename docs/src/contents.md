# [Contents](@id contents)

When you do `using DynamicalSystems` in your Julia session, the module
re-exports and brings into scope all packages listed in this session.
Of course, you could be using these packages directly instead of adding `DynamicalSystems`.
However, doing `using DynamicalSystems` provides the environment all these packages were designed to work together in, and so we recommend to simply install `DynamicalSystems` and use that.

Re-exported packages:

**Core**
- [`StateSpaceSets`](@ref)
- [`DynamicalSystemsBase`](@ref)

**For observed/measured data**
- [`ComplexityMeasures`](@ref)
- [`RecurrenceAnalysis`](@ref)
- [`DelayEmbeddings`](@ref)
- [`FractalDimensions`](@ref)
- [`TimeseriesSurrogates`](@ref)

**For dynamical system instances**
- [`PredefinedDynamicalSystems`](@ref)
- [`ChaosTools`](@ref)
- [`Attractors`](@ref)


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
```

## For dynamical system instances

```@docs
PredefinedDynamicalSystems
ChaosTools
Attractors
```
