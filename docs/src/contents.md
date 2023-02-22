# Contents

The _module_ `DynamicalSystems` re-exports all following packages.
Of course, you could be using these packages directly instead of adding `DynamicalSystems`.
However, doing `using DynamicalSystems` provides the environment all these packages were designed to work together in.

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
ChaosTools
Attractors
```
