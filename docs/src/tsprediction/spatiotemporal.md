# Spatiotemporal Timeseries Prediction

An application and extension of [local modeling](tsprediction/localmodels) to
spatiotemporal timeseries.

!!! tip "Examples"
    Several example scripts can be found in `TimeseriesPrediction/examples`. These examples are run in the [examples](stexamples.md) page.


## Spatio-Temporal Embeddings
some info here.
```@docs
AbstractSpatialEmbedding
SpatioTemporalEmbedding
cubic_shell_embedding
light_cone_embedding
PCAEmbedding
```
---

Boundary conditions
```@docs
ConstantBoundary
PeriodicBoundary
```

## Prediction functions
```@docs
temporalprediction
crossprediction
```
