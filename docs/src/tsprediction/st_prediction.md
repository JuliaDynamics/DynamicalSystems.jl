# Spatiotemporal Timeseries Prediction

An application and extension of [local modeling](tsprediction/localmodels) to
spatiotemporal timeseries.

Previously we have only considered relatively small systems.
If $B$-dim. multivariate timeseries were given, complete states would be embedded into a
reconstruction of dimension $D×B$. Spatiotemporal timeseries could in principle be
considered simply be considered in this way but it is far from being ideal.

One of the main reasons is that it does not utilize all available information.
All spatially extended physical systems posses a finite speed at which information travels.
Therefore the future value of any of the variables depends solely on its past and
its immediate spatial neighbors. Instead of trying to reconstruct the state of the whole
system into one vector, we limit ourselves to reconstructing small neighborhoods of all
points that carry enough information to predict one point one timestep into the future.

!!! tip "Examples"
    Several example scripts can be found in _.julia/v0.6/TimeseriesPrediction/examples_.

```@docs
localmodel_stts
crosspred_stts
STReconstruction
```
