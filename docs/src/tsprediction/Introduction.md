# Introduction

Predicting timeseries of chaotic systems can be a very difficult task. Most methods employed for such a feat are typically relying on large neural networks and machine learning. One does not need them though! In the package `TimeseriesPrediction` we are presenting methods that instead take advantage of dynamical systems theory. Many such methods exist, like for example [Cluster weighted modelling](https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527609970.ch3) or [networks of dynamical systems](https://github.com/JuliaDynamics/TimeseriesPrediction.jl/issues/40). The first method included in `TimeseriesPrediction`, which is also simplest one, is called "local modelling".

## Local Modelling

Local modelling predicts timeseries using a delay embedded state space reconstruction.
It finds the nearest neighbors of a query point within this reconstructed space and applies a local model to make a prediction. "Local" model refers to the fact that the images (future points) of the [`neighborhood`](@ref) of a point are the only component used to make a prediction.

In contrast to typical neural networks applications, there is no training happening in this approach. A given timeseries dataset constitutes a pool of points one uses to make predictions from.



## Available Functionality

### Local Models
`TimeseriesPrediction` has two local models (usable in any prediction scheme)
and uses the [`neighborhood`](@ref) types from `DynamicalSystemsBase`. See [`AbstractLocalModel`](@ref) for more details.

### Timeseries Prediction
[`localmodel_tsp`](@ref) predicts the future of one or many (univariate) timeseries. The details are in the [timeseries prediction page](localmodels).

### Spatio-Temporal Timeseries
One of the biggest strengths of `TimeseriesPrediction` is a robust, simple, and feature-rich interface that can predict the evolution of spatio-temporal systems (commonly represented by PDEs or by "map lattices" (coupled maps)). To see the full interface please visit the [spatio-temporal prediction page](spatiotemporal.md). In addition, the [spatio-temporal examples](stexamples.md) page is full of runnable code that displays the capabilities of spatio-temporal prediction of `TimeseriesPrediction`.
