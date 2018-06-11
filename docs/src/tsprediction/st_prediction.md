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
    Several example scripts can be found in `TimeseriesPrediction/examples`.

```@docs
localmodel_stts
crosspred_stts
STReconstruction
```

## Showcasing Results
### Chaotic Barkley Model cross-prediction
Here we cross-predicting the $V$ field from $U$ in the chaotic Barkley model.
It is defined by:

```math
\begin{align}
\frac{\partial u }{\partial t} =& \frac{1}{\epsilon} u (1-u)\left(u-\frac{v+b}{a}\right) +
 \nabla^2 u \\
\frac{\partial v }{\partial t} =& u - v
\end{align}
```

Embedding parameters are `D=30`, `τ=1`, `B=0` and training length `2000`
![Crossprediction U->V in Barkley](https://i.imgur.com/Q2yKRvB.png).

You can find the script that produced this figure in
`DynamicalSystems/coolanimations/barkley_crosspred.jl`.


### Periodic Barkley Model timeseries prediction
Using different parameters in the Barkley model can produce periodic behavior.
Here is an example of timeseries prediction with embedding parameters
`D=2`, `τ=1`, `B=2`, `k=1` and `c=200`.

Plotting the real evolution, prediction, and error side by side
with `Ttrain = 1000, p = 200` produces:

![Barkley prediction](https://i.imgur.com/ldChwOD.gif)

You can find the script that produced this animation in
`DynamicalSystems/coolanimations/barkley_stts.jl`.
