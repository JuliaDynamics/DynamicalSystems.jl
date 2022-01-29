# Nonlinear Timeseries Analysis

## Broomhead-King Coordinates
```@docs
broomhead_king
```
---
This alternative/improvement of the traditional delay coordinates can be a very
powerful tool. An example where it shines is noisy data where there is the effect
of superficial dimensions due to noise.

Take the following example where we produce noisy data from a system and then use
Broomhead-King coordinates as an alternative to "vanilla" delay coordinates:

```@example MAIN
using DynamicalSystems, CairoMakie

ds = Systems.gissinger()
data = trajectory(ds, 1000.0, Δt = 0.05)
x = data[:, 1]

L = length(x)
s = x .+ 0.5rand(L) #add noise

U, S = broomhead_king(s, 40)
summary(U)
```

Now let's simply compare the above result with the one you get from doing a "standard" call to [`embed`](@ref):
```@example MAIN
fig = Figure()
ax, = lines(fig[1,1], U[:, 1], U[:, 2])
ax.title = "Broomhead-King of s"

R = embed(s, 2, 30)
ax2, = lines(fig[1,2], columns(R)...)
ax2.title = "2D embedding of s"
fig
```

we have used the same system as in the [Delay Coordinates Embedding](@ref) example, and picked the optimal
delay time of `τ = 30` (for same `Δt = 0.05`). Regardless, the vanilla delay coordinates is much worse than the Broomhead-King coordinates.

## DyCA - Dynamical Component Analysis
```@docs
dyca
```

## Return time statistics
```@docs
mean_return_times
exit_entry_times
```

## Nearest Neighbor Prediction
Nearest neighbor timeseries prediction is a method commonly listed under nonlinear timeseries analysis.
This is not part of DynamicalSystems.jl, because in JuliaDynamics we have a dedicated package for this, [TimeseriesPrediction.jl](https://juliadynamics.github.io/TimeseriesPrediction.jl/dev/).

## Timeseries Surrogates
Timeseries surrogates is a method commonly listed under nonlinear timeseries analysis.
This is not part of DynamicalSystems.jl, because in JuliaDynamics we have a dedicated, small, and self-contained package for this, [TimeseriesSurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/dev/).
