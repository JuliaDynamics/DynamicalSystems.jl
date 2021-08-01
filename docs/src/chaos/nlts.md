# Nonlinear Timeseries Analysis

## Numerical Lyapunov Exponent
Given any timeseries, one can first [`embed`](@ref) it using
delay coordinates, and then calculate a maximum
Lyapunov exponent for it. This is done
with
```@docs
lyapunov_from_data
```
---
The function `lyapunov_from_data` has a total of 4 different approaches for the algorithmic process, by
combining 2 types of distances with 2 types of neighborhoods.

### Example of Numerical Lyapunov computation
```@example MAIN
using DynamicalSystems, PyPlot

ds = Systems.henon()
data = trajectory(ds, 100000)
x = data[:, 1] #fake measurements for the win!

ks = 1:20
ℜ = 1:10000
fig = figure(figsize=(10,6))

for (i, di) in enumerate([Euclidean(), Cityblock()])
    subplot(1, 2, i)
    ntype = NeighborNumber(2)
    title("Distance: $(di)", size = 18)
    for D in [2, 4, 7]
        R = embed(x, D, 1)
        E = lyapunov_from_data(R, ks;
        refstates = ℜ, distance = di, ntype = ntype)
        Δt = 1
        λ = linear_region(ks.*Δt, E)[2]
        # gives the linear slope, i.e. the Lyapunov exponent
        plot(ks .- 1, E .- E[1], label = "D=$D, λ=$(round(λ, digits = 3))")
        legend()
        tight_layout()
    end
end
tight_layout(pad=0.3); fig
```

### Bad Time-axis (`ks`) length

!!! danger "Large `ks`"
    This simply cannot be stressed enough! It is just too easy to overshoot
    the range at which the exponential expansion region is valid!

Let's revisit the example of the previous section:
```@example MAIN
ds = Systems.henon()
data = trajectory(ds, 100000)
x = data[:, 1]
length(x)
```
The timeseries of such length could be considered big. A time length of 100 seems
very small. Yet it turns out it is way too big! The following
```@example MAIN
ks = 1:100
R = embed(x, 2, 1)
E = lyapunov_from_data(R, ks, ntype = NeighborNumber(2))
fig = figure()
plot(ks .- 1, E .- E[1])
title("Lyapunov: $(linear_region(ks, E)[2])")
fig.tight_layout(pad=0.3); fig
```

Notice that even though this value
for the Lyapunov exponent is correct, it happened to be correct simply due to the
jitter of the saturated region. Since the saturated region is much bigger
than the linear scaling region, if it wasn't that jittery the function
[`linear_region`](@ref) would not give the scaling of the linear region, but instead
a slope near 0! (or if you were to give bigger tolerance as a keyword argument)

### Case of a Continuous system
The process for continuous systems works identically with discrete, but one must be
a bit more thoughtful when choosing parameters. The following example helps the users get familiar with the process:
```@example MAIN
using DynamicalSystems, PyPlot

ds = Systems.lorenz()
# create a timeseries of 1 dimension
Δt = 0.05
x = trajectory(ds, 1000.0; Δt)[:, 1]
```

We know that we have to use much bigger `ks` than `1:20`, because this is a continuous case! (See reference given in `lyapunov_from_dataspectrum`)
```@example MAIN
ks1 = 0:200
```
and in fact it is even better to not increment the `ks` one by one but instead do
```@example MAIN
ks2 = 0:4:200
```
Now we plot some example computations
```@example MAIN
fig = figure()
ntype = NeighborNumber(5) #5 nearest neighbors of each state

for d in [4, 8], τ in [7, 15]
    r = embed(x, d, τ)

    # E1 = lyapunov_from_data(r, ks1; ntype)
    # λ1 = linear_region(ks1 .* Δt, E1)[2]
    # plot(ks1,E1.-E1[1], label = "dense, d=$(d), τ=$(τ), λ=$(round(λ1, 3))")

    E2 = lyapunov_from_data(r, ks2; ntype)
    λ2 = linear_region(ks2 .* Δt, E2)[2]
    plot(ks2,E2.-E2[1], label = "d=$(d), τ=$(τ), λ=$(round(λ2, digits = 3))")
end

legend()
xlabel("k (0.05×t)")
ylabel("E - E(0)")
title("Continuous Reconstruction Lyapunov")
fig.tight_layout(pad=0.3); fig
```

As you can see, using `τ = 15` is not a great choice! The estimates with
`τ = 7` though are very good (the actual value is around `λ ≈ 0.89...`).

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
using DynamicalSystems, PyPlot

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
fig=figure(figsize= (10,6))
subplot(1,2,1)
plot(U[:, 1], U[:, 2])
title("Broomhead-King of s")

subplot(1,2,2)
R = embed(s, 2, 30)
plot(columns(R)...; color = "C3")
title("2D embedding of s")
fig.tight_layout(pad=0.3)
```

we have used the same system as in the [Delay Coordinates Embedding](@ref) example, and picked the optimal
delay time of `τ = 30` (for same `Δt = 0.05`). Regardless, the vanilla delay coordinates is much worse than the Broomhead-King coordinates.

## Nearest Neighbor Prediction
Nearest neighbor timeseries prediction is a method commonly listed under nonlinear timeseries analysis.
This is not part of DynamicalSystems.jl, because in JuliaDynamics we have a dedicated package for this, [TimeseriesPrediction.jl](https://juliadynamics.github.io/TimeseriesPrediction.jl/dev/).

## DyCA - Dynamical Component Analysis
```@docs
dyca
```
