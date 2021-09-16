# Lyapunov Exponents
Lyapunov exponents measure exponential rates of separation of nearby trajectories in the flow
of a dynamical system. The [Wikipedia](https://en.wikipedia.org/wiki/Lyapunov_exponent) and the [Scholarpedia](http://www.scholarpedia.org/article/Lyapunov_exponent) entries have a lot of valuable information about the history and usage of these quantities.

!!! info "Performance depends on the solver"
    Notice that the performance of functions that use `ContinuousDynamicalSystem`s depend crucially on the chosen solver. Please see the documentation page on [Choosing a solver](@ref) for an in-depth discussion.

## Concept of the Lyapunov exponent
Before providing the documentation of the offered functionality, it is good to demonstrate exactly *what* are the Lyapunov exponents.

For chaotic systems, nearby trajectories separate in time exponentially fast (while for stable systems they come close exponentially fast). This happens at least for small separations, and is demonstrated in the following sketch:

![](lyapunov.png).

In this sketch ``\lambda`` is the maximum Lyapunov exponent (and in general a system has as many exponents as its dimensionality).

Let's demonstrate these concepts using a real system, the Henon map:
```math
\begin{aligned}
x_{n+1} &= 1 - ax_n^2 + y_n \\
y_{n+1} &= bx_n
\end{aligned}
```
Let's get a trajectory
```@example MAIN
using DynamicalSystems, PyPlot
henon = Systems.henon()
tr1 = trajectory(henon, 100)
summary(tr1)
```
and create one more trajectory that starts very close to the first one
```@example MAIN
u2 = get_state(henon) + (1e-9 * ones(dimension(henon)))
tr2 = trajectory(henon, 100, u2)
summary(tr2)
```

We now want to demonstrate how the distance between these two trajectories increases with time:
```@example MAIN
using LinearAlgebra: norm

fig = figure()

# Plot the x-coordinate of the two trajectories:
ax1 = subplot(2,1,1)
plot(tr1[:, 1], alpha = 0.5)
plot(tr2[:, 1], alpha = 0.5)
ylabel("x")

# Plot their distance in a semilog plot:
ax2 = subplot(2,1,2, sharex = ax1)
d = [norm(tr1[i] - tr2[i]) for i in 1:length(tr2)]
ylabel("d"); xlabel("n"); semilogy(d);
fig.tight_layout(pad=0.3); fig
```

The *initial* slope of the `d` vs `n` plot (before the curve saturates) is approximately the maximum Lyapunov exponent!

## Lyapunov Spectrum

The function `lyapunovspectrum` calculates the entire spectrum of the Lyapunov
exponents of a system:
```@docs
lyapunovspectrum
```
---
As you can see, the documentation string is detailed and self-contained. For example,
the Lyapunov spectrum of the [folded towel map](http://www.scholarpedia.org/article/Hyperchaos)
is calculated as:
```@example MAIN
using DynamicalSystems

ds = Systems.towel()
λλ = lyapunovspectrum(ds, 10000)
```
Similarly, for a continuous system, e.g. the Lorenz system, you would do:
```@example MAIN
lor = Systems.lorenz(ρ = 32.0) #this is not the original parameter!
λλ = lyapunovspectrum(lor, 10000, Δt = 0.1)
```

`lyapunovspectrum` is also very fast:
```julia
using BenchmarkTools
ds = Systems.towel()
@btime lyapunovspectrum($ds, 2000);
```
```
  237.226 μs (45 allocations: 4.27 KiB)
```

Here is an example of plotting the exponents of the Henon map for various parameters:
```@example MAIN
using DynamicalSystems, PyPlot

he = Systems.henon()
as = 0.8:0.005:1.225; λs = zeros(length(as), 2)
for (i, a) in enumerate(as)
    set_parameter!(he, 1, a)
    λs[i, :] .= lyapunovspectrum(he, 10000; Ttr = 500)
end

fig = figure()
plot(as, λs); xlabel("\$a\$"); ylabel("\$\\lambda\$")
fig.tight_layout(pad=0.3); fig
```



## Maximum Lyapunov Exponent
It is possible to get only the maximum Lyapunov exponent simply by giving
`1` as the third argument of [`lyapunovspectrum`](@ref). However, there is a second algorithm that allows you to do the same thing, which is offered by the function `lyapunov`:
```@docs
lyapunov
```
---
For example:
```@example MAIN
using DynamicalSystems, PyPlot
henon = Systems.henon()
λ = lyapunov(henon, 10000, d0 = 1e-7, upper_threshold = 1e-4, Ttr = 100)
```

The same is done for continuous systems:
```@example MAIN
lor = Systems.lorenz(ρ = 32)
λ = lyapunov(lor, 10000.0, Δt = 10.0, Ttr = 100.0)
```

## Local Growth Rates
```@docs
local_growth_rates
```
Here is a simple example using the Henon map
```@example MAIN
using Statistics

ds = Systems.henon()
points = trajectory(ds, 2000; Ttr = 100)

λlocal = local_growth_rates(ds, points; Δt = 1)

λmeans = mean(λlocal; dims = 2)
λstds = std(λlocal; dims = 2)
fig = figure()
x, y = columns(points)
scatter(x, y; c=vec(λmeans), alpha = 0.5)
colorbar()
fig.tight_layout(pad=0.3); fig
```


## Lyapunov exponent from data
```@docs
lyapunov_from_data
```

### Example
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
