# Lyapunov exponents
Lyapunov exponents measure exponential rates of separation of nearby trajectories in the flow
of a dynamical system. The [Wikipedia](https://en.wikipedia.org/wiki/Lyapunov_exponent) and the [Scholarpedia](http://www.scholarpedia.org/article/Lyapunov_exponent) entries have a lot of valuable information about the history and usage of these quantities.

This page treats systems where the equations of motion are known. If instead
you have numerical data, see the [nonlinear timeseries analysis page](nlts).

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
```@example lyap
using DynamicalSystems, PyPlot
henon = Systems.henon()
tr1 = trajectory(henon, 100)
summary(tr1)
```
and create one more trajectory that starts very close to the first one
```@example lyap
u2 = get_state(henon) + (1e-9 * ones(dimension(henon)))
tr2 = trajectory(henon, 100, u2)
summary(tr2)
```

We now want to demonstrate how the distance between these two trajectories increases with time:
```@example lyap
using LinearAlgebra: norm

figure(figsize=(8,5))

# Plot the x-coordinate of the two trajectories:
ax1 = subplot(2,1,1)
plot(tr1[:, 1], alpha = 0.5)
plot(tr2[:, 1], alpha = 0.5)
ylabel("x")

# Plot their distance in a semilog plot:
ax2 = subplot(2,1,2, sharex = ax1)
d = [norm(tr1[i] - tr2[i]) for i in 1:length(tr2)]
ylabel("d"); xlabel("n"); semilogy(d);
tight_layout() # hide
savefig("demonstration.png"); nothing # hide
```
![](demonstration.png)

The *initial* slope of the `d` vs `n` plot (before the curve saturates) is approximately the maximum Lyapunov exponent!

## Lyapunov Spectrum

The function `lyapunovs` calculates the entire spectrum of the Lyapunov
exponents of a system:
```@docs
lyapunovs
```
---
As you can see, the documentation string is detailed and self-contained. For example,
the Lyapunov spectrum of the [folded towel map](http://www.scholarpedia.org/article/Hyperchaos)
is calculated as:
```@example lyap
using DynamicalSystems

ds = Systems.towel()
λλ = lyapunovs(ds, 10000)
```
Similarly, for a continuous system, e.g. the Lorenz system, you would do:
```@example lyap
lor = Systems.lorenz(ρ = 32.0) #this is not the original parameter!
λλ = lyapunovs(lor, 10000, dt = 0.1)
```

`lyapunovs` is also very fast:
```julia
using BenchmarkTools
ds = Systems.towel()
@btime lyapunovs($ds, 2000);
```
```
  237.226 μs (45 allocations: 4.27 KiB)
```

Here is an example of plotting the exponents of the Henon map for various parameters:
```@example lyap
using DynamicalSystems, PyPlot

he = Systems.henon()
as = 0.8:0.005:1.225; λs = zeros(length(as), 2)
for (i, a) in enumerate(as)
    set_parameter!(he, 1, a)
    λs[i, :] .= lyapunovs(he, 10000; Ttr = 500)
end

figure()
plot(as, λs); xlabel("\$a\$"); ylabel("\$\\lambda\$")
tight_layout() # hide
savefig("heλ.png"); nothing # hide
```
![](heλ.png)



## Maximum Lyapunov Exponent
It is possible to get only the maximum Lyapunov exponent simply by giving
`1` as the third argument of [`lyapunovs`](@ref). However, there is a second algorithm that allows you to do the same thing, which is offered by the function `lyapunov`:
```@docs
lyapunov
```
---
For example:
```@example lyap
using DynamicalSystems, PyPlot
henon = Systems.henon()
λ = lyapunov(henon, 10000, d0 = 1e-7, upper_threshold = 1e-4, Ttr = 100)
```

The same is done for continuous systems:
```@example lyap
lor = Systems.lorenz(ρ = 32)
λ = lyapunov(lor, 10000.0, dt = 10.0, Ttr = 100.0)
```
