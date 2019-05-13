# Lyapunov exponents
Lyapunov exponents measure exponential rates of separation of nearby trajectories in the flow
of a dynamical system. The [Wikipedia](https://en.wikipedia.org/wiki/Lyapunov_exponent) and the [Scholarpedia](http://www.scholarpedia.org/article/Lyapunov_exponent) entries have a lot of valuable information about the history and usage of these quantities.

This page treats systems where the equations of motion are known. If instead
you have numerical data, see the [nonlinear timeseries analysis page](nlts).

!!! info "Performance depends on the solver"
    Notice that the performance of functions that use `ContinuousDynamicalSystem`s depend crucially on the chosen solver. Please see the documentation page on [Choosing a solver](@ref) for an in-depth discussion.

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

Here is an example of plotting the exponents of the Roessler system for various parameters (using the advanced usage):
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
