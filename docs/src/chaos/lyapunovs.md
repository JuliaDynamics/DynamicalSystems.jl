# Lyapunov Exponents
Lyapunov exponents measure rates of separation of nearby trajectories in the flow
of a dynamical system. The [Wikipedia](https://en.wikipedia.org/wiki/Lyapunov_exponent) and the [Scholarpedia](http://www.scholarpedia.org/article/Lyapunov_exponent) entries have a lot of valuable information about the history and usage of these quantities.

This page treats systems where the equations of motion are known. If instead
you have numerical data, see the [nonlinear timeseries analysis page](nlts).

## Lyapunov Spectrum
The function `lyapunovs` calculates the entire spectrum of the Lyapunov
exponents of a system:
```@docs
lyapunovs
```
---
As you can see, the documentation string is detailed and self-contained. For example,
the lyapunov spectrum of the [folded towel map](http://www.scholarpedia.org/article/Hyperchaos)
is calculated as:
```julia
using DynamicalSystems

ds = Systems.towel()
λλ = lyapunovs(ds, 10000)
```
```julia
[0.432535, 0.372184, -3.29683]
```
Similarly, for a continuous system, e.g. the Lorenz system, you would do:
```julia
lor = Systems.lorenz(ρ = 32.0) #this is not the original parameter!
λλ = lyapunovs(lor, 10000, dt = 0.1)
```
```julia
[0.985688, 0.00271333, -14.6551]
```
`lyapunovs` is also very fast:
```julia
using BenchmarkTools
ds = Systems.towel()
@btime lyapunovs($ds, 1000);
```
```
  239.349 μs (178 allocations: 11.34 KiB)
```

## Maximum Lyapunov Exponent
The function `lyapunov` calculates the maximum lyapunov exponent of a system, more efficiently than getting the first result of `lyapunovs`:
```@docs
lyapunov
```
---
For example:
```julia
using DynamicalSystems
henon = Systems.henon()
λ = lyapunov(henon, 10000, d0 = 1e-7, upper_threshold = 1e-4, Ttr = 100)
```
```
0.42007471604734054
```

The same is done for continuous systems:
```julia
ross = Systems.roessler(a = 0.1, b = 0.1, c = 14.0) #not original parameters
λ = lyapunov(ross, 100000, dt = 10.0, Ttr = 100.0)
```
```
0.07127399326422117
```
