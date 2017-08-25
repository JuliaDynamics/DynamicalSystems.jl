# Nonlinear Timeseries Analysis
## Attractor Reconstruction
A timeseries can be used to gain information about the dynamics of the entire phase-space,
by reconstructing a phase-spaced from the timeseries. One method that can do this is
what is known as [delay coordinates embedding](https://en.wikipedia.org/wiki/Takens%27_theorem).

In `DynamicalSystems.jl` this is done through the `reconstruct` interface:
```@docs
reconstruct
```
The returned object can be immediately passed into e.g. a method that calculates the
attractor dimension, or simply be used as any matrix. For example:
```julia
using DynamicalSystems
he = Systems.henon()
ts = timeseries(he, 20000)
D1 = information_dim(ts) # around 1.20
x = view(ts, :, 1) # some "recorded" timeseries
R = reconstruct(x, 5, 3)
D2 = information_dim(R) #around 1.17
```
The 2 numbers are very close, which is not too bad given that there wasn't
any care for the value of `τ` that was chosen! More importantly:
```julia
using BenchmarkTools
@btime recostruct($x, 5, 3);
# @btime result:  94.675 ns (6 allocations: 304 bytes)
```
