# Lyapunov Exponents
Lyapunov exponents measure rates of separation of nearby trajectories in the flow
of a dynamical system. The [Wikipedia](https://en.wikipedia.org/wiki/Lyapunov_exponent) and the [Scholarpedia](http://www.scholarpedia.org/article/Lyapunov_exponent) entries have a lot of valuable information about the history and usage of these quantities.

The naming comes after Aleksandr M. Lyapunov, a Russian mathematician/physicist that had major impact on the analysis of the stability of systems.

## Dynamical Systems
For systems where the equations of motion are known, there two functions related with the computation of lyapunov exponents: `lyapunovs` and `lyapunov`.
### Lyapunov Spectrum
The function `lyapunovs` calculates the entire spectrum of the exponents of a system:
```@docs
lyapunovs
```
As you can see, the documentation string is detailed and self-contained. For example,
the lyapunov spectrum of the [folded towel map](http://www.scholarpedia.org/article/Hyperchaos)
is calculated as:
```julia
using DynamicalSystems

ds = Systems.towel()
λλ = lyapunovs(ds, 10000)
# result:
[0.432253, 0.371617, -3.29632]
```
Similarly, for a continuous system, e.g. the Lorenz system, you would do:
```julia
using DynamicalSystems

lor = Systems.lorenz(ρ = 32.0) #this is not the original parameter!
issubtype(typeof(ds), ContinuousDS) # true

λλ = lyapunovs(lor, 10000,
dt = 0.1, diff_eq_kwargs = Dict(:abstol => 1e-9, :reltol => 1e-9))
# result:
[0.999176, 0.000774754, -14.6666]
```

### Maximum Lyapunov Exponent
The function `lyapunov` calculates the maximum lyapunov exponent of a system, much
more efficiently than getting the first result of `lyapunovs`:
```@docs
lyapunov
```
For example:
```julia
using DynamicalSystems

henon = Systems.henon()
λ = lyapunov(henon, 10000, d0 = 1e-7, threshold = 1e-4, Ttr = 100)
# result:
0.42011626111385747
```
The same is done for continuous systems:
```julia
using DynamicalSystems, OrdinaryDiffEq

ross = Systems.roessler(a = 0.1, b = 0.1, c = 14.0) #not original parameters
λ = lyapunov(ross, 10000, dt = 0.5, diff_eq_kwargs = Dict(:solver => Vern8()))
# result:
0.06957484163052223
```

## Numerical Data
TBA
