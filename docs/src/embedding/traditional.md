# Traditional Optimal Embedding
This page discusses and provides algorithms for estimating optimal parameters to do Delay Coordinates Embedding (DCE) with.

The approaches can be grouped into two schools:
1. **Traditional**, where one tries to find the best value for a delay time `τ` and then an optimal embedding dimension `d`.
2. **Unified**, where at the same time an optimal combination of `τ, d` is found, and is discussed in the [Unified Optimal Embedding](@ref) page.

The independent approach is something "old school", while recent scientific research has shifted almost exclusively to unified approaches.

## Optimal delay time
```@docs
estimate_delay
exponential_decay_fit
```
### Self Mutual Information

```@docs
selfmutualinfo
```

Notice that mutual information between two *different* timeseries x, y exists in JuliaDynamics as well, but in the package [TransferEntropy.jl](https://github.com/JuliaDynamics/TransferEntropy.jl).
It is also trivial to define it yourself using [`genentropy`](@ref) by doing

```julia
function mutualinfo(x, y, est; base = 2, α = 1)
    X = genentropy(Dataset(x), est; base, α)
    Y = genentropy(Dataset(y), est; base, α)
    XY = genentropy(Dataset(x, y), est; base, α)
    return X + Y - XY
end
```

## Optimal embedding dimension
```@docs
optimal_traditional_de
delay_afnn
delay_ifnn
delay_fnn
delay_f1nn
DelayEmbeddings.stochastic_indicator
```

## Example
```@example MAIN
using DynamicalSystems, CairoMakie

ds = Systems.roessler()
# This trajectory is a chaotic attractor with fractal dim ≈ 2
# therefore the set needs at least embedding dimension of 3
tr = trajectory(ds, 1000.0; Δt = 0.05)
x = tr[:, 1]

dmax = 7
fig = Figure()
ax = Axis(fig[1,1])
for (i, method) in enumerate(["afnn", "fnn", "f1nn", "ifnn"])
    # Plot statistic used to estimate optimal embedding
    # as well as the automated output embedding
    𝒟, τ, E = optimal_traditional_de(x, method; dmax)
    lines!(1:dmax, E; label = method, marker = :circle, color = COLORS[i])
    optimal_d = size(𝒟, 2)
    scatter!(optimal_d, E[optimal_d]; marker = :rect, color = COLORS[i])
end
axislegend(ax)
ax.xlabel = "embedding dimension"
ax.ylabel = "estimator"
fig
```
