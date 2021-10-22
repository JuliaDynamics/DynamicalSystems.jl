# Delay Coordinates Embedding

A timeseries recorded in some manner from a dynamical system can be used to gain information about the dynamics of the entire state space of the system. This can be done by constructing a new state space from the timeseries. One method that can do this is what is known as [delay coordinates embedding](https://en.wikipedia.org/wiki/Takens%27_theorem) or delay coordinates *reconstruction*.

## Timeseries embedding

Delay embeddings are done through `embed`:
```@docs
embed
```

Here are some examples of embedding a 3D continuous chaotic system:
```@example MAIN
using DynamicalSystems, CairoMakie

ds = Systems.gissinger(ones(3))
data = trajectory(ds, 1000.0, Δt = 0.05)
xyz = columns(data)

fig = Figure(resolution = (1000, 800))
k = 1
for i in 1:3
    for (j,τ) in enumerate([5, 30, 100])
        R = embed(xyz[i], 2, τ)
        ax = Axis(fig[i,j])
        lines!(ax, R[:, 1], R[:, 2])
        ax.title = "var = $i, τ = $τ"
        global k+=1
    end
end

fig[:, 0] = Label("2D reconstructed space")
fig
```

!!! note "`τ` and `Δt`"
    Keep in mind that whether a value of `τ` is "reasonable" for continuous systems depends on `Δt`. In the above example the value `τ=30` is good, *only* for the case
    of using `Δt = 0.05`. For shorter/longer `Δt` one has to adjust properly `τ` so that their product `τ*Δt` is the same.

### Embedding Functors
The high level function [`embed`](@ref) utilize a low-level interface for creating embedded vectors on-the-fly. The high level interface simply loops over the low level interface. The low level interface is composed of the following two structures:
```@docs
DelayEmbedding
τrange
```

## Generalized embeddings
```@docs
genembed
GeneralizedEmbedding
```
