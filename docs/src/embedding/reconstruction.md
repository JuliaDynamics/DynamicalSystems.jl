# Delay Coordinates Embedding

A timeseries recorded in some manner from a dynamical system can be used to gain information about the dynamics of the entire state space of the system. This can be done by constructing a new state space from the timeseries. One method that can do this is what is known as [delay coordinates embedding](https://en.wikipedia.org/wiki/Takens%27_theorem) or delay coordinates *reconstruction*.

## Timeseries embedding

Delay embeddings are done through `embed`:
```@docs
embed
```

!!! note "Embedding discretized data values"
    If the data values are very strongly discretized (e.g., integers or floating-point numbers with very small bits), this can result to distances between points in the embedded space being 0. This is problematic for several library functions. Best practice here is to add noise to your original timeseries _before_ embedding, e.g., `s = s .+ 1e-15randn(length(s))`.

---

Here are some examples of embedding a 3D continuous chaotic system:
```@example MAIN
using DynamicalSystems, PyPlot

ds = Systems.gissinger(ones(3))
data = trajectory(ds, 1000.0, Δt = 0.05)

xyz = columns(data)

figure(figsize = (12,10))
k = 1
for i in 1:3
    for τ in [5, 30, 100]
        R = embed(xyz[i], 2, τ)
        ax = subplot(3,3,k)
        plot(R[:, 1], R[:, 2], color = "C$(k-1)", lw = 0.8)
        title("var = $i, τ = $τ")
        global k+=1
    end
end

tight_layout()
suptitle("2D reconstructed space")
subplots_adjust(top=0.9)
savefig("simple_reconstruction.png"); nothing # hide
```
![](simple_reconstruction.png)

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
