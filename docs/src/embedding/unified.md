# Unified Optimal Embedding
Unified approaches try to create an optimal embedding by in parallel optimizing what combination of delay times and embedding dimensions suits best.

In addition, the unified approaches are the only ones that can accommodate multi-variate inputs. This means that if you have multiple measured input timeseries, you should be able to take advantage of all of them for the best possible embedding of the dynamical system's set.

## An example

### Univariate input

In following we illustrate the most recent unified optimal embedding method, called PECUZAL, on three examples (see [`pecuzal_embedding`](@ref)).
We start with a univariate case, i.e. we only feed in one time series,
here the x-component of the Lorenz system.  
```@example MAIN
using DynamicalSystems

lo = Systems.lorenz([1.0, 1.0, 50.0])
tr = trajectory(lo, 100; Δt = 0.01, Ttr = 10)

s = vec(tr[:, 1]) # input timeseries = x component of Lorenz
theiler = estimate_delay(s, "mi_min") # estimate a Theiler window
Tmax = 100 # maximum possible delay

Y, τ_vals, ts_vals, Ls, εs = pecuzal_embedding(s; τs = 0:Tmax , w = theiler, econ = true)

println("τ_vals = ", τ_vals)
println("Ls = ", Ls)
println("L_total_uni: $(sum(Ls))")
```
The output reveals that PECUZAL suggests a 3-dimensional embedding out of the
un-lagged time series as the 1st component of the reconstruction, the time
series lagged by 18 samples as the 2nd component and the time series lagged by
9 samples as the 3rd component. In the third embedding cycle there is no *ΔL<0*
and the algorithm breaks. The result after two successful embedding cycles is
the 3-dimensional embedding `Y` which is also returned.
The total obtained decrease of *ΔL* throughout all encountered embedding cycles has been ~ -1.24.


We can also look at *continuity statistic*
```@example MAIN
using CairoMakie

fig = Figure()
ax = Axis(fig[1,1])
lines!(εs[:,1], label="1st emb. cycle")
scatter!([τ_vals[2]], [εs[τ_vals[2],1]])
lines!(εs[:,2], label="2nd emb. cycle")
scatter!([τ_vals[3]], [εs[τ_vals[3],2]])
lines!(εs[:,3], label="3rd emb. cycle")
ax.title = "Continuity statistics PECUZAL Lorenz"
ax.xlabel = "delay τ"
ax.ylabel = "⟨ε⋆⟩"
axislegend(ax)
fig
```

The picked delay values are marked with filled circles. As already mentioned, the
third embedding cycle did not contribute to the embedding, i.e. there has been
no delay value chosen.

### Multivariate input

Similar to the approach in the preceding example, we now highlight the capability
of the PECUZAL embedding method for a multivariate input. The idea is now to feed
in all three time series to the algorithm, even though this is a very
far-from-reality example. We already have an adequate representation of the
system we want to reconstruct, namely the three time series from the numerical
integration. But let us see what PECUZAL suggests for a reconstruction.

```julia
# compute Theiler window
w1 = estimate_delay(tr[:,1], "mi_min")
w2 = estimate_delay(tr[:,2], "mi_min")
w3 = estimate_delay(tr[:,3], "mi_min")
w = max(w1,w2,w3)
Y_m, τ_vals_m, ts_vals_m, = pecuzal_embedding(tr; τs = 0:Tmax , w = theiler, econ = true)

println(τ_vals_m)
println(ts_vals_m)
```
```
[0, 12, 0, 79, 64, 53]
[3, 1, 1, 1, 1, 1]
```

PECUZAL returns a 6-dimensional embedding using the un-lagged *z*- and *x*-component
as 1st and 3rd component of the reconstruction vectors, as well as the *x*-component
lagged by 12, 79, 64, and 53 samples. The total decrease of *ΔL* is ~-1.64, and
thus, way smaller compared to the univariate case, as we would expect it. Nevertheless,
the main contribution to this increase is made by the first two embedding cycles.
For surpressing embedding cycles, which yield negligible - but negative - *ΔL*-values
one can use the keyword argument *L_threshold*   
```@example MAIN

Y_mt, τ_vals_mt, ts_vals_mt, Ls_mt , εs_mt = pecuzal_embedding(tr; τs = 0:Tmax , L_threshold = 0.05, w = theiler, econ = true)

println(τ_vals_mt)
println(ts_vals_mt)
```
As you can see here the algorithm stopped already at 3-dimensional embedding.

Let's plot these three components:
```@example MAIN
ts_str = ["x", "y", "z"]

fig = Figure(resolution = (1000,500) )
ax1 = Axis3(fig[1,1], title = "PECUZAL reconstructed")
lines!(ax1, Y_mt[:,1], Y_mt[:,2], Y_mt[:,3]; linewidth = 1.0)
ax1.xlabel = "$(ts_str[ts_vals_mt[1]])(t+$(τ_vals_mt[1]))"
ax1.ylabel = "$(ts_str[ts_vals_mt[2]])(t+$(τ_vals_mt[2]))"
ax1.zlabel = "$(ts_str[ts_vals_mt[3]])(t+$(τ_vals_mt[3]))"

ax2 = Axis3(fig[1,2], title = "original")
lines!(ax2, tr[:,1], tr[:,2], tr[:,3]; linewidth = 1.0, color = COLORS[2])
ax2.xlabel = "x(t)"
ax2.ylabel = "y(t)"
ax2.zlabel = "z(t)"
ax2.azimuth = π/2 + π/4
fig
```

Finally we show what PECUZAL does with a non-deterministic source:

```@example MAIN
using Random

# Dummy input
Random.seed!(1234)
d1 = randn(1000)
d2 = rand(1000)
Tmax = 100
dummy_set = Dataset(d1,d2)

w1 = estimate_delay(d1, "mi_min")
w2 = estimate_delay(d2, "mi_min")
theiler = min(w1, w2)

Y_d, τ_vals_d, ts_vals_d, Ls_d , ε★_d = pecuzal_embedding(dummy_set; τs = 0:Tmax , w = theiler, econ = true)

size(Y_d)
```

So, no (proper) embedding is done.

## All unified algorithms

Several algorithms have been created to implement a unified approach to delay coordinates embedding. You can find some implementations below:
```@docs
pecora
uzal_cost
garcia_almeida_embedding
mdop_embedding
pecuzal_embedding
```

### Low-level functions of unified approach
```@docs
DelayEmbeddings.n_statistic
DelayEmbeddings.beta_statistic
DelayEmbeddings.mdop_maximum_delay
```
