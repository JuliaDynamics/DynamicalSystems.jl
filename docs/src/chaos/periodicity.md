# Fixed points & Periodicity

## Fixed points
```@docs
fixedpoints
```
A rather simple example of the fixed points can be demonstrated using E.g., the Lorenz-63 system, whose fixed points can be easily calculated analytically to be
the following three
```math
(0,0,0) \\
\left( \sqrt{\beta(\rho-1)}, \sqrt{\beta(\rho-1)}, \rho-1 \right) \\
\left( -\sqrt{\beta(\rho-1)}, -\sqrt{\beta(\rho-1)}, \rho-1 \right) \\ 
```

So, let's calculate
```@example MAIN
using DynamicalSystems
ρ, β = 30.0, 10/3
ds = Systems.lorenz(; ρ, β)
# Define the box within which to find fixed points:
x = -20..20
y = -20..20
z = 0.0..40
box = x × y × z

fp, eigs, stable = fixedpoints(ds, box)
fp
```
and compare this with the analytic ones:

```@example MAIN
lorenzfp(ρ, β) = [
    SVector(0, 0, 0.0),
    SVector(sqrt(β*(ρ-1)), sqrt(β*(ρ-1)), ρ-1),
    SVector(-sqrt(β*(ρ-1)), -sqrt(β*(ρ-1)), ρ-1),
]

lorenzfp(ρ, β)
```

## Stable and Unstable Periodic Orbits of Maps
Chaotic behavior
of low dimensional dynamical systems is affected by the position and the stability properties of the [periodic orbits](http://www.scholarpedia.org/article/Unstable_periodic_orbits) of a dynamical system.

Finding unstable (or stable) periodic orbits of a discrete mapping analytically
rapidly becomes impossible for higher orders of fixed points.
Fortunately there is a numeric algorithm due to
Schmelcher & Diakonos which allows such a computation. Notice that even though
the algorithm can find stable fixed points, it is mainly aimed at *unstable* ones.

The functions `periodicorbits` and `lambdamatrix` implement the algorithm:
```@docs
periodicorbits
lambdamatrix
lambdaperms
```

### Standard Map example
For example, let's find the fixed points of the [`Systems.standardmap`](@ref) of order 2, 3, 4, 5, 6
and 8. We will use all permutations for the `signs` but only one for the `inds`.
We will also only use one `λ` value, and a 21×21 density of initial conditions.

First, initialize everything
```@example MAIN
using DynamicalSystems

ds = Systems.standardmap()
xs = range(0, stop = 2π, length = 11); ys = copy(xs)
ics = [SVector{2}(x,y) for x in xs for y in ys]

# All permutations of [±1, ±1]:
singss = lambdaperms(2)[2] # second entry are the signs

# I know from personal research I only need this `inds`:
indss = [[1,2]] # <- must be container of vectors!

λs = 0.005 # <- only this allowed to not be vector (could also be vector)

orders = [2, 3, 4, 5, 6, 8]
ALLFP = Dataset{2, Float64}[];
```
Then, do the necessary computations for all orders

```@example MAIN
for o in orders
    FP = periodicorbits(ds, o, ics, λs, indss, singss)
    push!(ALLFP, FP)
end
```

Plot the phase space of the standard map
```@example MAIN
using CairoMakie
iters = 1000
dataset = trajectory(ds, iters)
for x in xs
    for y in ys
        append!(dataset, trajectory(ds, iters, SVector{2}(x, y)))
    end
end

fig = Figure()
ax = Axis(fig[1,1]; xlabel = L"\theta", ylabel = L"p",
    limits = ((xs[1],xs[end]), (xs[1],xs[end]))
)
scatter!(ax, dataset[:, 1], dataset[:, 2]; markersize = 1, color = "black")
fig
```

and finally, plot the fixed points
```@example MAIN
markers = [:diamond, :utriangle, :rect, :pentagon, :hexagon, :circle]

for i in 1:6
    FP = ALLFP[i]
    o = orders[i]
    scatter!(ax, columns(FP)...; marker=markers[i], color = Cycled(i),
        markersize = 30 - 2i, strokecolor = "grey", strokewidth = 1, label = "order $o"
    )
end
axislegend(ax)
fig
```

Okay, this output is great, and we can easily tell that it is correct for many reasons:

1. It is the same [fig. 12 of this publication](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.012914).
2. Fixed points of order $n$ are also fixed points of order $2n, 3n, 4n, ...$
3. Besides fixed points of previous orders, *original* fixed points of
   order $n$ come in (possible multiples of) $2n$-sized pairs (see e.g. order 5).
   This is a direct consequence of the Poincaré–Birkhoff theorem.

## Estimating the Period

The function [`estimate_period`](@ref) offers ways for estimating the period (either exact for periodic timeseries, or approximate for near-periodic ones) of a given timeseries.
We offer five methods to estimate periods, some of which work on evenly sampled data only, and others which accept any data.
The figure below summarizes this:
![](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/periodestimationmethods.png?raw=true)

```@docs
estimate_period
yin
ChaosTools.difference_function_original
```

### Example
Here we will use a modified FitzHugh-Nagumo system that results in periodic behavior, and then try to estimate its period. First, let's see the trajectory:
```@example MAIN
using DynamicalSystems, CairoMakie

function FHN(u, p, t)
    e, b, g = p
    v, w = u
    dv = min(max(-2 - v, v), 2 - v) - w
    dw = e*(v - g*w + b)
    return SVector(dv, dw)
end

g, e, b  = 0.8, 0.04, 0.0
p0 = [e, b, g]

fhn = ContinuousDynamicalSystem(FHN, SVector(-2, -0.6667), p0)
T, Δt = 1000.0, 0.1
v = trajectory(fhn, T; Δt)[:, 1]
t = 0:Δt:T

fig, ax = lines(0:Δt:T, v)
fig
```

Examining the figure, one can see that the period of the system is around `91` time units. To estimate it numerically let's use some of the methods:
```@example MAIN
estimate_period(v, :autocorrelation, t)
```
```@example MAIN
estimate_period(v, :periodogram, t)
```
```@example MAIN
estimate_period(v, :zerocrossing, t)
```
```@example MAIN
estimate_period(v, :yin, t; sr=round(Int, (1/Δt)), f0_min=0.01)
```