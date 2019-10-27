# Periodicity
In this page we describe methods related to the periodic behavior of dynamical systems or univariate timeseries.

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
---
### Standard Map example
For example, let's find the fixed points of the [Standard Map](system_definition/#DynamicalSystems.Systems.standardmap) of order 2, 3, 4, 5, 6
and 8. We will use all permutations for the `signs` but only one for the `inds`.
We will also only use one `λ` value, and a 21×21 density of initial conditions.

First, initialize everything
```@example sm
using DynamicalSystems, PyPlot, StaticArrays

ds = Systems.standardmap()
xs = range(0, stop = 2π, length = 21); ys = copy(xs)
ics = [SVector{2}(x,y) for x in xs for y in ys]

# All permutations of [±1, ±1]:
singss = lambdaperms(2)[2] # second entry are the signs

# I know from personal research I only need this `inds`:
indss = [[1,2]] # <- must be container of vectors!!!

λs = 0.005 # <- only this allowed to not be vector (could also be vector)

orders = [2, 3, 4, 5, 6, 8]
ALLFP = Dataset{2, Float64}[];
```
Then, do the necessary computations for all orders

```@example sm
for o in orders
    FP = periodicorbits(ds, o, ics, λs, indss, singss)
    push!(ALLFP, FP)
end
```

Plot the phase space of the standard map
```@example sm
iters = 1000
dataset = trajectory(ds, iters)
for x in xs
    for y in ys
        append!(dataset, trajectory(ds, iters, SVector{2}(x, y)))
    end
end
figure(figsize = (12,12))
m = Matrix(dataset)
PyPlot.scatter(view(m, :, 1), view(m, :, 2), s= 1, color = "black")
PyPlot.xlim(xs[1], xs[end])
PyPlot.ylim(ys[1], ys[end]);
```

and finally, plot the fixed points
```@example sm
markers = ["D", "^", "s", "p", "h", "8"]
colors = ["b", "g", "r", "c", "m", "grey"]

for i in 1:6
    FP = ALLFP[i]
    o = orders[i]
    PyPlot.plot(columns(FP)...,
    marker=markers[i], color = colors[i], markersize=10.0 + (8-o), linewidth=0.0,
    label = "order $o", markeredgecolor = "yellow", markeredgewidth = 0.5)
end
legend(loc="upper right", framealpha=0.9)
xlabel("\$\\theta\$")
ylabel("\$p\$")
savefig("fixedpoints.png"); nothing # hide
```
![Fixed points of the standard map](fixedpoints.png)

You can confirm for yourself that this is correct, for many reasons:

1. It is the same [fig. 12 of this publication](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.012914).
2. Fixed points of order $n$ are also fixed points of order $2n, 3n, 4n, ...$
3. Besides fixed points of previous orders, *original* fixed points of
   order $n$ come in (possible multiples of) $2n$-sized pairs (see e.g. order 5).
   This is a direct consequence of the Poincaré–Birkhoff theorem.

## Determining Periodicity

Once you have determined that your system is periodic, you might want to find out what its period is.  Fortunately, the function [`estimate_period`](@ref) from `ChaosTools` offers many ways for you to do this, given the system's timeseries as an input.

We offer five methods to estimate periods, some of which work on evenly sampled data only, and others which accept any data.  The figure below summarizes this:
![](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/periodestimationmethods.png?raw=true)

```@docs
estimate_period
```

### Example
Here we will use a modified FitzHugh-Nagumo system that results in periodic behavior, and then try to estimate its period. First, let's see the trajectory:
```@example sm
using DynamicalSystems, PyPlot

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
T, dt = 1000.0, 0.1
v = trajectory(fhn, T; dt = dt)[:, 1]
t = 0:dt:T

figure()
plot(0:dt:T, v)
savefig("fhn_trajectory.png"); nothing # hide
```
![A periodic trajectory](fhn_trajectory.png)

Examining the figure, one can see that the period of the system is around `91` time units. To estimate it numerically let's use some of the methods:
```@example sm
estimate_period(v, :autocorrelation, t)
```
```@example sm
estimate_period(v, :periodogram, t)
```
```@example sm
estimate_period(v, :zerocrossing, t)
```
