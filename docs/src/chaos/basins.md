# Basins of Attraction

In DynamicalSystems.jl we provide performant methods for estimating basins of attraction of various attractors.
The performance of these methods comes from a _constraint on a 2D plane_, as you will see below.

```@docs
basins_map2D
basins_general
```

## Stroboscopic map example
First define a dynamical system on the plane, for example with a *stroboscopic* map or Poincaré section. For example we can set up an dynamical system with a stroboscopic map using a periodically driven 2D continuous dynamical system, like the Duffing oscillator:
```@example MAIN
using DynamicalSystems
ω=1.; F = 0.2
ds =Systems.duffing([0.1, 0.25]; ω = ω, f = F, d = 0.15, β = -1)
integ = integrator(ds; reltol=1e-8)
```

Now we define the grid of ICs that we want to analyze and launch the procedure:

```@example MAIN
xg = range(-2.2,2.2,length=200)
yg = range(-2.2,2.2,length=200)
basin, attractors = basins_map2D(xg, yg, integ; T=2π/ω)
```

```@example MAIN
using PyPlot
figure()
pcolormesh(xg, yg, basin')
savefig("stroboscopic_basin.png"); nothing # hide
```
![](stroboscopic_basin.png)

## Poincaré map example

A Poincaré map of a 3D continuous system is a 2D discrete system and can be directly passed into [`basins_map2D`](@ref).
```@example MAIN
ds = Systems.rikitake(μ = 0.47, α = 1.0)
plane = (3, 0.0)
pmap = poincaremap(ds, (3, 0.), Tmax=1e6;
    idxs = 1:2, rootkw = (xrtol = 1e-8, atol = 1e-8), reltol=1e-9
)
```

Once the Poincaré map has been created, we simply call [`basins_map2D`](@ref)
```@example MAIN
xg=range(-6.,6.,length=200)
yg=range(-6.,6.,length=200)
basin, attractors = basins_map2D(xg, yg, pmap)

figure()
pcolormesh(xg, yg, basin')
savefig("poincare_basin.png"); nothing # hide
```
![](poincare_basin.png)


## Discrete system example
The process to compute the attraction basins of a discrete 2D dynamical system is trivial,
as one passes its integrator directly into [`basins_map2D`](@ref)

```@example MAIN
function newton_map(dz, z, p, n)
    z1 = z[1] + im*z[2]
    dz1 = f(z1, p[1])/df(z1, p[1])
    z1 = z1 - dz1
    dz[1]=real(z1)
    dz[2]=imag(z1)
    return
end
f(x, p) = x^p - 1
df(x, p)= p*x^(p-1)

# dummy Jacobian function to keep the initializator quiet
function newton_map_J(J,z0, p, n)
   return
end

ds = DiscreteDynamicalSystem(newton_map,[0.1, 0.2], [3.0], newton_map_J)
integ  = integrator(ds)

xg=range(-1.5,1.5,length=200)
yg=range(-1.5,1.5,length=200)

basin, attractors  = basins_map2D(xg, yg, integ)
figure
pcolormesh(xg, yg, basin')
savefig("newton_basin.png"); nothing # hide
```
![](newton_basin.png)

## Basins in Higher Dimensions
When it is not so simple to define a 2D stroboscopic map or Poincaré map, which is the case in continuous dynamical systems of higher dimensionality, you can always try the general method [`basins_general`](@ref). It is slower and may requires some tuning.
The algorithm looks for attractors on a 2D grid.
The initial conditions are set on this grid and all others variables are set to zero by default.

### Example

```@example MAIN
ds = Systems.magnetic_pendulum(γ=1, d=0.2, α=0.2, ω=0.8, N=3)
xg=range(-4, 4, length=150)
yg=range(-4, 4, length=150)
@time basin, attractors = basins_general(xg, yg, ds; idxs = 1:2, reltol = 1e-9)
attractors
```
Alright, so far so good, we found 3 attractors (the 3 points of the magnetic pendulum).
Let's visualize this beauty now

```@example MAIN
pcolormesh(xg, yg, basin')
savefig("magnetic_basin.png"); nothing # hide
```
![](magnetic_basin.png)
