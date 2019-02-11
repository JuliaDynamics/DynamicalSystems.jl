## Orbit Diagrams of Maps
An orbit diagram (also called bifurcation diagram) is a way to visualize the asymptotic
behavior of a map, when a parameter of the system is changed
```@docs
orbitdiagram
```
---

For example, let's compute the famous orbit diagram of the logistic map:
```@example orbit
using DynamicalSystems
using PyPlot

ds = Systems.logistic()
i = 1
pvalues = 3:0.001:4
ics = [rand() for m in 1:10]
n = 2000
Ttr = 2000
p_index = 1
output = orbitdiagram(ds, i, p_index, pvalues; n = n, Ttr = Ttr)

L = length(pvalues)
x = Vector{Float64}(undef, n*L)
y = copy(x)
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= pvalues[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end

figure()
PyPlot.title("total points: $(L*n)")
plot(x, y, ls = "None", ms = 0.5, color = "black", marker = "o", alpha = 0.05)
xlim(pvalues[1], pvalues[end]); ylim(0,1)
xlabel("\$r\$"); ylabel("\$x\$")
tight_layout()
savefig("logostic_od.png"); nothing # hide
```
![](logostic_od.png)

Notice that if you are using `PyPlot`, the plotting process will be slow, since it is slow at plotting big numbers of points.

The function is not limited to 1D maps, and can be applied just as well to any
discrete system.

## Poincaré Surface of Section
Also called [Poincaré map](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) is a
technique to reduce a continuous system into a discrete map with 1 less dimension.
We are doing this using the function:
```@docs
poincaresos
```
---

Here is an example of the [Henon-Heiles](https://en.wikipedia.org/wiki/H%C3%A9non%E2%80%93Heiles_system) system showing the mixed nature of the phase space
```@example orbit
using DynamicalSystems, PyPlot

hh = Systems.henonheiles()

plane = (1, 0.0)
u0s = [[0.0, -0.25, 0.42081, 0.0],
[0.0, -0.31596, 0.354461, 0.0591255],
[0.0, 0.1, 0.5, 0.0],
[0.0, -0.0910355, 0.459522, -0.173339],
[0.0, -0.205144, 0.449328, -0.0162098]]

figure()
for u0 in u0s
    psos = poincaresos(hh, plane, 20000.0; u0 = u0)
    scatter(psos[:, 2], psos[:, 4], s = 2.0)
end
xlabel("\$q_2\$"); ylabel("\$p_2\$")
savefig("hhpsos.png"); nothing # hide
```
![](hhpsos.png)

Here the surface of section was the (hyper-) plane that $q_1 = 0$. Some chaotic and regular orbits can be seen in the plot. You can tell the regular orbits apart because they look like a single connected curve. This is the result of cutting a 2-torus by a plane!

---
Finally here is one more example with a more complex hyperplane:
```@example orbit
gis = Systems.gissinger([2.32865, 2.02514, 1.98312])

# Define appropriate hyperplane for gissinger system
const ν = 0.1
const Γ = 0.9 # default parameters of the system

# I want hyperperplane defined by these two points:
Np(μ) = SVector{3}(sqrt(ν + Γ*sqrt(ν/μ)), -sqrt(μ + Γ*sqrt(μ/ν)), -sqrt(μ*ν))
Nm(μ) = SVector{3}(-sqrt(ν + Γ*sqrt(ν/μ)), sqrt(μ + Γ*sqrt(μ/ν)), -sqrt(μ*ν))

# Create hyperplane passing through Np, Nm and 0:
using LinearAlgebra
gis_plane(μ) = [cross(Np(μ), Nm(μ))..., 0]

μ = 0.119
set_parameter!(gis, 1, μ)
figure(figsize = (8,6))
psos = poincaresos(gis, gis_plane(μ), 10000.0, Ttr = 200.0,)
plot3D(columns(psos)..., marker = "o", ls = "None", ms = 2.0);
xlabel("Q"); ylabel("D"); zlabel("V");
savefig("gispsos.png"); nothing # hide
```
![](gispsos.png)


### Stroboscopic Map
A special case of a PSOS is a stroboscopic map, which is defined for non-autonomous
systems with periodic time dependence, like e.g. the [Duffing oscillator](/definition/predefined/#DynamicalSystemsBase.Systems.duffing).

A "cut" through the phase-space can be produced at every period $T = 2\pi/\omega$. There is no
reason to use `poincaresos` for this though, because you can simply use
[`trajectory`](@ref) and get the solution with a certain time sampling rate:
```@example orbit
ds = Systems.duffing(β = -1, ω = 1, f = 0.3) # non-autonomous chaotic system
a = trajectory(ds, 100000.0, dt = 2π) # every period T = 2π/ω
figure()
plot(a[:, 1], a[:, 2], lw = 0, marker ="o", ms = 1)
xlabel("\$x\$"); ylabel("\$\\dot{x}\$")
savefig("duffing.png"); nothing # hide
```
![](duffing.png)


## Producing Orbit Diagrams for Flows
The [`orbitdiagram`](@ref) does not make much sense for continuous systems, besides the
trivial case where the system is at a fixed point. In order for [`orbitdiagram`](@ref) to have meaning one must have a map.

If only there was a way to turn a continuous system into a map... **OH WAIT!** That is
what [`poincaresos`](@ref) does! By performing successive surfaces of section at different parameter values, one can indeed "produce" an orbit diagram for a flow.

We have bundled this process in the following function:
```@docs
produce_orbitdiagram
```
---

For example, we will calculate the orbit diagram of the Shinriki oscillator, a continuous system that undergoes a period doubling route to chaos, much like the logistic map!

```@example orbit
ds = Systems.shinriki([-2, 0, 0.2])

pvalues = range(19, stop = 22, length = 401)
i = 1
plane = (2, 0.0)
tf = 200.0
p_index = 1

output = produce_orbitdiagram(ds, plane, i, p_index, pvalues;
                              tfinal = tf, Ttr = 200.0)

figure()
for (j, p) in enumerate(pvalues)
    plot(fill(p, length(output[j])), output[j], lw = 0,
    marker = "o", ms = 0.2, color = "black")
end
xlabel("\$R_1\$"); ylabel("\$V_1\$")
tight_layout()
savefig("shinriki.png"); nothing # hide
```
![](shinriki.png)
