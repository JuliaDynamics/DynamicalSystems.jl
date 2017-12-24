## Orbit Diagrams of Maps
An orbit diagram (also called bifurcation diagram) is a way to visualize the asymptotic
behavior of a map, when a parameter of the system is changed
```@docs
orbitdiagram
```
---

For example, let's compute the famous orbit diagram of the logistic map:
```julia
using DynamicalSystems
using PyPlot

ds = Systems.logistic()
i = 1
pvalues = 2:0.001:4
ics = [rand() for m in 1:10]
n = 50
Ttr = 5000
output = orbitdiagram(ds, i, :r, pvalues; n = n, Ttr = Ttr)

figure()
for (j, p) in enumerate(pvalues)
    plot(p .* ones(output[j]), output[j], lw = 0,
    marker = "o", ms = 0.5, color = "black")
end
xlabel("\$r\$"); ylabel("\$x\$")
```
![Logistic diagram](https://i.imgur.com/BexsS9Y.png)

Notice that if you are using `PyPlot`, the plotting process will be slow, since it is slow at plotting big numbers of points.

The function is not limited to 1D maps, and can be applied just as well to any
discrete system.
```julia
ds = Systems.standardmap()
i = 2

pvalues = 0:0.005:2
ics = [0.001rand(2) for m in 1:10]
n = 50
Ttr = 5000
output = orbitdiagram(ds, i, :k, pvalues; n = n, Ttr = Ttr)

figure()
for (j, p) in enumerate(pvalues)
    plot(p .* ones(output[j]), output[j], lw = 0,
    marker = "o", ms = 0.5, color = "black")
end
xlabel("\$k\$"); ylabel("\$p\$")
```

![Standard map diagram](https://i.imgur.com/f97oYCx.png)

## Poincaré Surface of Section
Also called [Poincaré map](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) is a
technique to reduce a continuous system into a discrete map with 1 less dimension.
We are doing this using the function:
```@docs
poincaresos
```
---

An example of the [Henon-Helies](efinition/predefined/#DynamicalSystemsBase.Systems.henonhelies) system using a quasi-periodic solution
```julia
ds = Systems.henonhelies([0., 0.1, 0.5, 0.])
output = poincaresos(ds, 3, 1000.0, diff_eq_kwargs=Dict(:solver=>Vern9()))

figure()
plot(output[:, 2], output[:, 4], lw = 0.0, marker=".")
xlabel("\$q_2\$"); ylabel("\$p_2\$");
```

![Poincare SOS](https://i.imgur.com/AkAKOh4.png)

Here the surface of section was the (hyper-) plane that $p_1 = 0$. As expected the section is 1-dimensional, because the torus the solution lives in is 2-dimensional. if
we produced the PSOS for much longer times, the result would be a filled line instead
of individual points.

## Producing Orbit Diagrams for Continuous Flows
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

```julia
ds = Systems.shinriki([-2, 0, 0.2])

pvalues = linspace(19,22,201)
i = 1
j = 2
tf = 200.0

de = Dict(:abstol=>1e-9, :reltol => 1e-9) #necessary due to exponential function

output = produce_orbitdiagram(ds, j, i, :R1, pvalues; tfinal = tf,
Ttr = 200.0, diff_eq_kwargs = de, direction = -1, printparams = true)

figure()
for (j, p) in enumerate(pvalues)
    plot(p .* ones(output[j]), output[j], lw = 0,
    marker = "o", ms = 0.5, color = "black")
end
xlabel("\$R_1\$"); ylabel("\$V_1\$")
```

![shinriki period doubling](https://i.imgur.com/Yd1D3Ou.png)
