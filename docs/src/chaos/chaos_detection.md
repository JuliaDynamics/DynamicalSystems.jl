# Detecting & Categorizing Chaos
Being able to detect and distinguish chaotic from regular behavior is crucial in the study of dynamical systems.
Most of the time a positive maximum [`lyapunov`](@ref) exponent
and a bounded system indicate chaos.

However, the convergence of the Lyapunov exponent can be slow, or even misleading, as the types of chaotic behavior vary with respect to their predictability.
There are many alternatives, some more efficient and some more accurate in characterizing chaotic and regular motion. Some of these methods are included in **DynamicalSystems.jl**.

!!! info "Performance depends on the solver"
    Notice that the performance of functions that use `ContinuousDynamicalSystem`s depend crucially on the chosen solver. Please see the documentation page on [Choosing a solver](@ref) for an in-depth discussion.

## Generalized Alignment Index
"GALI" for sort, is a method that relies on the fact that initially orthogonal deviation vectors tend to align towards the direction of the maximum Lyapunov exponent for chaotic motion.
It is one of the most recent and cheapest methods for distinguishing chaotic and regular behavior, introduced first in 2007 by Skokos, Bountis & Antonopoulos.
```@docs
gali
```


### Discrete Example
We will use 3 coupled standard maps as an example for a discrete system:
```@example MAIN
using DynamicalSystems
M = 3; ks = 3ones(M); Γ = 0.1;
stable = [π, π, π, 0.01, 0, 0] .+ 0.1
chaotic = rand(2M)

ds = Systems.coupledstandardmaps(M, stable; ks, Γ)
```

For this example we will see the behavior of GALI for a stable orbit
```@example MAIN
using CairoMakie
tr = trajectory(ds, 100000)

fig, ax = scatter(tr[:,1], tr[:,1+M]; label="stable", markersize=4)
axislegend(ax)

ax = Axis(fig[1,2])
for k in [4, 5, 6]
    g, t = gali(ds, 1e5, k; threshold=1e-12)
    lines!(ax, log10.(t), log10.(g); label="GALI_$(k)")
end
lt = 2:0.5:5.5
lines!(ax, lt, -2(lt .- 3), label="slope -2")
lines!(ax, lt, -4(lt .- 3), label="slope -4")
lines!(ax, lt, -6(lt .- 3), label="slope -6")

xlims!(ax, 2, 5.5)
ylims!(ax, -12, 2)
fig
```

### Continuous Example
As an example of a continuous system, let's see the Henon-Heiles:
```@example MAIN
using DynamicalSystems
using CairoMakie, OrdinaryDiffEq
Δt = 1.0
diffeq = (abstol=1e-9, retol=1e-9, alg = Vern9(), maxiters = typemax(Int))
sp = [0, .295456, .407308431, 0] # stable periodic orbit: 1D torus
qp = [0, .483000, .278980390, 0] # quasiperiodic orbit: 2D torus
ch = [0, -0.25, 0.42081, 0]      # chaotic orbit
ds = Systems.henonheiles(sp)
```
Let's see what happens with a quasi-periodic orbit:
```@example MAIN
tr = trajectory(ds, 10000.0, qp; Δt, diffeq)
fig, ax = scatter(tr[:,1], tr[:,3]; label="qp", markersize=2)
axislegend(ax)

ax = Axis(fig[1,2]; yscale = log)
for k in [2,3,4]
    g, t = gali(ds, 10000.0, k; u0 = qp, Δt, diffeq)
    logt = log.(t)
    lines!(ax, logt, g; label="GALI_$(k)")
    if k == 2
        lines!(ax, logt, 1 ./ t.^(2k-4); label="slope -$(2k-4)")
    else
        lines!(ax, logt, 100 ./ t.^(2k-4); label="slope -$(2k-4)")
    end
end
ylims!(ax, 1e-12, 2)
fig
```

Finally, here is GALI of a continuous system with a chaotic orbit
```@example MAIN
tr = trajectory(ds, 10000.0, ch; Δt, diffeq)
fig, ax = scatter(tr[:,1], tr[:,3]; label="ch", markersize=2, color = (Main.COLORS[1], 0.5))
axislegend(ax)

ax = Axis(fig[1,2]; yscale = log)
ls = lyapunovspectrum(ds, 5000.0; Δt, u0 = ch, diffeq)
for k in [2,3,4]
    ex = sum(ls[1] - ls[j] for j in 2:k)
    g, t = gali(ds, 1000, k; u0 = ch, Δt = Δt, diffeq)
    lines!(t, exp.(-ex.*t); label="exp. k=$k")
    lines!(t, g; label="GALI_$(k)")
end
ylims!(ax, 1e-16, 1)
fig
```

As you can see, the results of both discrete and continuous systems match very well the theory described in [`gali`](@ref).

### Using GALI
No-one in their right mind would try to fit power-laws in order to distinguish between chaotic and regular behavior, like the above examples. These were just proofs that the method works as expected in all cases.

The most common usage of $\text{GALI}_k$ is to define a (sufficiently) small
amount of time and a (sufficiently) small threshold and see whether $\text{GALI}_k$
stays below it, for a (sufficiently) big $k$.

The following is an example of advanced usage (see [Advanced documentation](@ref)):
```julia
using DynamicalSystems, PyPlot

function main(k)
# Measure of chaoticity: final time of gali_2
dens = 201
chaoticity = zeros(Int, dens, dens)

θs = ps = range(0, stop = 2π, length = dens+1)
ds = Systems.standardmap(k = k)

tinteg = tangent_integrator(ds, 2)

for (i, θ) ∈ enumerate(θs[1:dens])
    println("i = $(i)")
    for (j, p) ∈ enumerate(ps[1:dens])

        # new initial state is the system initial state
        u0 = SVector{2}(θ, p)
        reinit!(tinteg, u0, orthonormal(2,2))

        # Low-level call signature of gali:
        #  gali(tinteg, tmax, Δt, threshold)
        chaoticity[i, j] = gali(tinteg, 500, 1, 1e-12)[2][end]
    end
end
figure()
pcolormesh(θs .- (θs[2] - θs[1])/2, ps .- (ps[2] - ps[1])/2,
chaoticity')
colorbar()
xlabel("\$\\theta\$")
ylabel("\$p\$")
return
end

main(0.9);
```
![](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/gali_standardmap.png?raw=true)

### Regular orbits in the Henon-Heiles system
In this example we use the [`poincaresos`](@ref) function to produce
surfaces of section of the [`Systems.henonheiles`](@ref) system
at different energies. At each energy [`gali`](@ref) is used to color-code
each initial condition according to how chaotic/regular it is, i.e. how much time
does it need to exceed the `threshold` of [`gali`](@ref).


```@raw html
<video width="100%" height="auto" controls loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/gali_psos_henonhelies.mp4?raw=true" type="video/mp4">
</video>
```

## Predictability of a chaotic system
Even if a system is "formally" chaotic, it can still be in phases where it is partially
predictable, because the correlation coefficient between nearby trajectories vanishes
very slowly with time.
[Wernecke, Sándor & Gros](https://www.nature.com/articles/s41598-017-01083-x) have
developed an algorithm that allows one to classify a dynamical system to one of three
categories: strongly chaotic, partially predictable chaos or regular 
(called *laminar* in their paper).

We have implemented their algorithm in the function [`predictability`](@ref).
Note that we set up the implementation to always return regular behavior for negative
Lyapunov exponent. You may want to override this for research purposes.

```@docs
predictability
```

### Example Hénon Map
We will create something similar to figure 2 of the paper, but for the Hénon map.

```@example MAIN
fig = Figure()
ax = Axis(fig[1,1]; xlabel = L"a", ylabel = L"x")
he = Systems.henon()
as = 0.8:0.01:1.225
od = orbitdiagram(he, 1, 1, as; n = 2000, Ttr = 2000)
colors = Dict(:REG => "blue", :PPC => "green", :SC => "red")
for (i, a) in enumerate(as)
    set_parameter!(he, 1, a)
    chaos_type, ν, C = predictability(he; T_max = 400000, Ttr = 2000)
    scatter!(ax, a .* ones(length(od[i])), od[i]; 
    color = (colors[chaos_type], 0.05), markersize = 2)
end
ax.title = "predictability of Hénon map"
fig
```

## The 0-1 test for chaos
The methods mentioned in this page so far require a `DynamicalSystem` instance.
But of course this is not always the case. The so-called "0 to 1" test for chaos, by 
Gottwald & Melbourne, takes as an input a timeseries and outputs a boolean `true` if 
the timeseries is chaotic or `false` if it is not.

Notice that the method does have a lot of caveats, so you should read the review paper 
before using. Also, it doesn't work for noisy data.

```@docs
testchaos01
```

## Expansion entropy
The expansion entropy is a quantity that is suggested by B. Hunt and E. Ott as a measure 
that can define chaos (so far no widely accepted definition of chaos exists). 
Positive expansion entropy means chaos.

```@docs
expansionentropy
boxregion
expansionentropy_sample
expansionentropy_batch
```
