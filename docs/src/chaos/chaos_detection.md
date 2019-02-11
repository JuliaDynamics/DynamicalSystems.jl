# Chaos Detection
Being able to detect and distinguish chaotic from regular behavior is crucial in the
study of dynamical systems. Most of the time a positive maximum [`lyapunov`](@ref) exponent
and a bounded system indicate chaos.

However, the convergence of the Lyapunov exponent is often very slow and
the computation costly. There are
many alternatives that are both more efficient and more accurate in characterizing
chaotic and regular motion, some of which are included in **DynamicalSystems.jl**.

!!! info "Performance depends on the solver"
    Notice that the performance of functions that use `ContinuousDynamicalSystem`s depend crucially on the chosen solver. Please see the documentation page on [Choosing a solver](@ref) for an in-depth discussion.

## Generalized Alignment Index
"GALI" for sort, is a method that relies on the fact that initially orthogonal deviation vectors tend to align towards the direction of the maximum Lyapunov exponent for chaotic
motion. It is one
of the most recent and cheapest methods for distinguishing chaos, introduced first in
2007 by Skokos, Bountis & Antonopoulos.
```@docs
gali
```
---
### Discrete Example
We will use 3 coupled standard maps as an example for a discrete system:
```@example gali
using DynamicalSystems
using PyPlot
M = 3; ks = 3ones(M); Γ = 0.1;
stable = [π, π, π, 0.01, 0, 0] .+ 0.1
chaotic = rand(2M)

ds = Systems.coupledstandardmaps(M, stable; ks=ks, Γ = Γ)
```

First, let's see the behavior of GALI for a stable orbit
```@example gali
figure(figsize = (8,4))
tr = trajectory(ds, 100000)

subplot(1,2,1)
plot(tr[:,1], tr[:,1+M], alpha = 0.5,
label="stable",marker="o", ms=1, linewidth=0)
legend()

subplot(1,2,2)
for k in [4, 5, 6]
    g, t = gali(ds, 1e5, k; threshold=1e-12)
    lt = log10.(t); lg = log10.(g)
    plot(lt, lg, label="GALI_$(k)")
end
lt = 2:0.5:5.5
plot(lt, -2(lt .- 3), label="slope -2")
plot(lt, -4(lt .- 3), label="slope -4")
plot(lt, -6(lt .- 3), label="slope -6")

xlim(2, 5.5)
ylim(-12, 2)
legend()
tight_layout()
savefig("gali_discrete_stable.png"); nothing # hide
```
![gali_discrete_stable](gali_discrete_stable.png)

Now do the same for a chaotic orbit

```@example gali
figure(figsize = (8,4))
tr = trajectory(ds, 100000, chaotic)
subplot(1,2,1)
plot(tr[:,1], tr[:,1+M], alpha = 0.5,
label="chaotic",marker="o", ms=1, linewidth=0)
legend()

subplot(1,2,2)
ls = lyapunovs(ds, 100000; u0 = chaotic)
for k in [2,3,6]
    ex = sum(ls[1] - ls[j] for j in 2:k)
    g, t = gali(ds, 1000, k; u0 = chaotic)
    semilogy(t, exp.(-ex.*t), label="exp. k=$k")
    semilogy(t, g, label="GALI_$(k)")
end
legend()
xlim(0,100)
ylim(1e-12, 1)
savefig("gali_discrete_chaos.png"); nothing # hide
```
![gali_discrete_chaos](gali_discrete_chaos.png)


### Continuous Example
As an example of a continuous system, let's see the [`henonheiles`](system_definition/#DynamicalSystems.Systems.henonheiles):
```@example gali
using DynamicalSystems
using PyPlot, OrdinaryDiffEq
sp = [0, .295456, .407308431, 0] # stable periodic orbit: 1D torus
qp = [0, .483000, .278980390, 0] # quasiperiodic orbit: 2D torus
ch = [0, -0.25, 0.42081, 0]      # chaotic orbit
ds = Systems.henonheiles(sp)
```
First, we see the behavior with a stable periodic orbit
```@example gali
figure(figsize = (8,4))
subplot(1,2,1)
dt = 1.0

diffeq = (abstol=1e-9, reltol=1e-9, alg = Tsit5(), maxiters = typemax(Int))
tr = trajectory(ds, 10000.0; dt=dt, diffeq...)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="sp",marker="o",markersize=2, linewidth=0)
legend()

subplot(1,2,2)
for k in [2,3,4]
    g, t = gali(ds, 10000.0, k; dt = dt, diffeq...)
    loglog(t, g, label="GALI_$(k)")
    if k < 4
        loglog(t, 100 ./ t.^(k-1), label="slope -$(k-1)")
    else
        loglog(t, 10000 ./ t.^(2k-4), label="slope -$(2k-4)")
    end
end
ylim(1e-12, 2)
legend();
savefig("gali_cont_stable.png"); nothing # hide
```
![gali_cont_stable](gali_cont_stable.png)

Next, let's see what happens with a quasi-periodic orbit.
Don't forget to change the `u0` arguments!
```@example gali
figure(figsize = (8,4))
subplot(1,2,1)
tr = trajectory(ds, 10000.0, qp; dt=dt, diffeq...)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="qp",marker="o",markersize=2, linewidth=0)
legend()

subplot(1,2,2)
for k in [2,3,4]
    g, t = gali(ds, 10000.0, k; u0 = qp, dt = dt, diffeq...)
    loglog(t, g, label="GALI_$(k)")
    if k == 2
        loglog(t, 1 ./ t.^(2k-4), label="slope -$(2k-4)")
    else
        loglog(t, 100 ./ t.^(2k-4), label="slope -$(2k-4)")
    end
end
ylim(1e-12, 2)
legend()
tight_layout()
savefig("gali_cont_quasi.png"); nothing # hide
```
![gali_cont_quasi](gali_cont_quasi.png)

Finally, here is GALI of a continuous system with a chaotic orbit
```@example gali
figure(figsize = (8,4))
tr = trajectory(ds, 10000.0, ch; dt=dt, diffeq...)
subplot(1,2,1)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="ch",marker="o",markersize=2, linewidth=0)
legend()

subplot(1,2,2)
ls = lyapunovs(ds, 5000.0; dt=dt, u0 = ch, diffeq...)
for k in [2,3,4]
    ex = sum(ls[1] - ls[j] for j in 2:k)
    g, t = gali(ds, 1000, k; u0 = ch, dt = dt, diffeq...)
    semilogy(t, exp.(-ex.*t), label="exp. k=$k")
    semilogy(t, g, label="GALI_$(k)")
end
legend()
ylim(1e-16, 1)
tight_layout()
savefig("gali_cont_chaos.png"); nothing # hide
```
![gali_cont_chaos](gali_cont_chaos.png)

As you can see, the results of both discrete and continuous systems match
very well the theory described in [`gali`](@ref).

## Using GALI
No-one in their right mind would try to fit power-laws in order to distinguish between chaotic and regular behavior, like the above examples. These were just demonstrations and proofs that the method works as expected in all cases.

The most common usage of $\text{GALI}_k$ is to define a (sufficiently) small
amount of time and a (sufficiently) small threshold and see whether $\text{GALI}_k$
stays below it, for a (sufficiently) big $k$.

The following is an example of [advanced usage](advanced):
```julia
using DynamicalSystems
using PyPlot, StaticArrays

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
            #  gali(tinteg, tmax, dt, threshold)
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
![](https://github.com/JuliaDynamics/JuliaDynamicsDocumentation.jl/tree/master/animations/chaos/gali_standardmap.png?raw=true)

### Regular orbits in the Henon-Heiles system
In this example we use the [`poincaresos`](@ref) function to produce
surfaces of section of the [`henonheiles`](system_definition/#DynamicalSystems.Systems.henonheiles) system
at different energies. At each energy [`gali`](@ref) is used to color-code
each initial condition according to how chaotic/regular it is, i.e. how much time
does it need to exceed the `threshold` of [`gali`](@ref).

<video width="100%" height="auto" controls>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/gali_psos_henonhelies.mp4?raw=true" type="video/mp4">
</video>

You can download the video using [this link](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/gali_psos_henonhelies.mp4?raw=true).

You can find the script that produced this animation in
`DynamicalSystems/docs/coolanimations/gali_psos_henonhelies.jl`.
