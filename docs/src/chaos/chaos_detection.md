# Chaos Detection
Being able to detect and distinguish chaotic from regular behavior is crucial in the
study of dynamical systems. Most of the time a positive maximum [`lyapunov`](@ref) exponent
and a bounded system indicate chaos.

However, the convergence of the Lyapunov exponent is often very slow and
the computation costly. There are
many alternatives that are both more efficient and more accurate in characterizing
chaotic and regular motion, some of which are included in **DynamicalSystems.jl**.

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
```julia
using DynamicalSystems
using PyPlot; figure()
M = 3; ks = 3ones(M); Γ = 0.1;
stable = [π, π, π, 0.01, 0, 0] .+ 0.1
chaotic = rand(2M)

ds = Systems.coupledstandardmaps(M, stable; ks=ks, Γ = Γ)
tr = trajectory(ds, 100000)

subplot(2,2,1)
plot(tr[:,1], tr[:,1+M], alpha = 0.5,
label="stable",marker="o", ms=1, linewidth=0)
legend()
#
subplot(2,2,2)
for k in [2,3,4, 5, 6]
    g, t = gali(ds, 1e5, k; threshold=1e-12)
    lt = log10.(t); lg = log10.(g)

    plot(lt, lg, label="GALI_$(k)")
end
lt = 2:0.5:5.5
plot(lt, zeros(lt), label="const")
plot(lt, -2(lt - 3), label="slope -2")
plot(lt, -4(lt - 3), label="slope -4")
plot(lt, -6(lt - 3), label="slope -6")

xlim(2, 5.5)
ylim(-12, 1)
legend(fontsize=12)
tight_layout()


ds = Systems.coupledstandardmaps(M, chaotic; ks=ks, Γ = Γ)
tr = trajectory(ds, 100000)
subplot(2,2,3)
plot(tr[:,1], tr[:,1+M], alpha = 0.5,
label="chaotic",marker="o", ms=1, linewidth=0)
legend()


subplot(2,2,4)
ls = lyapunovs(ds, 100000)
for k in [2,3,4,5 ,6]
    ex = sum(ls[1] - ls[j] for j in 2:k)
    g, t = gali(ds, 1000, k)
    semilogy(t, exp.(-ex.*t), label="exp. k=$k")
    semilogy(t, g, label="GALI_$(k)")
end
legend(fontsize=12)
xlim(0,50)
ylim(1e-12, 1)

```
![GALI Discrete](https://i.imgur.com/tzoaOqV.png)


### Continuous Example
As an example of a continuous system, let's see the [`henonheiles`](system_definition/#DynamicalSystems.Systems.henonheiles):
```julia
using DynamicalSystems
using PyPlot, OrdinaryDiffEq
figure(figsize=(10, 12))
sp = [0, .295456, .407308431, 0] #stable periodic orbit: 1D torus
qp = [0, .483000, .278980390, 0] #quasiperiodic orbit: 2D torus
ch = [0, -0.25, 0.42081, 0] # chaotic orbit
dt = 1.0

subplot(3,2,1)
ds = Systems.henonheiles(sp)
diffeq = Dict(:abstol=>1e-9, :reltol=>1e-9, :solver => Tsit5(), :maxiters => typemax(Int))
tr = trajectory(ds, 10000.0, dt=dt, diff_eq_kwargs = diffeq)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="sp",marker="o",markersize=2, linewidth=0)
legend()

subplot(3,2,2)
for k in [2,3,4]
    g, t = gali(ds, 50000.0, k; dt = dt, diff_eq_kwargs = diffeq)
    if k < 4
        loglog(t, 1./t.^(k-1), label="slope -$(k-1)")
    else
        loglog(t, 1./t.^(2k-4), label="slope -$(2k-4)")
    end
    loglog(t, g, label="GALI_$(k)")
end
legend(fontsize=12)

subplot(3,2,3)
ds = Systems.henonheiles(qp)
tr = trajectory(ds, 10000.0, dt=dt, diff_eq_kwargs = diffeq)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="qp",marker="o",markersize=2, linewidth=0)
legend()

subplot(3,2,4)
for k in [2,3,4]
    g, t = gali(ds, 10000.0, k; dt = dt, diff_eq_kwargs = diffeq)
    loglog(t, 1./t.^(2k-4), label="slope -$(2k-4)")
    loglog(t, g, label="GALI_$(k)")
end
legend(fontsize=12)
tight_layout()

ds = Systems.henonheiles(ch)
diffeq = Dict(:abstol=>1e-6, :reltol=>1e-6, :solver => Tsit5(), :maxiters => typemax(Int))
tr = trajectory(ds, 10000.0, dt=dt, diff_eq_kwargs = diffeq)
subplot(3,2,5)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="ch",marker="o",markersize=2, linewidth=0)
legend()

subplot(3,2,6)
ls = lyapunovs(ds, 5000.0, dt=dt)
for k in [2,3,4]
    ex = sum(ls[1] - ls[j] for j in 2:k)
    g, t = gali(ds, 1000, k; dt = dt)
    semilogy(t, exp.(-ex.*t), label="exp. k=$k")
    semilogy(t, g, label="GALI_$(k)")
end
legend(fontsize=12)
tight_layout()
```
![GALI Continuous](https://i.imgur.com/VJE6MpC.png)
As you can see, the results of both discrete and continuous match
very well the theory described in
[`gali`](@ref).

### Using GALI
No-one in their right mind would try to fit power-laws in order to distinguish between
chaotic and regular behavior, like the above examples. These were just demonstrations and proofs that the method works as expected in all cases.

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

    θs = ps = linspace(0, 2π, dens+1)
    ds = Systems.standardmap(k = k)

    tinteg = tangent_integrator(ds, 2)

    for (i, θ) ∈ enumerate(θs[1:dens])
        println("i = $(i)")
        for (j, p) ∈ enumerate(ps[1:dens])

            # new initial state is the system initial state
            u0 = SVector{2}(θ, p)
            reinit!(tinteg, u0; Q0 = orthonormal(2,2))

            # Low-level call signature of gali:
            #  _gali(tinteg, tmax, dt, threshold)
            chaoticity[i, j] = ChaosTools._gali(tinteg, 500, 1, 1e-12)[2][end]
        end
    end

    pcolormesh(θs .- (θs[2] - θs[1])/2, ps .- (ps[2] - ps[1])/2,
    chaoticity')
    colorbar()

end

main(0.9)
```
and after about a minute you will get:
![Chaos detection](https://i.imgur.com/z85KBRh.png)
