# Chaos Detection
Can detect chaos with maximum lyapunov exponent.

There are better ways.

## Generalized Alignment Index
GALI for short is the shit.
```@docs
gali
```
### Discrete Example
Let's see an example for a discrete:
```julia
using DynamicalSystems
using PyPlot

ds = Systems.towel()
ls = lyapunovs(ds, 10000)

figure(figsize=(10,5))
for (i, k) in enumerate([2,3])
    ax = subplot(1,2,i)
    ex = sum(ls[1] - ls[j] for j in 2:k)

    g, t = gali(ds, k, 1000; threshold = 1e-12)

    semilogy(t, g, label="GALI\$ _$(k)(t)\$")
    semilogy(t, exp.(-ex.*t), label="exp. k=$k")
    xlabel("t")
    legend()
end
tight_layout()
```
![GALI Discrete](https://i.imgur.com/Kl8bFIR.png)
The left figure demonstrates the drawback of $\text{SALI} \equiv \text{GALI}_2$ (which was one of the incentives for the invention of $\text{GALI}_k$) : If $\lambda_1$ and $\lambda_2$ are very close, the convergence is very slow.

### Continuous Example
As an example of a continuous system, let's see the [`henonhelies`](@ref):
```julia
using DynamicalSystems
using PyPlot
figure(figsize=(10, 12))
sp = [0, .295456, .407308431, 0] #stable periodic orbit: 1D torus
qp = [0, .483000, .278980390, 0] #quasiperiodic orbit: 2D torus
ch = [0, -0.25, 0.42081, 0] # chaotic orbit
dt = 1.0

subplot(3,2,1)
ds = Systems.henonhelies(sp)
diffeq = Dict(:abstol=>1e-9, :reltol=>1e-9, :solver => Tsit5())
tr = trajectory(ds, 10000.0, dt=dt, diff_eq_kwargs = diffeq)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="sp",marker="o",markersize=2, linewidth=0)
legend()

subplot(3,2,2)
for k in [2,3,4]
    g, t = gali(ds, k, 50000.0; dt = dt, diff_eq_kwargs = diffeq, threshold=1e-12)
    if k < 4
        loglog(t, 1./t.^(k-1), label="slope -$(k-1)")
    else
        loglog(t, 1./t.^(2k-4), label="slope -$(2k-4)")
    end
    loglog(t, g, label="GALI_$(k)")
end
legend(fontsize=12)

subplot(3,2,3)
ds = Systems.henonhelies(qp)
diffeq = Dict(:abstol=>1e-9, :reltol=>1e-9, :solver => Tsit5())
tr = trajectory(ds, 10000.0, dt=dt, diff_eq_kwargs = diffeq)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="qp",marker="o",markersize=2, linewidth=0)
legend()

subplot(3,2,4)
for k in [2,3,4]
    g, t = gali(ds, k, 10000.0; dt = dt, diff_eq_kwargs = diffeq, threshold=1e-12)
    loglog(t, 1./t.^(2k-4), label="slope -$(2k-4)")
    loglog(t, g, label="GALI_$(k)")
end
legend(fontsize=12)
tight_layout()

ds = Systems.henonhelies(ch)
diffeq = Dict(:abstol=>1e-6, :reltol=>1e-6, :solver => Tsit5())
tr = trajectory(ds, 50000.0, dt=dt, diff_eq_kwargs = diffeq)
subplot(3,2,5)
plot(tr[:,1], tr[:,3], alpha = 0.5,
label="ch",marker="o",markersize=2, linewidth=0)
legend()

subplot(3,2,6)
ls = lyapunovs(ds, 5000.0, dt=dt)
for k in [2,3,4]
    ex = sum(ls[1] - ls[j] for j in 2:k)
    g, t = gali(ds, k, 1000; dt = dt)
    semilogy(t, exp.(-ex.*t), label="exp. k=$k")
    semilogy(t, g, label="GALI_$(k)")
end
legend(fontsize=12)
tight_layout()
```
![GALI Continuous](https://i.imgur.com/VJE6MpC.png)
As you can see, the results match almost perfectly the theory described in
[`gali`](@ref0).

### Using GALI
No-one in their right mind would try to fit power-laws in order to distinguish between
chaotic and regular behavior, like the above examples. These were just demonstrations.

The most common usage of $\text{GALI}_k$ is to define a (sufficiently) small
amount of time and a (sufficiently) small threshold and see whether $\text{GALI}_k$
stays below it, for a (sufficiently) big $k$.

For example one could do
```julia
using DynamicalSystems
using PyPlot

# Default threshold
t = 500.0
ischaotic(ds) = gali(ds, 2, t)[2][end] < t ? 1.0 : 0.0

dens = 201
chaoticity = zeros(dens,dens)
θs = ps = linspace(0, 2π, dens+1)

for (i, θ) ∈ enumerate(θs[1:dens])
    for (j, p) ∈ enumerate(ps[1:dens])
        ds = Systems.standardmap([θ, p])
        chaoticity[i, j] = ischaotic(ds)
    end
end

pcolormesh(θs .- (θs[2] - θs[1])/2, ps .- (ps[2] - ps[1])/2,
chaoticity')
```
and after around `0.0005*201*201 ≈ 20` seconds you will get what looks
to be identical
with the phase space of the standard map:
![Chaos detection](https://i.imgur.com/lc1WfLF.png)
