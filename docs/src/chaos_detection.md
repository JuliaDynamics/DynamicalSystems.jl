# Chaos Detection
Can detect chaos with maximum lyapunov exponent.

There are better ways.

## Generalized Alignment Index
GALI for short is the shit.
```@docs
gali
```
### Examples
Let's see an example for a discrete:
```julia
using DynamicalSystems, PyPlot, LsqFit

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
```
![GALI Discrete](https://i.imgur.com/Kl8bFIR.png)
The left figure demonstrates the drawback of SALI≡GALI_2 (which was one of the incentives for the invention of GALI_k) : If $\lambda_1$ and $\lambda_2$ are very close, the convergence is very slow.

As an example of a continuous system, let's see the [`henonhelies`](@ref):
```julia
using DynamicalSystems, PyPlot

ds = Systems.henonhelies() # default i.c. is chaotic

tmax = 1000
dt = 0.5

tr = trajectory(ds, tmax, dt=dt)

subplot(2,1,1)
scatter(tr[:,1], tr[:,2], alpha = 0.5, s= 5, label="orbit")
legend()

subplot(2,1,2)
ls = lyapunovs(ds, 1000.0, dt=dt)

for k in [2,3,4]
    ex = sum(ls[1] - ls[j] for j in 2:k)
    g, t = gali(ds, k, tmax; dt = dt)
    semilogy(t, exp.(-ex.*t), label="exp. k=$k")
    semilogy(t, g, label="GALI_$(k)")
end
legend()
tight_layout()
```
![GALI Continuous](https://i.imgur.com/MML6i7M.png)

and the same process with a regular orbit:
```julia
using DynamicalSystems, PyPlot
figure()
ds = Systems.henonhelies([0, 0.1, 0.5, 0.0])
dt = 0.5

tr = trajectory(ds, 10000.0, dt=dt)

subplot(2,1,1)
plot(tr[:,1], tr[:,2], alpha = 0.5, label="orbit")
legend()

subplot(2,1,2)
ls = lyapunovs(ds, 1000.0, dt=dt)

for k in [2,3,4]
    ex = sum(ls[1] - ls[j] for j in 2:k)
    g, t = gali(ds, k, tmax; dt = dt)
    loglog(t, 1./t.^(2k-4), label="exp. k=$k")
    loglog(t, g, label="GALI_$(k)")
end
legend()
tight_layout()

```
