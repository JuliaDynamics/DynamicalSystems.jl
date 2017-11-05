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

    model(x,p)= @. exp(-p[1]*x)
    fite = curve_fit(model, t, g, [ex]).param[1]



    semilogy(t, g, label="GALI\$ _$(k)(t)\$")
    semilogy(t, exp.(-ex.*t), label="exp. k=$k")
    xlabel("t")
    legend()
end
```
The left figure demonstrates the drawback of SALI≡GALI_2 (which was one of the incentives for the invention of GALI_k) : If $\lambda_1$ and $\lambda_2$ are very close, the convergence is very slow.

Example of a continuous system:
