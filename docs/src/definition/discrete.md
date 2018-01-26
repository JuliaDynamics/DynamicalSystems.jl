# Discrete Systems
Discrete systems are of the form:
```math
\vec{x}_{n+1} = \vec{f}(\vec{x}_n).
```
**DynamicalSystems.jl** categorizes discrete systems in three cases, due to the
extremely performant handling that [`StaticArrays`](https://github.com/JuliaArrays/StaticArrays.jl) offers for small dimensionalities.

Handling of discrete systems is done exclusively from **DynamicalSystems.jl** and so
there is no interaction with DifferentialEquations.jl. This also means that the
definition of a discrete system may differ slightly from a continuous one.


!!! info "Non-autonomous systems"
    To define a discrete system that depends on "time" $n$, extend the
    equations of motion by introducing a new variable $\tau$ such that
    $\tau_{n+1} =  \tau_n + 1$ and use this in the equations for the other
    variables.


## High-Dimensional
At around `D=10` dimensions, Static Arrays start to become less efficient than Julia's
base Arrays, provided that the latter use in-place operations. For cases of
discrete systems with much high dimensionality, we offer a
type called `BigDiscreteDS`:
```@docs
BigDiscreteDS
```
---
The source code of the pre-defined [coupled standard maps](definition/systems#DynamicalSystemsBase.Systems.coupledstandardmaps) can
serve as an example of a `BigDiscreteDS` definition *(we do not show it here because it is very large and very complicated*).

Just keep in mind that the equations of motion for `BigDiscreteDS` are of the
form `eom!(xnew, x, p)`!

## Low-dimensional
The definition of low-dimensional discrete systems differs fundamentally from
high dimensional ones, because everything is *much* more efficiently done with
statically sized vectors. The `struct` representing such systems is called `DiscreteDS`:
```@docs
DiscreteDS
```
---
!!! note "Return form of the `eom` function"
    It is **heavily** advised that the equations of motion `eom` function returns an `SVector` from
    the julia package [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and similarly the `jacob` function returns an `SMatrix` in the case of `DiscreteDS`.

For example, here is the case if the pre-defined [henon map](definition/systems#DynamicalSystemsBase.Systems.henon):
```julia
function henon(u0=zeros(2); a = 1.4, b = 0.3)
    henon_eom(x, p) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
    henon_jacob(x, p) = @SMatrix [-2*p[1]*x[1] 1.0; p[2] 0.0]
    return DiscreteDS(u0, henon_eom, henon_jacob; parameters = [a, b])
end # should give lyapunov exponents [0.4189, -1.6229]
```

To change the parameter `a` you would use `ds.p[1] = 123` and the change
will affect *both* the equations of motion as well as the Jacobian!

## One-Dimensional
In the case of maps, there a special structure for one-dimensional systems.
The syntax is `DiscreteDS1D(state, eom [, deriv]; parameters = nothing)`.
In this one-dimensional case, you don't need to worry about StaticArrays.jl
because everything is in plain numbers.

For example, the [logistic map](definition/systems#DynamicalSystemsBase.Systems.henon)
is defined as:
```julia
function logistic(x0=rand(); r = 4.0)
    @inline logistic_eom(x, p) = p[1]*x*(1-x)
    @inline logistic_jacob(x, p) = p[1]*(1-2x)
    return DiscreteDS1D(x0, logistic_eom, logistic_jacob; parameters = [r])
end
```

Once again, if you skip the derivative functions it will be calculated automatically
using ForwardDiff.jl.
