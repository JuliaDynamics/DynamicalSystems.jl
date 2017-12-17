# Discrete Systems
Discrete systems are of the form:
```math
\vec{x}_{n+1} = \vec{f}(\vec{x}_n).
```
**DynamicalSystems.jl** categorizes discrete systems in three cases, due to the
extremely performant handling that [`StaticArrays`](https://github.com/JuliaArrays/StaticArrays.jl) offers for small dimensionalities.


## High-Dimensional
At around `D=10` dimensions, Static Arrays start to become less efficient than Julia's
base Arrays, provided that the latter use in-place operations. For cases of
discrete systems with much high dimensionality, we offer a
type called `BigDiscreteDS`:
```@docs
BigDiscreteDS
```
---
This system is identical to [`ContinuousDS`](@ref) as far as definition is concerned.
All operations are done in place, and the type is immutable. The same suggestions
about using Functors and initialized Jacobians also apply here.

See the source code of the pre-defined [coupled standard maps](definition/systems#DynamicalSystemsDef.Systems.coupledstandardmaps) for an example of a `BigDiscreteDS` definition.

---

## Low-dimensional
The definition of low-dimensional discrete systems differs fundamentally from
high dimensional ones, because everything is *much* more efficiently done with
statically sized vectors. The `struct` representing such systems is called `DiscreteDS`:
```@docs
DiscreteDS
```
---
The documentation string of the constructor is perfectly self-contained, but for the sake of clarity we will go through all the steps in the following.

`state` is simply the state the system starts (a.k.a. initial conditions) and
is always of type `SVector` from [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl).
`eom` is a *function* that takes a `state` as an input and returns the next state
as an output.
The `jacob` is also a *function* that takes a `state` as an input and returns the
Jacobian matrix of the system (at this state). So far this is actually
different than [`BigDiscreteDS`](@ref) where the functions where in-place.

!!! note "Return form of the `eom` function"
    It is **heavily** advised that the equations of motion `eom` function returns an `SVector` from
    the julia package [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and similarly the `jacob` function returns an `SMatrix`. [Numerous benchmarks](https://github.com/Datseris/DynamicalSystems.jl/tree/master/test/benchmarks) have been made in order to deduce the most efficient way to define
    a system, and this way was proved to be the best when the system's dimension is small.

Something important to note when defining a `DiscreteDS` using Functors: since
the function calls take only one argument (always a state), it is impossible to
use multiple dispatch to differentiate between a call to the e.o.m. or the Jacobian
functions.

However, it is very easy to still define both function calls
using a single `struct`, by using
a 2 argument function given to the constructor. For example:
```@example henon
using DynamicalSystems
using StaticArrays # only necessary when defining a system

mutable struct HénonMap
    a::Float64
    b::Float64
end
(f::HénonMap)(x::EomVector) = SVector{2}(1.0 - f.a*x[1]^2 + x[2], f.b*x[1])
(f::HénonMap)(x::EomVector, no::Void) = @SMatrix [-2*f.a*x[1] 1.0; f.b 0.0]

function henon(u0=zeros(2); a = 1.4, b = 0.3)
    he = HénonMap(a,b)
    # The jacobian function: (still uses the HénonMap)
    @inline jacob_henon(x) = he(x, nothing)
    return DiscreteDS(u0, he, jacob_henon)
end

ds = henon()
```
Here the example uses the type `Void` for dispatch, but you could use any other bittype
like e.g. `::Float64` and pass in `zero(T)`.

In this example case, doing `ds.eom.a = 2.5` would still affect *both* the equations
of motion as well as the Jacobian, making everything work perfectly!


### Defining a `DynamicalSystem` without Functors
As an example of defining a system without
using Functors, let's create another one of the [Predefined Systems](#predefined-systems) offered by this package, the Folded Towel map.

Because this map doesn't have any parameters, it is unnecessary to
associate a Functor with it.
```@example 2
using DynamicalSystems
using StaticArrays # only necessary when defining a system

@inline @inbounds function eom_towel(x)
x1, x2, x3 = x[1], x[2], x[3]
SVector( 3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
0.1*( (x2+0.35)*(1-2*x3) - 1 )*(1 - 1.9*x1),
3.78*x3*(1-x3)+0.2*x2 )
end

@inline @inbounds function jacob_towel(x)
    @SMatrix [3.8*(1 - 2x[1]) -0.05*(1-2x[3]) 0.1*(x[2] + 0.35);
    -0.19((x[2] + 0.35)*(1-2x[3]) - 1)  0.1*(1-2x[3])*(1-1.9x[1])  -0.2*(x[2] + 0.35)*(1-1.9x[1]);
    0.0  0.2  3.78(1-2x[3]) ]
end

u0=[0.085, -0.121, 0.075]
towel =  DiscreteDS(u0, eom_towel, jacob_towel)
```
If we did not want to write a Jacobian for it, we could do
```julia
towl_nojac = DiscreteDS(rand(3), eom_towel)
```
and the Jacobian function is created automatically.

## One-Dimensional
In the case of maps, there a special structure for one-dimensional systems.
The syntax is `DiscreteDS1D(state, eom [, deriv])`.
In this one-dimensional case, you don't need to worry about StaticArrays.jl
because everything is in plain numbers. For example:
```@example 3
using DynamicalSystems

@inline eom_logistic(r) = (x) -> r*x*(1-x)  # this is a closure
@inline deriv_logistic(r) = (x) -> r*(1-2x) # this is a closure
r = 3.7
logistic = DiscreteDS1D(rand(), eom_logistic(r), deriv_logistic(r))
```
Once again, if you skip the derivative functions it will be calculated automatically
using ForwardDiff.jl.
