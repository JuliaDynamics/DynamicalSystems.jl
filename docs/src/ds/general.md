# Dynamical System Definition
In **DynamicalSystems.jl** a Dynamical System can be either in continuous time
```math
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}, p, t),
```
or discrete time
```math
\vec{u}_{n+1} = \vec{f}(\vec{u}_n, p, n)
```
where $u$ is the state of the system and $p$ contains the parameters of the system.
The function $f$ is called the _dynamic rule_ of the system, also known as _equations of motion_.

In addition $f$, information about the Jacobian of the system $J_f$ is also used throughout the library.

!!! warning
    Keep in mind that almost all functions of **DynamicalSystems.jl** that use a `DynamicalSystem` assume that ``f`` is differentiable!

## Creating a Dynamical System
```@docs
DynamicalSystem
```

### Definition Table

Here is a handy table that summarizes in what form should be the functions required for the equations of motion and the Jacobian, for each system type:

|          System Type         |    equations of motion    |            Jacobian            |
|:----------------------------:|:-------------------------:|:------------------------------:|
| in-place (big systems)       | `eom!(du, u, p, t)`       | `jacobian!(J, u, p, t)`        |
| out-of-place (small systems) | `eom(u, p, t) -> SVector` | `jacobian(u, p, t) -> SMatrix` |

!!! tip "Use mutable containers for the parameters"
    It is highly suggested to use a subtype of `Array`,  [`LMArray`](https://github.com/JuliaDiffEq/LabelledArrays.jl) or a dictionary for the container
    of the model's parameters. Some functions offered by **DynamicalSystems.jl**,
    like e.g. [`orbitdiagram`](@ref),
    assume that the parameters can be first accessed by `p[x]` with `x` some qualifier
    as well as that this value can be set by `p[x] = newvalue`.

    The [Labelled Arrays](https://github.com/JuliaDiffEq/LabelledArrays.jl) package
    offers `Array` implementations that can be accessed both by index as
    well as by some name.


### Convenience functions
The following functions are defined for convenience for any dynamical system:
```@docs
dimension
jacobian
set_parameter!
```

## Example: continuous, out-of-place
Let's see an example for a small system, which is a case where out-of-place
equations of motion are preferred.
```julia
using DynamicalSystems # also exports relevant StaticArrays names
# Lorenz system
# Equations of motion:
@inline @inbounds function loop(u, p, t)
    σ = p[1]; ρ = p[2]; β = p[3]
    du1 = σ*(u[2]-u[1])
    du2 = u[1]*(ρ-u[3]) - u[2]
    du3 = u[1]*u[2] - β*u[3]
    return SVector{3}(du1, du2, du3)
end
# Jacobian:
@inline @inbounds function loop_jac(u, p, t)
    σ, ρ, β = p
    J = @SMatrix [-σ  σ  0;
    ρ - u[3]  (-1)  (-u[1]);
    u[2]   u[1]  -β]
    return J
end

ds = ContinuousDynamicalSystem(loop, rand(3), [10.0, 28.0, 8/3], loop_jac)
```
```
3-dimensional continuous dynamical system
 state:     [0.068248, 0.828095, 0.0743729]
 e.o.m.:    loop
 in-place?  false
 jacobian:  loop_jac
```

## Example: discrete, in-place
The following example is only 2-dimensional, and thus once again it is "correct" to
use out-of-place version with `SVector`. For the sake of example though, we use
the in-place version.
```julia
# Henon map.
# equations of motion:
function hiip(dx, x, p, n)
    dx[1] = 1.0 - p[1]*x[1]^2 + x[2]
    dx[2] = p[2]*x[1]
    return
end
# Jacobian:
function hiip_jac(J, x, p, n)
    J[1,1] = -2*p[1]*x[1]
    J[1,2] = 1.0
    J[2,1] = p[2]
    J[2,2] = 0.0
    return
end
ds = DiscreteDynamicalSystem(hiip, zeros(2), [1.4, 0.3], hiip_jac)
```
```
2-dimensional discrete dynamical system
 state:     [0.0, 0.0]
 e.o.m.:    hiip
 in-place?  true
 jacobian:  hiip_jac
```
Or, if you don't want to write a Jacobian and want to use the
auto-differentiation capabilities of **DynamicalSystems.jl**, which use the module
[`ForwardDiff`](http://www.juliadiff.org/ForwardDiff.jl/stable/index.html):
```julia
ds = DiscreteDynamicalSystem(hiip, zeros(2), [1.4, 0.3])
```
```
2-dimensional discrete dynamical system
 state:     [0.0, 0.0]
 e.o.m.:    hiip
 in-place?  true
 jacobian:  ForwardDiff
```

## Complex Example
In this example we will go through the implementation of the coupled standard maps
from our [Predefined Dynamical Systems](@ref). It is the most complex implementation
and takes full advantage of the flexibility of the constructors. The example will use a function-like-object as equations of motion, as well as a sparse matrix for the Jacobian.

Coupled standard maps is a big mapping that can have arbitrary number of
equations of motion, since you can couple `N` standard maps which are 2D maps, like:

```math
\theta_{i}' = \theta_i + p_{i}' \\
p_{i}' = p_i + k_i\sin(\theta_i) - \Gamma \left[\sin(\theta_{i+1} - \theta_{i}) + \sin(\theta_{i-1} - \theta_{i}) \right]
```

To model this, we will make a dedicated `struct`, which is parameterized on the
number of coupled maps:
```julia
struct CoupledStandardMaps{N}
    idxs::SVector{N, Int}
    idxsm1::SVector{N, Int}
    idxsp1::SVector{N, Int}
end
```
(what these fields are will become apparent later)

We initialize the struct with the amount of standard maps we want to couple,
and we also define appropriate parameters:
```julia
M = 5  # couple number
u0 = 0.001rand(2M) #initial state
ks = 0.9ones(M) # nonlinearity parameters
Γ = 1.0 # coupling strength
p = (ks, Γ) # parameter container

# Create struct:
SV = SVector{M, Int}
idxs = SV(1:M...) # indexes of thetas
idxsm1 = SV(circshift(idxs, +1)...)  #indexes of thetas - 1
idxsp1 = SV(circshift(idxs, -1)...)  #indexes of thetas + 1
# So that:
# x[i] ≡ θᵢ
# x[[idxsp1[i]]] ≡ θᵢ+₁
# x[[idxsm1[i]]] ≡ θᵢ-₁
csm = CoupledStandardMaps{M}(idxs, idxsm1, idxsp1);
```

We will now use this struct to define a [function-like-object](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects), a Type that also acts as a function.
```julia
function (f::CoupledStandardMaps{N})(xnew::AbstractVector, x, p, n) where {N}
    ks, Γ = p
    @inbounds for i in f.idxs

        xnew[i+N] = mod2pi(
            x[i+N] + ks[i]*sin(x[i]) -
            Γ*(sin(x[f.idxsp1[i]] - x[i]) + sin(x[f.idxsm1[i]] - x[i]))
        )

        xnew[i] = mod2pi(x[i] + xnew[i+N])
    end
    return nothing
end
```

We will use *the same* `struct` to create a function for the Jacobian:
```julia
function (f::CoupledStandardMaps{M})(
    J::AbstractMatrix, x, p, n) where {M}

    ks, Γ = p
    # x[i] ≡ θᵢ
    # x[[idxsp1[i]]] ≡ θᵢ+₁
    # x[[idxsm1[i]]] ≡ θᵢ-₁
    @inbounds for i in f.idxs
        cosθ = cos(x[i])
        cosθp= cos(x[f.idxsp1[i]] - x[i])
        cosθm= cos(x[f.idxsm1[i]] - x[i])
        J[i+M, i] = ks[i]*cosθ + Γ*(cosθp + cosθm)
        J[i+M, f.idxsm1[i]] = - Γ*cosθm
        J[i+M, f.idxsp1[i]] = - Γ*cosθp
        J[i, i] = 1 + J[i+M, i]
        J[i, f.idxsm1[i]] = J[i+M, f.idxsm1[i]]
        J[i, f.idxsp1[i]] = J[i+M, f.idxsp1[i]]
    end
    return nothing
end
```
The only reason that this is possible, is because the `eom` always takes
a `AbstractVector` as first argument, while the Jacobian always
takes an `AbstractMatrix`. Therefore we can take advantage of multiple dispatch!

Notice in addition, that the Jacobian function accesses *only half the elements of the matrix*. This is intentional, and takes advantage of the fact that the
other half is constant. We can leverage this further, by making the Jacobian a sparse matrix. Because the `DynamicalSystem` constructors allow us to give in a pre-initialized Jacobian matrix, we take advantage of that and create:
```julia
J = zeros(eltype(u0), 2M, 2M)
# Set ∂/∂p entries (they are eye(M,M))
# And they dont change they are constants
for i in idxs
    J[i, i+M] = 1
    J[i+M, i+M] = 1
end
sparseJ = sparse(J)

csm(sparseJ, u0, p, 0) # apply Jacobian to initial state
```
And finally, we are ready to create our dynamical system:
```julia
ds = DiscreteDynamicalSystem(csm, u0, p, csm, sparseJ)
```
```
10-dimensional discrete dynamical system
 state:       [0.000803001, 0.00092095, 0.000313022, …, 3.07769e-5, 0.000670152]
 e.o.m.:      CoupledStandardMaps
 in-place?    true
 jacobian:    CoupledStandardMaps
 parameters:  Tuple
```

## Automatic Jacobians
Notice that if you are using automatic differentiation for the Jacobian, you should take care to NOT define your equations of motion so that they explicitly use, or return, `Float64` numbers.
This is because `ForwardDiff` uses `DualNumbers` for differentiation.
For example, if you did
```julia
function lorenz(u,p,t)
    σ, ρ, β = p
    dx = zeros(3)
    du1 = σ*(u[2] - u[1]) +
    du2 = u[1]*(ρ - u[3]) - u[2]
    du3 = u[1]*u[2] - β*u[3]
    return SVector{Float64, 3}(du1, du2, du3)
end
```
this function could not be used to auto-differentiate, as you would get an error when adding dual numbers to `SVector{Float64}`. Instead, leave the number type untyped, or use `eltype(u)` as the number type.

## Time Evolution of Systems

**DynamicalSystems.jl** provides a convenient function for getting a trajectory
of a system at equally spaced time points:
```@docs
trajectory
```

Notice that if you want to do repeated evolutions of different states of a
continuous system, you should use the [`integrator`](@ref) interface instead.

### Example
```@example traject
using DynamicalSystems
ds = Systems.towel()
```

```@example traject
tr = trajectory(ds, 100)
```

To get every 3-rd point of the trajectory, do
```@example traject
tr = trajectory(ds, 100; dt = 3)
```

Identical syntax is used for continuous systems
```@example traject
ds = Systems.lorenz()
```
```@example traject
tr = trajectory(ds, 10.0; dt = 0.01)
```

And a final example controlling the integrator accuracy:
```@example traject
ds = Systems.lorenz()
tr = trajectory(ds, 10.0; dt = 0.1, abstol = 1e-9, reltol = 1e-9)
```

### Solution precision for continuous systems
A numerical solution of an ODE is not the "true" solution, uniquely defined by a (well-defined) ODE and an initial condition. Especially for chaotic systems, where deviations are amplified exponentially, one is left worried if the numerical solutions truly are part of the system and can truly give insight in understanding the system.

DifferentialEquations.jl offers a tool, called [Uncertainty Quantification](https://docs.juliadiffeq.org/latest/analysis/uncertainty_quantification/),
which allows users to asses up to what time-scales the numerical solution is close
to the "true" solution. For example, using the default solving parameters of
**DynamicalSystems.jl**, the Lorenz system is accurate up to time `t = 50.0`.

However, fortunately for us, there is not too much worry about the numerical solution diverging from the true solution. That is because of the [shadowing theorem](http://mathworld.wolfram.com/ShadowingTheorem.html) (or
[shadowing lemma](http://www.scholarpedia.org/article/Shadowing_lemma_for_flows)):

!!! quote "Shadowing Theorem"
    Although a numerically computed chaotic trajectory diverges exponentially from the true trajectory with the same initial coordinates, there exists an errorless trajectory with a slightly different initial condition that stays near ("shadows") the numerically computed one.

This simply means that one can always numerically study chaos not only qualitatively but also quantitatively. For more information, see the book *Chaos in Dynamical Systems* by E. Ott, or the
[scholarpedia](http://www.scholarpedia.org/article/Shadowing_lemma_for_flows) entry.
