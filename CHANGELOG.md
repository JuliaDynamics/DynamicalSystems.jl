The package DynamicalSystems.jl does not have any actual code.
It re-exports other packages. The changelog here only lists major changes to the overarching DynamicalSystems.jl ecosystem.

# DynamicalSystems.jl v3.0

Super duper major update.

1. [Modularization](#modularization)
2. [Re-write of `DynamicalSystem`](#re-write-of-dynamicalsystem)
3. [Default integrator: `Tsit5`](#default-integrator-tsit5)
4. [New package StateSpaceSets.jl](#new-package-statespacesetsjl)

## Modularization

The packages that compose DynamicalSystems.jl have been split up into more packages, making a more modular library, and also allowing a much smaller dependency tree for users that want only a specific part of the functionality

## Re-write of `DynamicalSystem`

The inner library DynamicalSystemsBase.jl has had a complete re-write. Furthermore, now `DynamicalSystem` defines a proper, extendable interface. The set of functions that the interface satisfies are listed in `DynamicalSystem` documentation string. This has two major impacts:

1. Concept of what is a "dynamical system" is much more general. More variants of a dynamical system are offered. Any arbitrary steppable thing, like an agent based model, can become a dynamical system.

2. Now **_all concrete implementations of `DynamicalSystem` can be iteratively evolved in time via the `step!` function_**. For example, `ContinuousDynamicalSystem` is now a light wrapper around `ODEIntegrator`. After working with this library for years I realized that the "middle man" of representing a dynamical system or an ode problem was not really necessary in downstream functions, since all of them exclusively use an integrator.


## Default integrator: `Tsit5`


As a result it has no real changelog and one should see the individual packages' changelogs. In

# Major release notes
For reference, version listing of DynamicalSystems.jl v3.0 is at the end.

## New package StateSpaceSets.jl
It has all functionality surrounding `Dataset` and similar types. It therefore detaches the infastructure for handing numeric data from delay embeddings.
