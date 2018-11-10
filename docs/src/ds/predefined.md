# Predefined Systems
Predefined systems exist in the `Systems` submodule exported by DynamicalSystemsBase.jl, in the form of functions that return a `DynamicalSystem`. They are accessed
like:
```julia
using DynamicalSystems
ds = Systems.lorenz(ρ = 32.0)
ts = trajectory(ds, 10.0)
```

So far, the predefined systems that exist in the `Systems` sub-module are:
```@autodocs
Modules = [Systems]
Order   = [:function]
```
