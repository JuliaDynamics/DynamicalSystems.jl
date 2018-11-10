# Predefined Systems
Predefined systems exist in the `Systems` submodule in the form of functions that return a `DynamicalSystem`. They are accessed
like:
```julia
using DynamicalSystems # or DynamicalSystemsBase
ds = Systems.lorenz(ρ = 32.0)
```

So far, the predefined systems that exist in the `Systems` sub-module are:
```@autodocs
Modules = [Systems]
Order   = [:function]
```
