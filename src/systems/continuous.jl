using OrdinaryDiffEq, Requires, ForwardDiff
import OrdinaryDiffEq.ODEProblem
import OrdinaryDiffEq.ODEIntegrator

export ContinuousDS, ODEProblem, ODEIntegrator

#######################################################################################
#                                     Constructors                                    #
#######################################################################################
"""
    ContinuousDS(state, eom [, jacob]) <: DynamicalSystem
`D`-dimensional continuous dynamical system.
## Fields:
* `state::SVector{D}` : Current state-vector of the system, stored in the data format
  of `StaticArray`'s `SVector`.
* `eom` (function) : The function that represents the system's equations of motion
  (also called vector field). The function is of the format: `eom(u) -> SVector`
  which means that given a state-vector `u` it returns an `SVector` containing the
  next state.
* `jacob` (function) : The function that represents the Jacobian of the system,
  given in the format: `jacob(u) => J` (i.e. returns a matrix). If the matrix is
  an `SMatrix` from the module `StaticArrays` there are major performance gains.

If the `jacob` is not provided by the user, it is created efficiently
using the module `ForwardDiff`.
"""
mutable struct ContinuousDS{D, T<:Number, F, J} #<: DynamicalSystem
  state::SVector{D, T}
  eom::F
  jacob::J
end
function ContinuousDS(u0::AbstractVector, eom, jac)
  su0 = SVector{length(u0)}(u0)
  return ContinuousDS(su0, eom, jac)
end
# Constructor without Jacobian (nothing in the field)
function ContinuousDS(state, eom)
  jc = ForwardDiff.JacobianConfig(eom, state)
  jac = (u) -> ForwardDiff.jacobian(eom, u, jc)
  return ContinuousDS(state, eom, jac)
end

dimension{D, T, F, J}(ds::ContinuousDS{D, T, F, J}) = D
#######################################################################################
#                         Interface to DifferentialEquations                          #
#######################################################################################
"""
```julia
ODEProblem(ds::ContinuousDS, tfinal)
```
Return an `ODEProblem` with the given
system information (t0 is zero).
This can be passed directly into `solve` from `DifferentialEquations.jl`.
"""
function ODEProblem(ds::ContinuousDS, t)
  odef = (t, u) -> ds.eom(u)
  OrdinaryDiffEq.ODEProblem(odef, ds.state, (zero(t), t))
end

"""
```julia
ODEIntegrator(ds::ContinuousDS, t; diff_eq_kwargs)
```
Return an `ODEIntegrator`, by first creating an `ODEProblem(ds, t)`.
This can be used directly with the interfaces of `DifferentialEquations`.

`diff_eq_kwargs = Dict()` is a dictionary `Dict{Symbol, ANY}`
of keyword arguments
passed into the `init` of the `DifferentialEquations` package,
for example `Dict(:abstol => 1e-9)`. If you want to specify a solver,
do so by using the symbol `:solver`, e.g.:
`Dict(:solver => DP5(), :tstops => 0:0.01:t)`. This requires you to have been first
`using OrdinaryDiffEq` to access the solvers.
"""
function ODEIntegrator(ds::ContinuousDS, t; diff_eq_kwargs = Dict())
  prob = ODEProblem(ds, t)
  # Check if there is a solver in the keywords:
  if haskey(diff_eq_kwargs, :solver)
    solver = diff_eq_kwargs[:solver]
    pop!(diff_eq_kwargs, :solver)
    integrator = init(prob, solver; diff_eq_kwargs...,
    save_first=false, save_everystep=false)
  else
    integrator = init(prob, Tsit5(); diff_eq_kwargs...,
    save_first=false, save_everystep=false)
  end
  return integrator
end

"""
    get_sol(prob::ODEProblem, diff_eq_kwargs::Dict = Dict())
Solve the `prob` using `solve` and return the solution.
"""
function get_sol(prob::ODEProblem, diff_eq_kwargs::Dict = Dict())
  # Check if there is a solver in the keywords:
  if haskey(diff_eq_kwargs, :solver)
    solver = diff_eq_kwargs[:solver]
    pop!(diff_eq_kwargs, :solver)
    sol = solve(prob, solver; diff_eq_kwargs..., save_everystep=false)
  else
    sol = solve(prob, Tsit5(); diff_eq_kwargs..., save_everystep=false)
  end
  return sol.u
end

#######################################################################################
#                                Evolution of System                                  #
#######################################################################################
# See discrete.jl for the documentation string
function evolve(ds::ContinuousDS, t::Real = 1.0; diff_eq_kwargs = Dict())
  prob = ODEProblem(ds, t)
  return get_sol(prob, diff_eq_kwargs)[end]
end
function evolve!(ds::ContinuousDS, t::Real = 1.0; diff_eq_kwargs = Dict())
  ds.state = evolve(ds, t, diff_eq_kwargs = diff_eq_kwargs)
  return ds.state
end

# See discrete.jl for the documentation string
function timeseries(ds::ContinuousDS, T::Real;
  dt::Real=0.05, diff_eq_kwargs = Dict())

  # Necessary due to DifferentialEquations:
  if !issubtype(typeof(T), AbstractFloat)
    T = convert(Float64, T)
  end
  T<=0 && throw(ArgumentError("Total time `T` must be positive."))

  D = dimension(ds)
  t = zero(T):dt:T #time vector
  prob = ODEProblem(ds, T)
  kw = Dict{Symbol, Any}(diff_eq_kwargs) #nessesary conversion to add :saveat
  kw[:saveat] = t
  return get_sol(prob, kw)
end

#######################################################################################
#                                 Pretty-Printing                                     #
#######################################################################################
import Base.show
function Base.show{D, T, F, J}(io::IO, s::ContinuousDS{D, T, F, J})
  print(io, "$D-dimensional continuous dynamical system:\n",
  " state: $(s.state)\n", " e.o.m.: $F\n", " jacobian: $J")
end

@require Juno begin
  function Juno.render{D, T, F, J}(i::Juno.Inline, s::ContinuousDS{D, T, F, J})
    t = Juno.render(i, Juno.defaultrepr(s))
    t[:head] = Juno.render(i, Text("$D-dimensional continuous dynamical system"))
    # t[:children][1][:child] = Juno.Text("SVector($(s.state))")
    t
  end
end
