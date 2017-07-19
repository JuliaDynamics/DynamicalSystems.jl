using OrdinaryDiffEq, ForwardDiff
import OrdinaryDiffEq.ODEProblem

export ContinuousDS, ODEProblem

#######################################################################################
#                                     Constructors                                    #
#######################################################################################
"""
    ContinuousDS(state, eom [, jacob]) <: DynamicalSystem
`D`-dimensional continuous dynamical system (used for `D ≤ 10`).
# Fields:
* `state::SVector{D}` : Current state-vector of the system, stored in the data format
  of `StaticArray`'s `SVector`.
* `eom::F` (function) : The function that represents the system's equations of motion
  (also called vector field). The function is of the format: `eom(u) -> SVector`
  which means that given a state-vector `u` it returns an `SVector` containing the
  derivatives `du` at the current state.
* `jacob::J` (function) : A function that calculates the system's jacobian matrix,
  based on the format: `jacob(u) -> SMatrix` which means that given a state-vector
  `u` it returns an `SMatrix` containing the Jacobian at that state.
  If the `jacob` is not provided by the user, it is created with *tremendous* efficiency
  using the module `ForwardDiff`. Most of the time, for low dimensional systems, this
  Jacobian is within a few % of speed of a user-defined one.
"""
mutable struct ContinuousDS{D, T<:Real, F, J} <: DynamicalSystem
  state::SVector{D,T}
  eom::F
  jacob::J
end
# constsructor without jacobian (uses ForwardDiff)
function ContinuousDS(u0::AbstractVector, eom)
  su0 = SVector{length(u0)}(u0)
  @inline ForwardDiff_jac(x) = ForwardDiff.jacobian(eom, x)
  return ContinuousDS(su0, eom, ForwardDiff_jac)
end
function ContinuousDS(u0::AbstractVector, eom, jac)
  su0 = SVector{length(u0)}(u0)
  return ContinuousDS(su0, eom, jac)
end

dimension(::ContinuousDS{D, T, F, J})  where {D<:ANY, T<:ANY, F<:ANY, J<:ANY} = D
#######################################################################################
#                                 System Evolution                                    #
#######################################################################################
"""
```julia
ODEProblem(ds::ContinuousDS, t)
```
Return a type `ODEProblem` with the given
system information (t0 is zero).
This can be passed directly into `solve` from `DifferentialEquations`.
"""
function ODEProblem(ds::ContinuousDS, t)
  odef = (t, u) -> ds.eom(u)
  OrdinaryDiffEq.ODEProblem(odef, ds.state, (zero(t), t))
end

"""
    get_sol(prob::ODEProblem, diff_eq_kwargs::Dict = Dict())
Solve the `prob` using `solve` and return the solution.
"""
function get_sol(prob::ODEProblem, diff_eq_kwargs::Dict)
  # Check if there is a solver in the keywords:
  if haskey(diff_eq_kwargs, :solver)
    solver = diff_eq_kwargs[:solver]
    pop!(diff_eq_kwargs, :solver)
    sol = solve(prob, solver; diff_eq_kwargs..., save_everystep=false)
  else
    sol = solve(prob, Tsit5(); diff_eq_kwargs..., save_everystep=false)
  end
  return sol
end
function get_sol(prob::ODEProblem)
  sol = solve(prob, Tsit5(); save_everystep=false)
end

# See discrete.jl for the documentation string
function evolve(ds::ContinuousDS, t::Real = 1.0; diff_eq_kwargs::Dict=Dict())
  prob = ODEProblem(ds, t)
  return get_sol(prob, diff_eq_kwargs)[end]
end
function evolve(
  state::AbstractVector, ds::ContinuousDS, t::Real=1.0;
  diff_eq_kwargs::Dict=Dict())

  prob = ODEProblem(ds, t)
  prob.u0 = state
  return get_sol(prob, diff_eq_kwargs)[end]
end

function evolve!(ds::ContinuousDS, t::Real = 1.0; diff_eq_kwargs::Dict = Dict())
  ds.state = evolve(ds, t, diff_eq_kwargs = diff_eq_kwargs)
  return ds
end

"""
```julia
evolve!(prob::ODEProblem [, t::Real, diff_eq_kwargs::Dict]) -> final_state
```
Evolve the problem using the solvers of `DifferentialEquations`, update
the problem's state as the final state of the solution and return that state.

If `t` is given, the problem is evolved for that much time (else the existing
`tspan` is used). Notice that in this function, `diff_eq_kwargs` is *not* a keyword
argument.
"""
function evolve!(prob::ODEProblem, t::Real, diff_eq_kwargs::Dict)
  prob.tspan = (zero(t), t)
  if haskey(diff_eq_kwargs, :solver)
    solver = diff_eq_kwargs[:solver]
    pop!(diff_eq_kwargs, :solver)
    state = solve(prob, solver; diff_eq_kwargs..., save_everystep=false)[end]
  else
    state = solve(prob, Tsit5(); diff_eq_kwargs..., save_everystep=false)[end]
  end
  prob.u0 = state
  return state
end
function evolve!(prob::ODEProblem, diff_eq_kwargs::Dict)
  if haskey(diff_eq_kwargs, :solver)
    solver = diff_eq_kwargs[:solver]
    pop!(diff_eq_kwargs, :solver)
    state = solve(prob, solver; diff_eq_kwargs..., save_everystep=false)[end]
  else
    state = solve(prob, Tsit5(); diff_eq_kwargs..., save_everystep=false)[end]
  end
  prob.u0 = state
  return state
end
function evolve!(prob::ODEProblem)
  state = solve(prob, Tsit5(); save_everystep=false)[end]
  prob.u0 = state
  return state
end
evolve!(prob::ODEProblem, t::Real) = (prob.tspan = (zero(t), t); evolve!(prob))

function timeseries(ds::ContinuousDS, T::Real;
  dt::Real=0.05, diff_eq_kwargs::Dict=Dict(), mutate::Bool = false)

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
  sol = get_sol(prob, kw)
  TS = zeros(eltype(ds.state), length(t), D)
  for j in 1:D
    TS[:, j] .= sol[j,:]
  end
  # using: transpose(hcat(get_sol(prob, diff_eq_kwargs).u...))
  # is so absurdly tragically slower.
  if mutate
    ds.state = TS[end, :]
  end
  return TS
end

#######################################################################################
#                                 Tangent Space                                       #
#######################################################################################
function tangentbundle_setup(ds::ContinuousDS, dt)
  D = dimension(ds)

  # S is the matrix that keeps the system state in the first column
  # and tangent dynamics (Jacobian of the Flow) in the rest of the columns
  S = [ds.state eye(eltype(ds.state), D)]
  f = ds.eom
  jac = ds.jacob

  # the equations of motion `tbeom` evolve the system and the tangent dynamics
  # The e.o.m. for the tangent dynamics is simply:
  # dY/dt = J ⋅ Y
  # with J the Jacobian of the system (NOT the flow) at the current state
  function tbeom(t, u, du)
    du[:, 1] .= f(u)
    A_mul_B!(view(du, :, 2:D+1), jac(view(u, :, 1)), view(u, :, 2:D+1))
  end
  tbprob = ODEProblem(tbeom, S, (zero(dt), dt))
  return tbprob
end

#######################################################################################
#                                 Pretty-Printing                                     #
#######################################################################################
import Base.show
function Base.show(io::IO, s::ContinuousDS{N, S, F, J}) where {N<:ANY, S<:ANY, F<:ANY, J<:ANY}
  print(io, "$N-dimensional continuous dynamical system:\n",
  "state: $(s.state)\n", "e.o.m.: $F\n", "jacobian: $J")
end

@require Juno begin
  function Juno.render(i::Juno.Inline, s::ContinuousDS{N, S, F, J}) where {N<:ANY, S<:ANY, F<:ANY, J<:ANY}
    t = Juno.render(i, Juno.defaultrepr(s))
    t[:head] = Juno.render(i, Text("$N-dimensional continuous dynamical system"))
    t
  end
end
