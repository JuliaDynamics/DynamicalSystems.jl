using DynamicalSystems, OrdinaryDiffEq, StaticArrays, BenchmarkTools, Requires, JLD
#######################################################################################
#                                 Tangent Space                                       #
#######################################################################################
mutable struct TBJacobian{J, T<:Real}
  jac!::J
  J::AbstractArray{T}
end

function tangentbundle_setup_eom(ds::ContinuousDS, dt;
  diff_eq_kwargs=Dict(), jacob_inplace = true)

  D = dimension(ds)
  f! = ds.eom!
  jac! = ds.jacob

  # the equations of motion `tbeom` evolve the system and the tangent dynamics
  # The e.o.m. for the system is f!(t, u , du).
  # The e.o.m. for the tangent dynamics is simply:
  # dY/dt = J(u) ⋅ Y
  # with J the Jacobian of the system (NOT the flow), at the current state

  # Now there are 2 versions, because I want different approach when Jacobian is
  # in-place (it becomes a Functor).
  # If Jacobian is not in-place then it creates new array each time

  if jacob_inplace
    JAC = TBJacobian(jac!, eye(eltype(ds.state), D))
    tbeom!_inplaceJAC = (t, u, du) -> begin
      f!(view(du, :, 1), u)
      JAC.jac!(JAC.J, u)
      A_mul_B!(
      view(du, :, 2:D+1),
      JAC.J,
      view(u, :, 2:D+1)
      )
    end
    return tbeom!_inplaceJAC
  else
    tbeom! = (t, u, du) -> begin
      f!(view(du, :, 1), u)
      A_mul_B!(
      view(du, :, 2:D+1),
      jac!(view(u, :, 1)),
      view(u, :, 2:D+1)
      )
    end
    return tbeom!
  end
end

function tangentbundle_setup_integrator(ds::ContinuousDS, dt;
  diff_eq_kwargs=Dict(), jacob_inplace = true)

  tbeom! = tangentbundle_setup_eom(ds::ContinuousDS, dt;
  diff_eq_kwargs=diff_eq_kwargs, jacob_inplace = jacob_inplace)

  # S is the matrix that keeps the system state in the first column
  # and tangent dynamics (Jacobian of the Flow) in the rest of the columns
  D = dimension(ds)
  S = [ds.state eye(eltype(ds.state), D)]

  tbprob = ODEProblem(tbeom!, S, (zero(dt), dt))
  if haskey(diff_eq_kwargs, :solver)
    solver = diff_eq_kwargs[:solver]
    pop!(diff_eq_kwargs, :solver)
    tb_integ = init(tbprob, solver; diff_eq_kwargs..., save_everystep=false)
  else
    tb_integ = init(tbprob, Tsit5(); diff_eq_kwargs..., save_everystep=false)
  end
  return tb_integ
end
#######################################################################################
#         TESTING
#######################################################################################
function lor!(du, u)
  σ = 10.0; ρ = 28.0; β = 8/3
  du[1] = σ*(u[2]-u[1])
  du[2] = u[1]*(ρ-u[3]) - u[2]
  du[3] = u[1]*u[2] - β*u[3]
end
using StaticArrays
@inline @inbounds function jacob_ret(u)
  σ = 10.0; ρ = 28.0; β = 8/3
  i = one(eltype(u))
  o = zero(eltype(u))
  @SMatrix [-σ*i           σ*i    zero(i);
            (ρ*i - u[3])   (-i)   (-u[1]);
            u[2]           u[1]   (-β*i) ]
end
@inline @inbounds function jacob_inp(J, u)
  σ = 10.0; ρ = 28.0; β = 8/3
  i = one(eltype(u))
  o = zero(eltype(u))
  J[1,1] = -σ*i; J[1,2] = σ*i; J[1,3] = zero(i)
  J[2,1] = ρ*i - u[3]; J[2,2] = -i; J[2,3] = -u[1]
  J[3,1] = u[2]; J[3,2] = u[1]; J[3,3]= -β*i
end

ds_ret = ContinuousDS([0, 10.0, 0], lor!, jacob_ret)
ds_inp = ContinuousDS([0, 10.0, 0], lor!, jacob_inp)


###################
# Lyapunov
###################

function lyapunovs(ds::ContinuousDS, N::Real;
  Ttr::Real = 1.0, diff_eq_kwargs::Dict = Dict(), dt::Real = 1.0)

  tstops = dt:dt:N*dt
  D = dimension(ds)

  integ = tangentbundle_setup_integrator(ds, tstops[end]; jacob_inplace = false)
  #integ = tangentbundle_setup_integrator(ds, tstops[end]; jacob_inplace = true)

  integ.opts.advance_to_tstop=true

  λ = zeros(eltype(integ.u), D)
  Q = eye(eltype(integ.u), D)

  # Main algorithm
  for τ in tstops
    integ.u[:, 2:end] .= Q # update tangent dynamics state
    push!(integ.opts.tstops, τ)
    step!(integ)

    # Perform QR (on the tangent flow):
    Q, R = qr(view(integ.u, :, 2:D+1))
    # Add correct (positive) numbers to Lyapunov spectrum
    for j in 1:D
      λ[j] += log(abs(R[j,j]))
    end
  end
  λ./(N*dt) #return spectrum
end

function lyapunovs_inp(ds::ContinuousDS, N::Real;
  Ttr::Real = 1.0, diff_eq_kwargs::Dict = Dict(), dt::Real = 1.0)

  tstops = dt:dt:N*dt
  D = dimension(ds)

  #integ = tangentbundle_setup_integrator(ds, tstops[end]; jacob_inplace = false)
  integ = tangentbundle_setup_integrator(ds, tstops[end]; jacob_inplace = true)

  integ.opts.advance_to_tstop=true

  λ = zeros(eltype(integ.u), D)
  Q = eye(eltype(integ.u), D)

  # Main algorithm
  for τ in tstops
    integ.u[:, 2:end] .= Q # update tangent dynamics state
    push!(integ.opts.tstops, τ)
    step!(integ)

    # Perform QR (on the tangent flow):
    Q, R = qr(view(integ.u, :, 2:D+1))
    # Add correct (positive) numbers to Lyapunov spectrum
    for j in 1:D
      λ[j] += log(abs(R[j,j]))
    end
  end
  λ./(N*dt) #return spectrum
end

println(lyapunovs(ds_ret, 10000; dt = 0.1))
@btime lyapunovs(ds_ret, 10000; dt = 0.1)

println(lyapunovs_inp(ds_inp, 10000; dt = 0.1))
@btime lyapunovs_inp(ds_inp, 10000; dt = 0.1)
