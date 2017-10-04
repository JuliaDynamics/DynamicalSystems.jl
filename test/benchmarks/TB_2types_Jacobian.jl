using OrdinaryDiffEq, StaticArrays, BenchmarkTools

const σ, ρ, β = 10.0, 28.0, 8/3

mutable struct system1{T, F, J}
  state::T
  eom!::F
  jacob::J
end

println("---Benchmarking in-place version of eom---")
# in-place approach:
@inline @inbounds function eom_lorenz!(du, u)
  du[1] = σ*(u[2]-u[1])
  du[2] = u[1]*(ρ-u[3]) - u[2]
  du[3] = u[1]*u[2] - β*u[3]
end
@inline @inbounds function jacob_lorenz(u)
  i = one(eltype(u))
  o = zero(eltype(u))
  @SMatrix [-σ*i           σ*i    zero(i);
            (ρ*i - u[3])   (-i)   (-u[1]);
            u[2]           u[1]   (-β*i) ]
end# should give exponents [0.9056, 0, -14.5723]
u0 = [0, 10.0, 0]
lor1 = system1(u0, eom_lorenz!, jacob_lorenz)

function tangentbundle_setup_integrator(ds::system1, t_final)
  D = 3
  f! = ds.eom!
  jac = ds.jacob

  tbeom! = (t, u, du) -> begin
    f!(view(du, :, 1), u)
    A_mul_B!(
      view(du, :, 2:D+1),
      jac(view(u, :, 1)),
      view(u, :, 2:D+1)
    )
  end

  S = [ds.state eye(eltype(ds.state), D)]
  tbprob = ODEProblem(tbeom!, S, (zero(t_final), t_final))
  tb_integ = init(tbprob, Tsit5(); save_everystep=false)
  return tb_integ
end

integ1 = tangentbundle_setup_integrator(lor1, 100.0)
integ1.opts.advance_to_tstop=true

t = 1.0:1.0:100.0
#compile:
push!(integ1.opts.tstops, t[1]); step!(integ1)
state1 = integ1.u
#time:
for i in 2:6
  println("time for one step:")
  push!(integ1.opts.tstops, t[i]);
  @time step!(integ1);
end
#time for solve:
println("time for solve:")
@time solve!(integ1)


println("---Benchmarking SMatrix version of eom---")
# SVector approach:
mutable struct system2{T, F, J}
  state::T
  eom::F
  jacob::J
end

@inline @inbounds function eom_lorenz(u)
  @SVector [σ*(u[2]-u[1]), u[1]*(ρ-u[3]) - u[2], u[1]*u[2] - β*u[3]]
end

u0 = @SVector [0, 10.0, 0]
lor2 = system2(u0, eom_lorenz, jacob_lorenz)

function tangentbundle_setup_integrator(ds::system2, t_final)
  D = 3
  f = ds.eom
  jac = ds.jacob

  # function tbeom2(t, u)
  #   return SMatrix{3, 4}(f(u)..., (jac(view(u, :, 1))*view(u, :, 2:4))...)
  # end
  #
  # # Different approach that makes an SMatrix first for multiplication:
  # function tbeom3(t, u)
  #   mm = SMatrix{3,3}(u[:, 2:4])
  #   return SMatrix{3, 4}(f(u)..., jac(u[:, 1])*mm...)
  # end

  # Approach 3 without splatting
  function tbeom4(t, u)
    mm = SMatrix{3,3}(u[:, 2:4])
    return hcat(f(u), jac(u[:, 1])*mm)
  end

  S = SMatrix{3, 4}(ds.state..., eye(eltype(ds.state), 3)...)
  tbprob = ODEProblem(tbeom4, S, (zero(t_final), t_final))
  tb_integ = init(tbprob, Tsit5(); save_everystep=false)
  return tb_integ
end

integ2 = tangentbundle_setup_integrator(lor2, 100.0)
integ2.opts.advance_to_tstop=true
# compile:
push!(integ2.opts.tstops, t[1]); step!(integ2)
state2 = integ2.u
#time:
for i in 2:6
  println("time for one step:")
  push!(integ2.opts.tstops, t[i]);
  @time step!(integ2);
end
#time for solve:
println("time for solve:")
@time solve!(integ2)
