#σ = 10.0; ρ = 28.0; β = 8/3;
function lorentz(u0; σ = 10.0, ρ = 28.0, β = 8/3)
  function eom!(du, u)
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
    return du
  end
  function jacob!(J, u)
    i = one(eltype(u))
    J[1,1] = -σ*i; J[1,2] = σ*i; J[1,3] = zero(i)
    J[2,1] = ρ*i - u[3]; J[2,2] = -i; J[2,3] = -u[1]
    J[3,1] = u[2]; J[3,2] = u[1]; J[3,3] = -β*i
  end# should give exponents [0.9056, 0, -14.5723]
  return ContinuousDynamicalSystem(u0, eom!, jacob!)
end
