using DynamicalSystems
using Base.Test
using ParameterizedFunctions

σ = 10.0; ρ = 28.0; β = 8/3;
f = @ode_def_noinvjac LorentzParameterized begin
  dx = σp*(y-x)
  dy = x*(ρp-z) - y
  dz = x*y - βp*z
end σp=σ ρp=ρ βp=β
lorentz_jacobian_par!(J, u) = f(Val{:jac},0.0,u,J)
lorentz_eom_par!(du, u) = f(0.0, u, du)

function lorentz_eom_forw!(du, u)
  du[1] = σ*(u[2]-u[1])
  du[2] = u[1]*(ρ-u[3]) - u[2]
  du[3] = u[1]*u[2] - β*u[3]
  return du
end

@testset "Construction of ContinuousDynamicalSystem" begin
  @testset "Lorentz System" begin
    @test typeof(DynamicalSystems.lorentz(rand(3))) ==
    ContinuousDynamicalSystem
    @test typeof(DynamicalSystems.lorentz(rand(3), σ=15.0, ρ = 6//13)) ==
    ContinuousDynamicalSystem
    @test typeof(ContinuousDynamicalSystem(rand(3), lorentz_eom_forw!)) ==
    ContinuousDynamicalSystem
    @test typeof(ContinuousDynamicalSystem(rand(3), lorentz_eom_par!, lorentz_jacobian_par!)) ==
    ContinuousDynamicalSystem
  end
end
