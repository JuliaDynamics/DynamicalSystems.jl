#towel map:
println("Testing discrete system evolution...")
using StaticArrays, Base.Test, DynamicalSystems

@testset "Logistic Map" begin

  d1 = Systems.logistic(0.1)
  d2 = DiscreteDS1D(0.1, d1.eom)
  d3 = DiscreteDS1D(big(0.1), d1.eom, d1.deriv)

  @testset "Evolution & Timeseries" begin
    dd1 = evolve!(d1)
    dd2 = evolve!(d2)
    dd3 = evolve!(d3)
    @test dd1.state == dd2.state
    @test dd1.state ≈ dd3.state
    @test typeof(dd3.state) == BigFloat
    ts1 = timeseries(dd1, 100, mutate = false)
    ts3 = timeseries(dd3, 100, mutate = false)
    @test ts1[10] ≈ ts3[10]
    @test eltype(ts3) == BigFloat
  end
  @testset "Derivatives" begin
    f1 = d1.deriv(d1.state)
    f2 = d2.deriv(d2.state)
    f3 = d3.deriv(d3.state)
    @test isapprox(f1, f2;rtol = 1e-12)
    @test isapprox(f1, f3;rtol = 1e-12)
    @test typeof(f3) == BigFloat
  end
end

@testset "Folded-Towel Map" begin

  s1 = Systems.towel(0.1ones(3))
  s2 = DiscreteDS(0.1ones(3), s1.eom)
  s4 = DiscreteDS(round.(big.(0.1ones(3)),3), s1.eom, s1.jacob)

  @testset "Evolution & Timeseries" begin
    st1 = evolve(s1)
    st2 = evolve(s2)
    st4 = evolve(s4)

    @test isapprox.(st1, st2; rtol = 1e-12) == trues(s1.state)
    @test isapprox.(st1, st4; rtol = 1e-12) == trues(s1.state)

    evolve!(s1); evolve!(s2); evolve!(s4);
    @test s1.state ≈ s2.state
    @test s2.state ≈ s4.state

    ts = timeseries(s1, 100; mutate = false)
    @test size(ts) == (100, 3)
    ts4 = timeseries(s4, 100; mutate = false)
    @test size(ts4) == (100, 3)
    @test eltype(ts4) == BigFloat
    @test isapprox.(ts[10, :],ts4[10, :]) == trues(3)
  end
  @testset "Jacobians" begin

    J1 = jacobian(s1)
    @test typeof(J1) <: SMatrix
    J2 = jacobian(s2)
    J4 = jacobian(s4)
    @test typeof(J4) <: SMatrix

    @test isapprox.(J1, J2; rtol = 1e-6) == trues(J1)
    @test isapprox.(J1, J4; rtol = 1e-6) == trues(J1)
    @test eltype(J4) == BigFloat
  end
end
