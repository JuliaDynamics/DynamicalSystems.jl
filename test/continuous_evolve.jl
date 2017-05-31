using DynamicalSystems
using Base.Test

@testset "Evolution of ContinuousDynamicalSystem" begin
  @testset "Lorentz system evolution" begin
    system = DynamicalSystems.lorentz(rand(3))
    @testset "initialization" begin
      update!(system, ones(3))
    end
    @testset "evolution" begin
      evolve!(system, 1.0); evolve!(system, 1.0); evolve!(system, 1.0)
    end
    @testset "accuracy of evolution" begin
      system = DynamicalSystems.lorentz(10ones(3))
      update!(system, 10ones(3))
      x = Float64[]; y = Float64[]; z = Float64[]
      for i in 1:100
        evolve!(system, 0.01, Dict(:abstol => 1e-9, :reltol => 1e-9))
        push!(x, system.u[1])
        push!(y, system.u[2])
        push!(z, system.u[3])
      end
      update!(system, 10ones(3))
      x2 = Float64[]; y2 = Float64[]; z2 = Float64[]
      for i in 1:100
        evolve!(system, 0.01, Dict(:abstol => 1e-9, :reltol => 1e-9))
        push!(x2, system.u[1])
        push!(y2, system.u[2])
        push!(z2, system.u[3])
      end
      @test !in(false, isapprox.(x, x2))
      @test !in(false, isapprox.(y, y2))
      @test !in(false, isapprox.(z, z2))
    end
  end
end
