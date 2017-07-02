using DynamicalSystems, Base.Test, OrdinaryDiffEq

println("\nTesting continuous system lyapunov exponents...")

@testset "Lorenz system" begin
  ds = Systems.lorenz()
  ds2 = ContinuousDS(ds.state, ds.eom)
  @testset "lyapunovs" begin
    λ = lyapunovs(ds, 2e4)
    @test 0.89 < λ[1] < 0.91
    @test 0.0 < λ[2] < 0.01
    @test -14.6 < λ[3] < -14.5

    λ = lyapunovs(ds, 2e4; dt = 0.1, Ttr = 10.0,
    diff_eq_kwargs = Dict(:abstol=>1e-9, :solver => DP5()))
    @test 0.89 < λ[1] < 0.91
    @test 0.0 < λ[2] < 0.01
    @test -14.6 < λ[3] < -14.5
  end

  @testset "lyapunovs ForwardDiff" begin
    λ = lyapunovs(ds2, 2e4)
    @test 0.89 < λ[1] < 0.91
    @test 0.0 < λ[2] < 0.01
    @test -14.6 < λ[3] < -14.5
  end

  @testset "lyapunov" begin
    λ1 = lyapunov(ds, 2000, dt =  0.1)
    λ2 = lyapunov(ds2, 2000, dt = 0.1,
    diff_eq_kwargs = Dict(:solver => DP5(), :abstol => 1e-9))
    @test 0.89 < λ1[1] < 0.91
    @test 0.89 < λ2[1] < 0.91
  end
end

@testset "Roessler system" begin
  ds = Systems.roessler()
  ds2 = ContinuousDS(ds.state, ds.eom)
  @testset "lyapunovs" begin
    λ = lyapunovs(ds, 5e4)
    @test 0.06 < λ[1] < 0.08
    @test -0.01 < λ[2] < 0.01
    @test -5.4 < λ[3] < -5.39

    λ = lyapunovs(ds, 5e4; dt = 0.1, Ttr = 10.0,
    diff_eq_kwargs = Dict(:abstol=>1e-9, :solver => DP5()))
    @test 0.06 < λ[1] < 0.08
    @test -0.01 < λ[2] < 0.01
    @test -5.4 < λ[3] < -5.39
  end

  @testset "lyapunovs ForwardDiff" begin
    λ = lyapunovs(ds2, 5e4)
    @test 0.06 < λ[1] < 0.08
    @test -0.01 < λ[2] < 0.01
    @test -5.4 < λ[3] < -5.39
  end

  @testset "lyapunov" begin
    λ1 = lyapunov(ds, 5000, dt =  0.1)
    λ2 = lyapunov(ds2, 5000, dt = 0.1,
    diff_eq_kwargs = Dict(:solver => DP5(), :abstol => 1e-9))
    @test 0.069 < λ1[1] < 0.075
    @test 0.069 < λ2[1] < 0.075
  end
end
