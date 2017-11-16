if current_module() != DynamicalSystems
  using DynamicalSystems
end
using Base.Test, OrdinaryDiffEq
println("\nTesting continuous system lyapunov exponents...")

@testset "Lorenz system" begin
  ds = Systems.lorenz()
  # ds2 = ContinuousDS(ds.state, ds.eom!)
  @testset "lyapunovs" begin
    λ = lyapunovs(ds, 1e5)
    @test 0.9 < λ[1] < 1.0
    @test -0.1 < λ[2] < 0.1
    @test -15 < λ[3] < -14

    λ = lyapunovs(ds, 1e5; dt = 0.1, Ttr = 10.0,
    diff_eq_kwargs = Dict(:abstol=>1e-9, :solver => DP5()))
    @test 0.9 < λ[1] < 1.0
    @test -0.1 < λ[2] < 0.1
    @test -14.7 < λ[3] < -14.1
  end

  # @testset "lyapunovs ForwardDiff" begin
  #   λ = lyapunovs(ds2, 1e5)
  #   @test 0.89 < λ[1] < 0.91
  #   @test -0.001 < λ[2] < 0.01
  #   @test -14.6 < λ[3] < -14.5
  # end

  @testset "lyapunov" begin
    T = 10000.0
    λ1 = lyapunov(ds, T)
    @test 0.89 < λ1 < 0.92



    # test convergence:
    ls, ts = lyapunov(ds, T, Val{true}; Ttr = 100.0)
    @test size(ls, 1) == size(ls, 1)
    @test 0.89 < ls[end] < 0.92
    @test ts[end] <= T

  end
end

@testset "Roessler system" begin
  ds = Systems.roessler()
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

  @testset "lyapunov" begin
    λ1, ts = lyapunov(ds, 10000.0, Val{true}, dt =  1.0, Ttr = 10.0)
    @test 0.06 < λ1[end] < 0.08
  end
end
