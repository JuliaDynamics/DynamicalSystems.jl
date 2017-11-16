if current_module() != DynamicalSystems
  using DynamicalSystems
end
using Base.Test
println("\nTesting discrete system lyapunov exponents...")


@testset "Towel Map" begin
  @testset "lyapunovs" begin
    ds = Systems.towel()
    λ1 = lyapunovs(ds, 1e5)
    @test 0.42 < λ1[1] < 0.44
    @test 0.36 < λ1[2] < 0.38
    @test -3.4 < λ1[3] < -3.2
  end

  @testset "lyapunovs ForwardDiff" begin
    ds = Systems.towel()
    ds = DiscreteDS(ds.state, ds.eom)
    λ1 = lyapunovs(ds, 1e5)
    @test 0.42 < λ1[1] < 0.44
    @test 0.36 < λ1[2] < 0.38
    @test -3.4 < λ1[3] < -3.1
  end

  @testset "lyapunov" begin
    ds1 = Systems.towel()
    ds2 = DiscreteDS(ds1.state, ds1.eom)
    λ1 = lyapunov(ds1, 100000)
    λ2 = lyapunov(ds2, 100000)
    @test 0.42 < λ1[1] < 0.44
    @test 0.42 < λ2[1] < 0.44
    @test isapprox(λ1, λ2; rtol = 1e-3)
  end
end

@testset "Logistic Map" begin
  lg1 = Systems.logistic(0.1, r=4)
  lg2 = DiscreteDS1D(lg1.state, lg1.eom)
  lg3 = Systems.logistic(big(0.1), r=4)
  @test typeof(lg3.state) == BigFloat
  λ1 = lyapunov(lg1, 10000000; Ttr = 100)
  λ2 = lyapunovs(lg2, 10000000; Ttr = 100)
  λ3 = lyapunovs(lg3, 100000; Ttr = 100)
  @test typeof(λ3) == BigFloat
  @test isapprox(λ1, log(2); rtol = 1e-3)
  @test isapprox(λ2, log(2); rtol = 1e-3)
  @test isapprox(λ3, log(2); rtol = 1e-3)
end

@testset "Henon Map" begin
  @testset "lyapunovs" begin
    ds = Systems.henon()
    λ1 = lyapunovs(ds, 1e6)
    @test 0.418 < λ1[1] < 0.422
    @test -1.63 < λ1[2] < -1.61
  end

  @testset "lyapunovs ForwardDiff" begin
    ds = Systems.henon()
    ds = DiscreteDS(ds.state, ds.eom)
    λ1 = lyapunovs(ds, 1e6)
    @test 0.418 < λ1[1] < 0.422
    @test -1.63 < λ1[2] < -1.61
  end

  @testset "lyapunov" begin
    ds1 = Systems.henon()
    ds2 = DiscreteDS(ds1.state, ds1.eom)
    λ1 = lyapunov(ds1, 1000000)
    λ2 = lyapunov(ds2, 1000000)
    @test 0.418 < λ1[1] < 0.422
    @test 0.418 < λ2[1] < 0.422
    @test isapprox(λ1, λ2; rtol = 1e-3)

    ls, ts = lyapunov(ds1, 1000000, Val{true})
    @test 0.418 < ls[end] < 0.422
    @test length(ts) == length(ls)
  end
end

@testset "Coupled Standard Maps" begin
  M = 5;
  ds = Systems.coupledstandardmaps(5)
  @testset "lyapunovs" begin
    ls = lyapunovs(ds, 10000)
    for i in 1:M
      @test abs(ls[i] + ls[2M - i + 1]) < 0.01
    end
  end
  @testset "lyapunov" begin
    l = lyapunov(ds, 10000)
    @test l > 0
    ll, tt = lyapunov(ds, 10000, Val{true})
    @test ll[end] ≈ l
  end
end
