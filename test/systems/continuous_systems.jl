using DynamicalSystems, Base.Test
println("\nTesting continuous system evolution...")

@testset "Lorenz System" begin

  lo11 = Systems.lorenz()
  lo22 = ContinuousDS(copy(lo11.state), lo11.eom!)
  lo33 = Systems.lorenz(big.([0.0, 10.0, 0.0]))

  @testset "Construction" begin
    @test typeof(lo11) <: ContinuousDS
    @test typeof(lo33) <: ContinuousDS
  end
  @testset "Evolve & kwargs" begin
    # test evolve(system):
    st1 = evolve(lo11, 1.0)
    st3 = evolve(lo33, 1.0)
    @test st1 ≈ st3
    # Test evolve(system,keywords):
    st1 = evolve(lo11, 1.0;
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    st3 = evolve(lo33, 1.0;
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    @test st1 ≈ st3
    # Test evolve!(system, keywords):
    evolve!(lo11, 1.0;
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    evolve!(lo33, 1.0;
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    @test lo11.state ≈ lo33.state
  end

  @testset "Timeseries" begin
    # timeseries pure:
    ts1 = timeseries(lo11, 2.0)
    ts3 = timeseries(lo33, 2.0)
    @test eltype(ts3[1]) == BigFloat
    @test size(ts1) == size(ts3)
    @test ts1[end, :] ≈ ts3[end,:]
    # timeseries with diff_eq_kwargs and dt:
    ts1 = timeseries(lo11, 2.0; dt=0.1,
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    ts3 = timeseries(lo33, 2.0; dt=0.1,
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    @test size(ts1) == size(ts3)
    @test ts1[end, :] ≈ ts3[end,:]
  end

  lo11 = Systems.lorenz()
  lo22 = ContinuousDS(lo11.state, lo11.eom!)
  lo33 = Systems.lorenz(big.([0.0, 10.0, 0.0]))

  @testset "Jacobians" begin
    j1 = lo11.jacob(lo11.state);
    # j2 = lo22.jacob(lo22.state);
    j3 = lo33.jacob(lo33.state)
    @test eltype(j3) == BigFloat
    # @test j1 ≈ j2
    @test j1 ≈ j3
    s1 = evolve(lo11, 1.0)
    # s2 = evolve(lo22, 1.0)
    s3 = evolve(lo33, 1.0)
    j1 = lo11.jacob(s1);
    # j2 = lo22.jacob(s2);
    j3 = lo33.jacob(s3)
    @test eltype(j3) == BigFloat
    # @test j1 ≈ j2
    @test j1 ≈ j3
  end

end
