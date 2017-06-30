using DynamicalSystems
using Base.Test

println("\nTesting continuous system evolution...")

@testset "Lorenz System" begin

  lo11 = Systems.lorenz()
  lo22 = ContinuousDS([0.0, 10.0, 0.0], lo11.eom)
  lo33 = Systems.lorenz(big.([0.0, 10.0, 0.0]))
  @testset "Construction" begin
    @test typeof(lo11) <: ContinuousDS
    @test typeof(lo22) <: ContinuousDS
    @test typeof(lo33) <: ContinuousDS
  end
  @testset "Evolve & kwargs" begin
    # test evolve(state):
    st1 = evolve(lo11.state, lo11, 1.0)
    st2 = evolve(lo22.state, lo22, 1.0)
    st3 = evolve(lo33.state, lo33, 1.0)
    @test st1 ≈ st2
    @test st1 ≈ st3
    # Test evolve(state,keywords):
    st1 = evolve(lo11.state, lo11, 1.0; diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    st2 = evolve(lo22.state, lo22, 1.0, diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    st3 = evolve(lo33.state, lo33, 1.0, diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    @test st1 ≈ st2
    @test st1 ≈ st3
    # Test evolve(system):
    s1 = deepcopy(lo11); s2 = deepcopy(lo22); s3 = deepcopy(lo33)
    evolve!(s1, 1.0)
    evolve!(s2, 1.0)
    evolve!(s3, 1.0)
    @test s1.state ≈ s2.state
    @test s1.state ≈ s3.state
  end

  @testset "Timeseries" begin
    # timeseries pure:
    ts1 = timeseries(lo11, 5.0)
    ts2 = timeseries(lo22, 5.0)
    ts3 = timeseries(lo33, 5.0)
    @test eltype(ts3) == BigFloat
    @test size(ts1) == size(ts2) == size(ts3)
    @test ts1[end, :] ≈ ts2[end,:]
    @test ts1[end, :] ≈ ts3[end,:]
    # timeseries with dt:
    ts1 = timeseries(lo11, 5.0, 0.1; mutate = false)
    ts2 = timeseries(lo22, 5.0, 0.1; mutate = false)
    ts3 = timeseries(lo33, 5.0, 0.1; mutate = false)
    @test size(ts1) == size(ts2) == size(ts3)
    @test ts1[end, :] ≈ ts2[end,:]
    @test ts1[end, :] ≈ ts3[end,:]
    # timeseries with diff_eq_kwargs:
    ts1 = timeseries(lo11, 5.0, 0.1, diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9), mutate = false)
    ts2 = timeseries(lo22, 5.0, 0.1, diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9), mutate = false)
    ts3 = timeseries(lo33, 5.0, 0.1, diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9), mutate = false)
    @test size(ts1) == size(ts2) == size(ts3)
    @test ts1[end, :] ≈ ts2[end,:]
    @test ts1[end, :] ≈ ts3[end,:]
    # timeseries with kwargs but without dt:
    ts1 = timeseries(lo11, 5.0, diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9), mutate = false)
    ts2 = timeseries(lo22, 5.0, diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9), mutate = false)
    ts3 = timeseries(lo33, 5.0, diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9), mutate = false)
    @test size(ts1) == size(ts2) == size(ts3)
    @test ts1[end, :] ≈ ts2[end,:]
    @test ts1[end, :] ≈ ts3[end,:]
  end

  @testset "Jacobians" begin
    j1 = jacobian(lo11); j2 = jacobian(lo22); j3 = jacobian(lo33)
    @test eltype(j3) == BigFloat
    @test j1 ≈ j2
    @test j1 ≈ j3
    s1 = evolve(lo11, 1.0)
    s2 = evolve(lo22, 1.0)
    s3 = evolve(lo33, 1.0)
    j1 = lo11.jacob(s1); j2 = lo22.jacob(s2); j3 = lo33.jacob(s3)
    @test eltype(j3) == BigFloat
    @test j1 ≈ j2
    @test j1 ≈ j3
  end

end
