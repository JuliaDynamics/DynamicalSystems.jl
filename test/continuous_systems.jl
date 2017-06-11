using DynamicalSystems
using Base.Test

@inline eom_lorenz(u) =
SVector{3}(10.0*(u[2]-u[1]), u[1]*(28.0-u[3]) - u[2], u[1]*u[2] - 8/3*u[3])
function jacob_lorenz(u)
  i = one(eltype(u))
  @SMatrix [-10.0*i           10.0*i    zero(i);
            (28.0*i - u[3])   (-i)   (-u[1]);
            u[2]           u[1]   (-i*8/3) ]
end# should give exponents [0.9056, 0, -14.5723]

@testset "Lorenz System" begin

  lo11 = Systems.lorenz()
  lo22 = ContinuousDS([0.0, 10.0, 0.0], eom_lorenz, jacob_lorenz)
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
    st1 = evolve(lo11.state, lo11, 1.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
    st2 = evolve(lo22.state, lo22, 1.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
    st3 = evolve(lo33.state, lo33, 1.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
    @test st1 ≈ st2
    @test st1 ≈ st3
    # Test evolve(system):
    s1 = evolve(lo11, 1.0)
    s2 = evolve(lo22, 1.0)
    s3 = evolve(lo33, 1.0)
    @test s1.state ≈ s2.state
    @test s1.state ≈ s3.state
    # Test evolve(system, keywords):
    s1 = evolve(lo11, 1.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
    s2 = evolve(lo22, 1.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
    s3 = evolve(lo33, 1.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
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
    ts1 = timeseries(lo11, 5.0, 0.1)
    ts2 = timeseries(lo22, 5.0, 0.1)
    ts3 = timeseries(lo33, 5.0, 0.1)
    @test size(ts1) == size(ts2) == size(ts3)
    @test ts1[end, :] ≈ ts2[end,:]
    @test ts1[end, :] ≈ ts3[end,:]
    # timeseries with diff_eq_kwargs:
    ts1 = timeseries(lo11, 5.0, 0.1, Dict(:abstol=>1e-9, :reltol=>1e-9))
    ts2 = timeseries(lo22, 5.0, 0.1, Dict(:abstol=>1e-9, :reltol=>1e-9))
    ts3 = timeseries(lo33, 5.0, 0.1, Dict(:abstol=>1e-9, :reltol=>1e-9))
    @test size(ts1) == size(ts2) == size(ts3)
    @test ts1[end, :] ≈ ts2[end,:]
    @test ts1[end, :] ≈ ts3[end,:]
    # timeseries with kwargs but without dt:
    ts1 = timeseries(lo11, 5.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
    ts2 = timeseries(lo22, 5.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
    ts3 = timeseries(lo33, 5.0, Dict(:abstol=>1e-9, :reltol=>1e-9))
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
    j1 = jacobian(s1); j2 = jacobian(s2); j3 = jacobian(s3)
    @test eltype(j3) == BigFloat
    @test j1 ≈ j2
    @test j1 ≈ j3
  end

end
