if current_module() != DynamicalSystems
  using DynamicalSystems
end
using Base.Test, StaticArrays

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

  end

  @testset "trajectory" begin
    # trajectory pure:
    ts1 = trajectory(lo11, 2.0)
    ts3 = trajectory(lo33, 2.0)
    @test eltype(ts3[1]) == BigFloat
    @test size(ts1) == size(ts3)
    @test ts1[end, :] ≈ ts3[end,:]
    # trajectory with diff_eq_kwargs and dt:
    ts1 = trajectory(lo11, 2.0; dt=0.1,
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    ts3 = trajectory(lo33, 2.0; dt=0.1,
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    @test size(ts1) == size(ts3)
    @test ts1[end, :] ≈ ts3[end,:]
  end

  lo11 = Systems.lorenz()
  lo22 = ContinuousDS(lo11.state, lo11.eom!)
  lo33 = Systems.lorenz(big.([0.0, 10.0, 0.0]))

  @testset "Jacobians" begin
    j1 = lo11.J
    j2 = lo22.J
    j3 = lo33.J
    @test eltype(j3) == BigFloat
    @test j1 ≈ j2
    @test j1 ≈ j3
    s1 = evolve(lo11, 1.0)
    s2 = evolve(lo22, 1.0)
    s3 = evolve(lo33, 1.0)
    j1 = (lo11.jacob!(lo11.J, s1); lo11.J)
    j2 = (lo22.jacob!(lo22.J, s2); lo22.J)
    j3 = (lo33.jacob!(lo33.J, s3); lo33.J)
    @test eltype(j3) == BigFloat
    @test j1 ≈ j2
    @test j1 ≈ j3
  end

end

@testset "Lorenz96" begin
    u = ones(5)
    lo11 = Systems.lorenz96(5, u)
    lo33 = Systems.lorenz96(5, big.(ones(5)))
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
    end

    @testset "trajectory" begin
      # trajectory pure:
      ts1 = trajectory(lo11, 2.0)
      ts3 = trajectory(lo33, 2.0)
      @test eltype(ts3[1]) == BigFloat
      @test size(ts1) == size(ts3)
      @test ts1[end, :] ≈ ts3[end,:]
      # trajectory with diff_eq_kwargs and dt:
      ts1 = trajectory(lo11, 2.0; dt=0.1,
      diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
      ts3 = trajectory(lo33, 2.0; dt=0.1,
      diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
      @test size(ts1) == size(ts3)
      @test ts1[end, :] ≈ ts3[end,:]
    end
end

@testset "Old Continuous" begin
  function lorenzo(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8//3)
    @inline @inbounds function eom_lorenz!(du, u)
        du[1] = σ*(u[2]-u[1])
        du[2] = u[1]*(ρ-u[3]) - u[2]
        du[3] = u[1]*u[2] - β*u[3]
    end
    i = one(eltype(u0))
    o = zero(eltype(u0))
    J = zeros(eltype(u0), 3, 3)
    J[1,:] .= [-σ*i,        σ*i,    o]
    J[2,:] .= [ρ*i - u0[3],   -i,   -u0[1]]
    J[3,:] .= [u0[2],        u0[1],   (-β*i)]

    @inline @inbounds function jacob_lorenz!(J, u)
        J[2,1] = ρ - u[3]
        J[2,3] = -u[1]
        J[3,1] = u[2]; J[3,2] = u[1]
    end
    name = "Old Lorenz63 system"
    return ContinuousDS(u0, eom_lorenz!, jacob_lorenz!, J; name = name)
  end

  ds = lorenzo()
  lo33 = lorenzo(big.([0.0, 10.0, 0.0]))

  @testset "Construction" begin
    @test typeof(ds) <: ContinuousDS
    @test typeof(lo33) <: ContinuousDS
  end
  @testset "Evolve & kwargs" begin
    # test evolve(system):
    st1 = evolve(ds, 1.0)
    st3 = evolve(lo33, 1.0)
    @test st1 ≈ st3
    # Test evolve(system,keywords):
    st1 = evolve(ds, 1.0;
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    st3 = evolve(lo33, 1.0;
    diff_eq_kwargs=Dict(:abstol=>1e-9, :reltol=>1e-9))
    @test st1 ≈ st3

  end
end
