@testset "Discrete Lyapunovs" begin
  @testset "Towel Map" begin
    @testset "λspectrum" begin
      ds = DynamicalSystems.Systems.towel()
      λ1 = λspectrum(ds, 1e6)
      @test 0.42 < λ1[1] < 0.44
      @test 0.36 < λ1[2] < 0.38
      @test -3.4 < λ1[3] < -3.2
    end

    @testset "λspectrum ForwardDiff" begin
      ds = DynamicalSystems.Systems.towel()
      ds = DiscreteDS(ds.state, ds.eom)
      λ1 = λspectrum(ds, 1e6)
      @test 0.42 < λ1[1] < 0.44
      @test 0.36 < λ1[2] < 0.38
      @test -3.4 < λ1[3] < -3.1
    end

    # test lmax
    # @testset "λmax" begin
    #   #aasd
    # end
    #
    # @testset "λspectrum benchmarks" begin
    #   #asdf
    # end
  end
end
