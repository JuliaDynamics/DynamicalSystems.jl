using DynamicalSystems, Base.Test

test_value(val, vmin, vmax) = @test vmin <= val <= vmax

println("\nTesting delay coordinates reconstruction...")
@testset "Delay Coords" begin
  he = Systems.henon()
  ts = timeseries(he, 20000)
  x = ts[:, 1] # some "recorded" timeseries
  R = reconstruct(x, 3, 3)
  @test length(size(R)) == 2
  D2 = information_dim(R)
  test_value(D2, 1, 2)
end
