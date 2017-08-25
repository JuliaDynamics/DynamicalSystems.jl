using DynamicalSystems, Base.Test

test_value(val, vmin, vmax) = @test vmin <= val <= vmax

println("\nTesting delay coordinates reconstruction...")
@testset "Delay Coords" begin
  he = Systems.henon()
  ts = timeseries(he, 20000)
  x = ts[:, 1] # some "recorded" timeseries
  R = reconstruct(x, 5, 3)
  @test size(R) == (19990, 3)
  D2 = information_dim(R) #around 1.17
  test_value(D2, 1, 2)

  x = view(ts, :, 1)
  R = reconstruct(x, 5, 3)
  @test size(R) == (19990, 3)
  D2 = information_dim(R) #around 1.17
end
