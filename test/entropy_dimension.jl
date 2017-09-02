using Base.Test, DynamicalSystems

test_value(val, vmin, vmax) = @test vmin <= val <= vmax

println("\nTesting generalized entropy (renyi) & linear scaling...")
@testset "Generalized Dimensions" begin
  @testset "Henon Map" begin
    ds = Systems.henon()
    ts = timeseries(ds, 200000)
    mat = convert(Matrix, ts)
    # Test call with dataset
    renyi(1, 0.001, ts)
    es = logspace(-0, -3, 7)
    dd = zeros(es)
    for q in [0,2,1, 2.56]
      for (i, ee) in enumerate(es)
        dd[i] = renyi(q, ee, mat)
      end
      linr, dim = linear_region(-log.(es), dd)
      test_value(dim, 1.1, 1.3)
    end
  end
  @testset "Lorenz System" begin
    ds = Systems.lorenz()
    ts = timeseries(ds, 5000)
    es = logspace(1, -3, 11)
    dd = zeros(es)
    for q in [0,1,2]
      for (i, ee) in enumerate(es)
        dd[i] = renyi(q, ee, ts)
      end
      linr, dim = linear_region(-log.(es), dd)
      if q == 0
        test_value(dim, 1.7, 1.9)
      else
        test_value(dim, 1.85, 2.0)
      end
    end
  end
end

println("\nTesting dimension calls (all names)...")
@testset "Dimension calls" begin
  ds = Systems.henon()
  ts = timeseries(ds, 20000)
  # Test call with dataset
  test_value(generalized_dim(1.32, ts), 1.1, 1.3)
  test_value(capacity_dim(ts), 1.1, 1.3)
  test_value(information_dim(ts), 1.1, 1.3)
end
