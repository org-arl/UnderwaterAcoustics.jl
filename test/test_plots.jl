using TestItems

@testitem "plots-fields" begin
  using RecipesBase, Colors
  # normally defined by Plots
  hasmethod(RecipesBase.is_key_supported, Tuple{Symbol}) || @eval RecipesBase.is_key_supported(::Symbol) = true
  series(args...; kw...) = only(RecipesBase.apply_recipe(Dict{Symbol,Any}(kw...), args...)).args
  bathy = SampledField([100.0, 200.0]; x=[0.0, 1000.0])
  x, y = series(bathy, 0, 1000; npts=3)
  @test x ≈ [0.0, 500.0, 1000.0]
  @test y ≈ [100.0, 150.0, 200.0]
  x, y = series(bathy, 0:1000; npts=3)
  @test x ≈ [0.0, 500.0, 1000.0]
  @test y ≈ [100.0, 150.0, 200.0]
  ssp = SampledField([1500.0 1510.0 1520.0; 1530.0 1540.0 1550.0]; x=[0.0, 1000.0], z=[-100.0, -50.0, 0.0])
  x, z, v = series(ssp; npts=3)
  @test v ≈ [1500.0 1515.0 1530.0; 1510.0 1525.0 1540.0; 1520.0 1535.0 1550.0]
  x, z, v = series(ssp, [0.0, 1000.0], [-100.0, -50.0, 0.0])
  @test (x, z) == ([0.0, 1000.0], [-100.0, -50.0, 0.0])
  @test v ≈ [1500.0 1530.0; 1510.0 1540.0; 1520.0 1550.0]
  @test last(series(bathy, [0.0, 500.0], [-10.0, 0.0])) ≈ [100.0 150.0; 100.0 150.0]
  fld = SampledField([1.0 2.0 3.0; 4.0 5.0 6.0]; x=[0.0, 1000.0], y=[-10.0, 0.0, 10.0])
  x, y, v = series(fld; npts=2)
  @test v ≈ [1.0 4.0; 3.0 6.0]
end
