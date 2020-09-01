using Test
using UnderwaterAcoustics
using UnderwaterAcoustics: amp2db, db2amp

@testset "basic" begin

  @test soundspeed() ≈ 1539.0 atol=0.1
  @test soundspeed(; voidfrac=1e-5) ≈ 1402.1 atol=0.1
  @test soundspeed(; voidfrac=1.0) ≈ 340.0

  @test amp2db(absorption(10000, 1000.0, 35.0, 15.0)) ≈ -1.0 atol=0.5
  @test amp2db(absorption(50000)) ≈ -11.0 atol=0.5
  @test amp2db(absorption(100000)) ≈ -36.0 atol=0.5
  @test amp2db(absorption(100000, 1000.0, 38.5, 14.0, 0.0)) ≈ -40.0 atol=0.5
  @test amp2db(absorption(100000, 1000.0, 38.5, 14.0, 2000.0)) ≈ -30.0 atol=0.5
  @test amp2db(absorption(100000, 1000.0, 38.5, 14.0, 6000.0)) ≈ -16.0 atol=0.5

  @test waterdensity() ≈ 1022.7 atol=0.1

  @test reflectioncoef(0.0, 1200.0/1023.0, 1600.0/1540.0, 0.0) isa Complex
  @test reflectioncoef(0.0, 1200.0/1023.0, 1600.0/1540.0, 0.0) ≈ 0.0986 atol=0.0001
  @test amp2db(abs(reflectioncoef(0.0, 2.5, 2.5, 0.01374))) ≈ -2.8 atol=0.1
  @test amp2db(abs(reflectioncoef(0.0, 2.492, 1.3370, 0.01705))) ≈ -5 atol=0.5
  @test amp2db(abs(reflectioncoef(0.0, 1.195, 1.0179, 0.02158))) ≈ -20 atol=0.5
  @test amp2db(abs(reflectioncoef(1.22, 1.149, 0.9873, 0.00386))) ≈ -32 atol=0.5

  @test amp2db(surfaceloss(15.0, 20000.0, 80°)) ≈ -6.5 atol=0.1
  @test amp2db(surfaceloss(10.0, 20000.0, 80°)) ≈ -3.4 atol=0.1
  @test amp2db(surfaceloss(5.0, 20000.0, 80°)) ≈ -0.5 atol=0.1

  @test doppler(0.0, 50000.0) == 50000.0
  @test doppler(10.0, 50000.0) ≈ 50325 atol=0.5
  @test doppler(-10.0, 50000.0) ≈ 49675 atol=0.5

  @test bubbleresonance(100e-6) ≈ 32500.0 atol=0.1
  @test bubbleresonance(32e-6) ≈ 101562.5 atol=0.1
  @test bubbleresonance(100e-6, 10.0) ≈ 45962.0 atol=0.1
end
