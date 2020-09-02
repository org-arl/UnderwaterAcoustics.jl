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

@testset "pm-core" begin

  env = UnderwaterEnvironment()

  @test models() isa AbstractArray
  @test models(env) isa AbstractArray

  @test env isa UnderwaterAcoustics.BasicUnderwaterEnvironment
  @test altimetry(env) isa Altimetry
  @test bathymetry(env) isa Bathymetry
  @test ssp(env) isa SoundSpeedProfile
  @test salinity(env) isa Real
  @test seasurface(env) isa ReflectionModel
  @test seabed(env) isa ReflectionModel
  @test noise(env) isa NoiseModel

  @test altitude(altimetry(env), 0.0, 0.0) isa Real
  @test depth(bathymetry(env), 0.0, 0.0) isa Real
  @test soundspeed(ssp(env), 0.0, 0.0, 0.0) isa Real
  @test reflectioncoef(seasurface(env), 1000.0, 0.0) isa Complex
  @test reflectioncoef(seabed(env), 1000.0, 0.0) isa Complex
  @test length(record(noise(env), 1.0, 44100.0)) == 44100

  @test AcousticReceiver(0.0, 0.0) isa AcousticReceiver
  @test location(AcousticReceiver(100.0, 10.0, -50.0)) == (100, 10.0, -50.0)
  @test location(AcousticReceiver(100.0, -50.0)) == (100, 0.0, -50.0)

  src = AcousticSource(100.0, 10.0, -50.0, 4410.0; sourcelevel=1.0)
  sig = record(src, 1.0, 44100.0)
  @test src isa NarrowbandAcousticSource
  @test location(src) == (100.0, 10.0, -50.0)
  @test location(AcousticSource(100.0, -50.0, 1000.0)) == (100.0, 0.0, -50.0)
  @test nominalfrequency(src) == 4410.0
  @test phasor(src) == complex(1.0, 0.0)
  @test length(sig) == 44100
  @test eltype(sig) <: Complex
  @test sum(abs2.(sig))/44100.0 ≈ 1.0 atol=1e-4
  @test sig[1] != sig[2]
  @test sig[1] ≈ sig[11]
  @test sig[2] ≈ sig[12]

  @test IsoSSP(1500.0) isa SoundSpeedProfile
  @test soundspeed(IsoSSP(1500.0), 100.0, 10.0, -50.0) == 1500.0
  @test MunkSSP() isa SoundSpeedProfile
  @test soundspeed(MunkSSP(), 0.0, 0.0, -2000.0) ≈ 1505.0 atol=1.0
  @test soundspeed(MunkSSP(), 0.0, 0.0, -3000.0) ≈ 1518.0 atol=1.0
  s = SampledSSP(0.0:500.0:1000.0, [1500.0, 1520.0, 1510.0])
  @test s isa SoundSpeedProfile
  @test soundspeed(s, 0.0, 0.0, 0.0) ≈ 1500.0
  @test soundspeed(s, 0.0, 0.0, -250.0) ≈ 1510.0
  @test soundspeed(s, 0.0, 0.0, -500.0) ≈ 1520.0
  @test soundspeed(s, 0.0, 0.0, -750.0) ≈ 1515.0
  @test soundspeed(s, 0.0, 0.0, -1000.0) ≈ 1510.0
  s = SampledSSP(0.0:500.0:1000.0, [1500.0, 1520.0, 1510.0], :cubic)
  @test s isa SoundSpeedProfile
  @test soundspeed(s, 0.0, 0.0, 0.0) ≈ 1500.0
  @test soundspeed(s, 0.0, 0.0, -250.0) > 1510.0
  @test soundspeed(s, 0.0, 0.0, -500.0) ≈ 1520.0
  @test soundspeed(s, 0.0, 0.0, -750.0) > 1515.0
  @test soundspeed(s, 0.0, 0.0, -1000.0) ≈ 1510.0

  @test ConstantDepth(20.0) isa Bathymetry
  @test depth(ConstantDepth(20.0), 0.0, 0.0) == 20.0
  @test maxdepth(ConstantDepth(20.0)) == 20.0
  b = SampledDepth(0.0:500.0:1000.0, [20.0, 25.0, 15.0])
  @test b isa Bathymetry
  @test depth(b, 0.0, 0.0) ≈ 20.0
  @test depth(b, 250.0, 0.0) ≈ 22.5
  @test depth(b, 500.0, 0.0) ≈ 25.0
  @test depth(b, 750.0, 0.0) ≈ 20.0
  @test depth(b, 1000.0, 0.0) ≈ 15.0
  @test maxdepth(b) ≈ 25.0
  b = SampledDepth(0.0:500.0:1000.0, [20.0, 25.0, 15.0], :cubic)
  @test b isa Bathymetry
  @test depth(b, 0.0, 0.0) ≈ 20.0
  @test depth(b, 250.0, 0.0) > 22.5
  @test depth(b, 500.0, 0.0) ≈ 25.0
  @test depth(b, 750.0, 0.0) > 20.0
  @test depth(b, 1000.0, 0.0) ≈ 15.0
  @test maxdepth(b) ≈ 25.0

  @test FlatSurface() isa Altimetry
  @test altitude(FlatSurface(), 0.0, 0.0) == 0.0

  @test ReflectionCoef(0.5 + 0.3im) isa ReflectionModel
  @test reflectioncoef(ReflectionCoef(0.5 + 0.3im), 1000.0, 0.0) == 0.5 + 0.3im
  @test Rayleigh(1.0, 1.0) isa ReflectionModel
  @test Rayleigh(1.0, 1.0, 0.0) isa ReflectionModel
  @test reflectioncoef(Rayleigh(1.0, 1.0, 0.0), 1000.0, 0.0) ≈ 0.0
  @test reflectioncoef(Rayleigh(0.0, 1.0, 0.0), 1000.0, 0.0) ≈ -1.0
  @test Rock isa Rayleigh
  @test Pebbles isa Rayleigh
  @test SandyGravel isa Rayleigh
  @test CoarseSand isa Rayleigh
  @test MediumSand isa Rayleigh
  @test FineSand isa Rayleigh
  @test VeryFineSand isa Rayleigh
  @test ClayeySand isa Rayleigh
  @test CoarseSilt isa Rayleigh
  @test SandySilt isa Rayleigh
  @test Silt isa Rayleigh
  @test FineSilt isa Rayleigh
  @test SandyClay isa Rayleigh
  @test SiltyClay isa Rayleigh
  @test Clay isa Rayleigh
  @test Vacuum isa ReflectionModel
  @test reflectioncoef(Vacuum, 1000.0, 0.0) ≈ -1.0
  @test 0.0 < abs(reflectioncoef(Rock, 1000.0, 0.0)) < 1.0
  @test abs(reflectioncoef(Silt, 1000.0, 0.0)) < abs(reflectioncoef(SandyGravel, 1000.0, 0.0))
  @test abs(reflectioncoef(CoarseSilt, 1000.0, 0.0)) < abs(reflectioncoef(Rock, 1000.0, 0.0))
  @test SurfaceLoss(5.0) isa ReflectionModel
  @test SeaState0 isa SurfaceLoss
  @test SeaState1 isa SurfaceLoss
  @test SeaState2 isa SurfaceLoss
  @test SeaState3 isa SurfaceLoss
  @test SeaState4 isa SurfaceLoss
  @test SeaState5 isa SurfaceLoss
  @test SeaState6 isa SurfaceLoss
  @test SeaState7 isa SurfaceLoss
  @test SeaState8 isa SurfaceLoss
  @test SeaState9 isa SurfaceLoss
  @test 0.0 < abs(reflectioncoef(SeaState2, 1000.0, 0.0)) < 1.0
  @test abs(reflectioncoef(SeaState2, 1000.0, 0.0)) > abs(reflectioncoef(SeaState5, 1000.0, 0.0))
  @test abs(reflectioncoef(SeaState5, 1000.0, 0.0)) > abs(reflectioncoef(SeaState5, 2000.0, 0.0))

  # Pinger
  # RedNoise details
  # AcousticReceiverGrid2D, AcousticReceiverGrid3D
end

@testset "pekeris" begin

  @test PekerisRayModel in models()

  env = UnderwaterEnvironment()
  #pm = PekerisRayModel(env, )

  # impulse response

end
