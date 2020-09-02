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

@testset "pm-core-basic" begin

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
  sig = record(noise(env), 1.0, 44100.0)
  @test length(sig) == 44100
  @test sum(abs2.(sig)) > 0.0

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
  src = Pinger(10.0, -10.0, 1000.0; sourcelevel=1.0, interval=0.5)
  sig = record(src, 1.0, 44100.0)
  @test src isa AcousticSource
  @test location(src) == (10.0, 0.0, -10.0)
  @test nominalfrequency(src) == 1000.0
  @test phasor(src) == complex(1.0, 0.0)
  @test length(sig) == 44100
  @test eltype(sig) <: Complex
  @test sum(abs2.(sig))/44100.0 ≈ 0.04 atol=1e-4
  @test maximum(abs2.(sig)) ≈ 1.0 atol=1e-4

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

  rx = AcousticReceiverGrid2D(100.0, 2.0, 100, -100.0, 5.0, 10)
  @test rx isa AcousticReceiverGrid2D
  @test rx isa AbstractMatrix
  @test size(rx) == (100, 10)
  @test rx[1,1] == AcousticReceiver(100.0, -100.0)
  @test rx[end,end] == AcousticReceiver(298.0, -55.0)
  rx = AcousticReceiverGrid3D(100.0, 2.0, 100, 0.0, 1.0, 100, -100.0, 5.0, 10)
  @test rx isa AcousticReceiverGrid3D
  @test rx isa AbstractArray
  @test size(rx) == (100, 100, 10)
  @test rx[1,1,1] == AcousticReceiver(100.0, 0.0, -100.0)
  @test rx[end,end,end] == AcousticReceiver(298.0, 99.0, -55.0)

end

@testset "pekeris" begin

  @test PekerisRayModel in models()

  env = UnderwaterEnvironment()
  pm = PekerisRayModel(env, 7)
  @test pm isa PekerisRayModel

  arr = arrivals(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(arr) == 7
  @test arr[1].time ≈ 0.0650 atol=0.0001
  @test arr[2].time ≈ 0.0657 atol=0.0001
  @test arr[3].time ≈ 0.0670 atol=0.0001
  @test all([arr[j].time > arr[j-1].time for j ∈ 2:7])
  @test abs(arr[1].phasor) ≈ 0.01 atol=0.001
  @test real(arr[2].phasor) < 0.0
  @test imag(arr[2].phasor) ≈ 0.0
  @test all([abs(arr[j].phasor) < abs(arr[j-1].phasor) for j ∈ 2:7])
  @test [(arr[j].surface, arr[j].bottom) for j ∈ 1:7] == [(0,0), (1,0), (0,1), (1,1), (1,1), (2,1), (1,2)]
  @test all([abs(arr[j].arrivalangle) == abs(arr[j].launchangle) for j ∈ 1:7])

  r = eigenrays(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(r) == 7
  @test all([abs(r[j].arrivalangle) == abs(r[j].launchangle) for j ∈ 1:7])
  @test all([r[j].raypath[1] == (0.0, 0.0, -5.0) for j ∈ 1:7])
  @test all([r[j].raypath[end] == (100.0, 0.0, -10.0) for j ∈ 1:7])
  @test all([length(r[j].raypath) == r[j].surface + r[j].bottom + 2 for j ∈ 1:7])

  r = rays(pm, AcousticSource(0.0, -5.0, 1000.0), -60°:15°:60°, 100.0)
  @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(r) == 9
  @test all([r[j].launchangle for j ∈ 1:9] .≈ -60°:15°:60°)
  @test all([abs(r[j].arrivalangle) == abs(r[j].launchangle) for j ∈ 1:9])
  @test all([r[j].raypath[1] == (0.0, 0.0, -5.0) for j ∈ 1:9])
  @test all([r[j].raypath[end][1] ≥ 100.0 for j ∈ 1:9])
  @test all([length(r[j].raypath) == r[j].surface + r[j].bottom + 2 for j ∈ 1:9])

  ir1 = impulseresponse(arr, 10000.0; reltime=true)
  ir2 = impulseresponse(arr, 10000.0; reltime=false)
  @test length(ir2) ≈ length(ir1) + round(Int, 10000.0 * arr[1].time) atol=1
  @test length(ir2) == round(Int, 10000.0 * arr[end].time) + 1
  @test sum(ir1 .!= 0.0) == 7
  @test sum(ir2 .!= 0.0) == 7
  ndx = findall(abs.(ir1) .> 0)
  @test (ndx .- 1) ./ 10000.0 ≈ [arr[j].time - arr[1].time for j ∈ 1:7] atol=1e-4
  ndx = findall(abs.(ir2) .> 0)
  @test (ndx .- 1) ./ 10000.0 ≈ [arr[j].time for j ∈ 1:7] atol=1e-4

  env = UnderwaterEnvironment(ssp=IsoSSP(1500.0))
  pm = PekerisRayModel(env, 2)
  d = (√1209.0)/4.0
  x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
  @test x isa Complex
  @test abs(x) ≈ 0.0 atol=0.0002
  x′ = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
  @test x′ isa Complex
  @test imag(x′) == 0.0
  @test abs(x′) > 1/100.0
  d = (√2409.0)/8.0
  x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
  @test abs(x) > abs(x′)
  y = transmissionloss(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
  @test abs(x) ≈ abs(x′) atol=0.0001
  y = transmissionloss(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x1 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
  x2 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  x3 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]
  x1 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
  x2 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  x3 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]

  env = UnderwaterEnvironment()
  pm = PekerisRayModel(env, 7)
  tx = AcousticSource(0.0, -5.0, 1000.0)
  rx = AcousticReceiver(100.0, -10.0)
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
  rx = AcousticReceiver(100.0, -10.0)
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
  rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  tx = AcousticSource(0.0, -5.0, 1000.0)
  rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)

  env = UnderwaterEnvironment(noise=missing)
  pm = PekerisRayModel(env, 7)
  tx = Pinger(0.0, -5.0, 1000.0; interval=0.3)
  rx = AcousticReceiver(100.0, -10.0)
  sig1 = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig1) == (44100,)
  sig2 = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig2) == (44100,)
  @test sig1[22051:end] ≈ sig2[1:22050]
  rx = AcousticReceiver(100.0, -11.0)
  sig3 = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig3) == (44100,)
  @test !(sig1 ≈ sig3)
  rx = AcousticReceiver(100.0/√2, 100.0/√2, -10.0)
  sig3 = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig3) == (44100,)
  @test sig1 ≈ sig3
  tx = [Pinger(0.0, -5.0, 1000.0; interval=0.3), Pinger(1.0, -5.0, 2000.0; interval=0.5)]
  rx = AcousticReceiver(100.0, 0.0, -10.0)
  sig1 = record(pm, tx, rx, 1.0, 44100.0)
  rx = AcousticReceiver(100.0/√2, 100.0/√2, -10.0)
  rx = AcousticReceiver(-100.0, 0.0, -10.0)
  sig2 = record(pm, tx, rx, 1.0, 44100.0)
  @test !(sig1 ≈ sig2)

end

@testset "bellhop" begin

  if Bellhop in models()

    env = UnderwaterEnvironment(seasurface=Vacuum)
    pm = Bellhop(env)
    @test pm isa Bellhop

    arr = arrivals(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
    @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
    r = eigenrays(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
    @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
    r = rays(pm, AcousticSource(0.0, -5.0, 1000.0), -60°:15°:60°, 100.0)
    @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
    x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
    @test x isa Complex
    y = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
    @test -10 * log10(abs2(x)) ≈ y atol=0.1
    x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0); mode=:incoherent)
    @test x isa Complex
    @test imag(x) == 0.0
    y = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0); mode=:incoherent)
    @test -10 * log10(abs2(x)) ≈ y atol=0.1
    x1 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
    x2 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
    x3 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
    x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
    @test x isa AbstractVector
    @test [x1, x2, x3] == x
    x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (1, 3)
    @test [x1 x2 x3] == x
    x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (3, 3)
    @test [x1, x2, x3] == x[1,:]
    x1 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
    x2 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
    x3 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
    x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
    @test x isa AbstractVector
    @test [x1, x2, x3] == x
    x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (1, 3)
    @test [x1 x2 x3] == x
    x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (3, 3)
    @test [x1, x2, x3] == x[1,:]

    tx = AcousticSource(0.0, -5.0, 1000.0)
    rx = AcousticReceiver(100.0, -10.0)
    sig = record(pm, tx, rx, 1.0, 44100.0)
    @test size(sig) == (44100,)
    sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,)
    sig = recorder(pm, tx, rx)(1.0, 44100.0)
    @test size(sig) == (44100,)
    sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,)
    tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
    rx = AcousticReceiver(100.0, -10.0)
    sig = record(pm, tx, rx, 1.0, 44100.0)
    @test size(sig) == (44100,)
    sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,)
    sig = recorder(pm, tx, rx)(1.0, 44100.0)
    @test size(sig) == (44100,)
    sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,)
    tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
    rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
    sig = record(pm, tx, rx, 1.0, 44100.0)
    @test size(sig) == (44100,2)
    sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,2)
    sig = recorder(pm, tx, rx)(1.0, 44100.0)
    @test size(sig) == (44100,2)
    sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,2)
    tx = AcousticSource(0.0, -5.0, 1000.0)
    rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
    sig = record(pm, tx, rx, 1.0, 44100.0)
    @test size(sig) == (44100,2)
    sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,2)
    sig = recorder(pm, tx, rx)(1.0, 44100.0)
    @test size(sig) == (44100,2)
    sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,2)

  else
    @test_skip true
  end

end

@testset "raysolver" begin

  @test RaySolver in models()

  env = UnderwaterEnvironment()
  pm = RaySolver(env)
  @test pm isa RaySolver

  arr = arrivals(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(arr) >= 7
  @test arr[1].time ≈ 0.0650 atol=0.0001
  @test arr[2].time ≈ 0.0657 atol=0.0001
  @test arr[3].time ≈ 0.0670 atol=0.0001
  @test all([arr[j].time > arr[j-1].time for j ∈ 2:7])
  @test abs(arr[1].phasor) ≈ 0.01 atol=0.001
  @test real(arr[2].phasor) < 0.0
  @test imag(arr[2].phasor) ≈ 0.0
  @test all([abs(arr[j].phasor) < abs(arr[j-1].phasor) for j ∈ 2:7])
  @test [(arr[j].surface, arr[j].bottom) for j ∈ 1:7] == [(0,0), (1,0), (0,1), (1,1), (1,1), (2,1), (1,2)]
  @test abs.([a.arrivalangle for a ∈ arr]) ≈ abs.([a.launchangle for a ∈ arr])

  r = eigenrays(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(r) >= 7
  @test abs.([a.arrivalangle for a ∈ r]) ≈ abs.([a.launchangle for a ∈ r])
  @test all([r[j].raypath[1] == (0.0, 0.0, -5.0) for j ∈ 1:7])
  #@test all([r[j].raypath[end][k] .≈ (100.0, 0.0, -10.0)[k] for j ∈ 1:7, k ∈ 1:3])

  r = rays(pm, AcousticSource(0.0, -5.0, 1000.0), -60°:15°:60°, 100.0)
  @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(r) == 9
  @test all([r[j].launchangle for j ∈ 1:9] .≈ -60°:15°:60°)
  @test abs.([a.arrivalangle for a ∈ r]) ≈ abs.([a.launchangle for a ∈ r])
  @test all([r[j].raypath[1] == (0.0, 0.0, -5.0) for j ∈ 1:9])
  @test r[4].raypath[end][1] ≥ 99.9
  @test r[5].raypath[end][1] ≥ 99.9
  @test r[6].raypath[end][1] ≥ 99.9
  @test r[7].raypath[end][1] ≥ 99.9

  env = UnderwaterEnvironment(ssp=IsoSSP(1500.0), seabed=Rayleigh(1.0, 1.0))
  pm = RaySolver(env)
  d = (√1209.0)/4.0
  x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
  @test x isa Complex
  @test abs(x) ≈ 0.0 atol=0.0002
  x′ = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
  @test x′ isa Complex
  @test imag(x′) == 0.0
  @test abs(x′) > 1/100.0
  d = (√2409.0)/8.0
  x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
  @test abs(x) > abs(x′)
  y = transmissionloss(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
  @test abs(x) ≈ abs(x′) atol=0.0001
  y = transmissionloss(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x1 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
  x2 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  x3 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] ≈ x atol=0.01
  y = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
  @test y isa AbstractMatrix
  @test size(y) == (3, 3)
  @test x' ≈ y[1,:] atol=0.05
  x1 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
  x2 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  x3 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1, x2] ≈ x[1:2] atol=1.0
  y = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
  @test y isa AbstractMatrix
  @test size(y) == (3, 3)
  @test x' ≈ y[1,:] atol=0.1

  tx = AcousticSource(0.0, -5.0, 1000.0)
  rx = AcousticReceiver(100.0, -10.0)
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
  rx = AcousticReceiver(100.0, -10.0)
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
  rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  tx = AcousticSource(0.0, -5.0, 1000.0)
  rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)

end
