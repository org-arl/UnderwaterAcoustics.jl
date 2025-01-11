using Test
using UnderwaterAcoustics

# using Statistics
# import ForwardDiff, ReverseDiff, Zygote

# units used in the tests
const ° = u"°"

# defined here to avoid DSP.jl dependency
amp2db(x) = 20 * log10(x)
pow2db(x) = 10 * log10(x)
db2amp(x) = 10^(x/20)

@testset "basic" begin

  @test soundspeed() ≈ 1539.0 atol=0.1
  @test soundspeed(; γ=1e-5) ≈ 1402.1 atol=0.1
  @test soundspeed(; γ=1.0) ≈ 340.0

  @test amp2db(absorption(10000, 1000.0, 35.0, 15.0)) ≈ -1.0 atol=0.5
  @test amp2db(absorption(50000)) ≈ -11.0 atol=0.5
  @test amp2db(absorption(100000)) ≈ -36.0 atol=0.5
  @test amp2db(absorption(100000, 1000.0, 38.5, 14.0, 0.0)) ≈ -40.0 atol=0.5
  @test amp2db(absorption(100000, 1000.0, 38.5, 14.0, 2000.0)) ≈ -30.0 atol=0.5
  @test amp2db(absorption(100000, 1000.0, 38.5, 14.0, 6000.0)) ≈ -16.0 atol=0.5

  @test water_density() ≈ 1022.7 atol=0.1

  @test reflection_coef(0.0, 1200.0/1023.0, 1600.0/1540.0, 0.0) isa Complex
  @test reflection_coef(0.0, 1200.0/1023.0, 1600.0/1540.0, 0.0) ≈ 0.0986 atol=0.0001
  @test amp2db(abs(reflection_coef(0.0, 2.5, 2.5, 0.01374))) ≈ -2.8 atol=0.1
  @test amp2db(abs(reflection_coef(0.0, 2.492, 1.3370, 0.01705))) ≈ -5 atol=0.5
  @test amp2db(abs(reflection_coef(0.0, 1.195, 1.0179, 0.02158))) ≈ -20 atol=0.5
  @test amp2db(abs(reflection_coef(1.22, 1.149, 0.9873, 0.00386))) ≈ -32 atol=0.5

  @test amp2db(surface_reflection_coef(15.0, 20000.0, 80°)) ≈ -6.5 atol=0.1
  @test amp2db(surface_reflection_coef(10.0, 20000.0, 80°)) ≈ -3.4 atol=0.1
  @test amp2db(surface_reflection_coef(5.0, 20000.0, 80°)) ≈ -0.5 atol=0.1

  @test doppler(0.0, 50000.0) == 50000.0
  @test doppler(10.0, 50000.0) ≈ 50325 atol=0.5
  @test doppler(-10.0, 50000.0) ≈ 49675 atol=0.5

  @test bubble_resonance(100e-6) ≈ 32465.562964 atol=0.1
  @test bubble_resonance(32e-6) ≈ 101454.8842 atol=0.1
  @test bubble_resonance(100e-6, 10.0) ≈ 45796.45437634176 atol=0.1

end

# @testset "pm-core-basic" begin

#   env = UnderwaterEnvironment()

#   @test models() isa AbstractArray
#   @test models(env) isa AbstractArray

#   @test env isa UnderwaterAcoustics.BasicUnderwaterEnvironment
#   @test altimetry(env) isa Altimetry
#   @test bathymetry(env) isa Bathymetry
#   @test ssp(env) isa SoundSpeedProfile
#   @test salinity(env) isa Real
#   @test seasurface(env) isa ReflectionModel
#   @test seabed(env) isa ReflectionModel
#   @test noise(env) isa NoiseModel

#   @test altitude(altimetry(env), 0.0, 0.0) isa Real
#   @test depth(bathymetry(env), 0.0, 0.0) isa Real
#   @test soundspeed(ssp(env), 0.0, 0.0, 0.0) isa Real
#   @test reflectioncoef(seasurface(env), 1000.0, 0.0) isa Complex
#   @test reflectioncoef(seabed(env), 1000.0, 0.0) isa Complex
#   sig = record(noise(env), 1.0, 44100.0)
#   @test length(sig) == 44100
#   @test sum(abs2.(sig)) > 0.0

#   @test AcousticReceiver(0.0, 0.0) isa AcousticReceiver
#   @test location(AcousticReceiver(100.0, 10.0, -50.0)) == (100, 10.0, -50.0)
#   @test location(AcousticReceiver(100.0, -50.0)) == (100, 0.0, -50.0)

#   src = AcousticSource(100.0, 10.0, -50.0, 4410.0; sourcelevel=1.0)
#   sig = record(src, 1.0, 44100.0)
#   @test src isa NarrowbandAcousticSource
#   @test location(src) == (100.0, 10.0, -50.0)
#   @test location(AcousticSource(100.0, -50.0, 1000.0)) == (100.0, 0.0, -50.0)
#   @test nominalfrequency(src) == 4410.0
#   @test phasor(src) == complex(1.0, 0.0)
#   @test length(sig) == 44100
#   @test eltype(sig) <: Complex
#   @test sum(abs2.(sig))/44100.0 ≈ 1.0 atol=1e-4
#   @test sig[1] != sig[2]
#   @test sig[1] ≈ sig[11]
#   @test sig[2] ≈ sig[12]
#   src = Pinger(10.0, -10.0, 1000.0; sourcelevel=1.0, interval=0.5)
#   sig = record(src, 1.0, 44100.0)
#   @test src isa AcousticSource
#   @test location(src) == (10.0, 0.0, -10.0)
#   @test nominalfrequency(src) == 1000.0
#   @test phasor(src) == complex(1.0, 0.0)
#   @test length(sig) == 44100
#   @test eltype(sig) <: Complex
#   @test sum(abs2.(sig))/44100.0 ≈ 0.04 atol=1e-4
#   @test maximum(abs2.(sig)) ≈ 1.0 atol=1e-4
#   src = SampledAcousticSource(10.0, -10.0, ones(1000); fs=1000.0, frequency=100.0)
#   @test src isa SampledAcousticSource
#   @test location(src) == (10.0, 0.0, -10.0)
#   @test nominalfrequency(src) == 100.0
#   sig = record(src, 1.0, 1000.0)
#   @test length(sig) == 1000
#   @test eltype(sig) === Float64
#   @test all(sig .== 1.0)
#   sig = record(src, 1.0, 1000.0; start=-0.5)
#   @test length(sig) == 1000
#   @test eltype(sig) === Float64
#   @test all(sig[1:500] .== 0.0)
#   @test all(sig[501:1000] .== 1.0)
#   sig = record(src, 1.0, 1000.0; start=0.5)
#   @test length(sig) == 1000
#   @test eltype(sig) === Float64
#   @test all(sig[1:500] .== 1.0)
#   @test all(sig[501:1000] .== 0.0)
#   sig = record(src, 1.0, 1000.0; start=-2.0)
#   @test length(sig) == 1000
#   @test eltype(sig) === Float64
#   @test all(sig .== 0.0)
#   sig = record(src, 1.0, 1000.0; start=2.0)
#   @test length(sig) == 1000
#   @test eltype(sig) === Float64
#   @test all(sig .== 0.0)
#   @test_throws ArgumentError record(src, 1.0, 2000.0)
#   src = SampledAcousticSource(10.0, 5.0, -10.0, cis.(2π * 1000 * (0:999) ./ 10000.0); fs=10000.0)
#   @test src isa SampledAcousticSource
#   @test location(src) == (10.0, 5.0, -10.0)
#   @test nominalfrequency(src) == 1000.0
#   sig = record(src, 0.1, 10000.0)
#   @test length(sig) == 1000
#   @test eltype(sig) === ComplexF64

#   @test IsoSSP(1500.0) isa SoundSpeedProfile
#   @test soundspeed(IsoSSP(1500.0), 100.0, 10.0, -50.0) == 1500.0
#   @test MunkSSP() isa SoundSpeedProfile
#   @test soundspeed(MunkSSP(), 0.0, 0.0, -2000.0) ≈ 1505.0 atol=1.0
#   @test soundspeed(MunkSSP(), 0.0, 0.0, -3000.0) ≈ 1518.0 atol=1.0
#   s = SampledSSP(0.0:500.0:1000.0, [1500.0, 1520.0, 1510.0])
#   @test s isa SoundSpeedProfile
#   @test soundspeed(s, 0.0, 0.0, 0.0) ≈ 1500.0
#   @test soundspeed(s, 0.0, 0.0, -250.0) ≈ 1510.0
#   @test soundspeed(s, 0.0, 0.0, -500.0) ≈ 1520.0
#   @test soundspeed(s, 0.0, 0.0, -750.0) ≈ 1515.0
#   @test soundspeed(s, 0.0, 0.0, -1000.0) ≈ 1510.0
#   s = SampledSSP(0.0:500.0:1000.0, [1500.0, 1520.0, 1510.0], :smooth)
#   @test s isa SoundSpeedProfile
#   @test soundspeed(s, 0.0, 0.0, 0.0) ≈ 1500.0
#   @test soundspeed(s, 0.0, 0.0, -250.0) > 1510.0
#   @test soundspeed(s, 0.0, 0.0, -500.0) ≈ 1520.0
#   @test soundspeed(s, 0.0, 0.0, -750.0) > 1515.0
#   @test soundspeed(s, 0.0, 0.0, -1000.0) ≈ 1510.0
#   s = SampledSSP(0.0:50.0:100.0, 0.0:500.0:1000.0, [1400.0 1420.0 1410.0; 1500.0 1520.0 1510.0; 1300.0 1320.0 1310.0])
#   @test s isa SoundSpeedProfile
#   @test soundspeed(s, 50.0, 0.0, 0.0) ≈ 1500.0
#   @test soundspeed(s, 50.0, 0.0, -250.0) ≈ 1510.0
#   @test soundspeed(s, 50.0, 0.0, -500.0) ≈ 1520.0
#   @test soundspeed(s, 50.0, 0.0, -750.0) ≈ 1515.0
#   @test soundspeed(s, 50.0, 0.0, -1000.0) ≈ 1510.0
#   @test soundspeed(s, 0.0, 0.0, 0.0) ≈ 1400.0
#   @test soundspeed(s, 0.0, 0.0, -250.0) ≈ 1410.0
#   @test soundspeed(s, 0.0, 0.0, -500.0) ≈ 1420.0
#   @test soundspeed(s, 0.0, 0.0, -750.0) ≈ 1415.0
#   @test soundspeed(s, 0.0, 0.0, -1000.0) ≈ 1410.0

#   @test ConstantDepth(20.0) isa Bathymetry
#   @test depth(ConstantDepth(20.0), 0.0, 0.0) == 20.0
#   @test maxdepth(ConstantDepth(20.0)) == 20.0
#   b = SampledDepth(0.0:500.0:1000.0, [20.0, 25.0, 15.0])
#   @test b isa Bathymetry
#   @test depth(b, 0.0, 0.0) ≈ 20.0
#   @test depth(b, 250.0, 0.0) ≈ 22.5
#   @test depth(b, 500.0, 0.0) ≈ 25.0
#   @test depth(b, 750.0, 0.0) ≈ 20.0
#   @test depth(b, 1000.0, 0.0) ≈ 15.0
#   @test maxdepth(b) ≈ 25.0
#   b = SampledDepth(0.0:500.0:1000.0, [20.0, 25.0, 15.0], :smooth)
#   @test b isa Bathymetry
#   @test depth(b, 0.0, 0.0) ≈ 20.0
#   @test depth(b, 250.0, 0.0) > 22.5
#   @test depth(b, 500.0, 0.0) ≈ 25.0
#   @test depth(b, 750.0, 0.0) > 20.0
#   @test depth(b, 1000.0, 0.0) ≈ 15.0
#   @test maxdepth(b) ≈ 25.0

#   @test FlatSurface() isa Altimetry
#   @test altitude(FlatSurface(), 0.0, 0.0) == 0.0
#   a = SampledAltitude(0.0:500.0:1000.0, [0.0, 1.0, -1.0])
#   @test a isa Altimetry
#   @test altitude(a, 0.0, 0.0) ≈ 0.0
#   @test altitude(a, 250.0, 0.0) ≈ 0.5
#   @test altitude(a, 500.0, 0.0) ≈ 1.0
#   @test altitude(a, 750.0, 0.0) ≈ 0.0
#   @test altitude(a, 1000.0, 0.0) ≈ -1.0
#   a = SampledAltitude(0.0:500.0:1000.0, [0.0, 1.0, -1.0], :smooth)
#   @test a isa Altimetry
#   @test altitude(a, 0.0, 0.0) ≈ 0.0 atol=1e-6
#   @test altitude(a, 250.0, 0.0) > 0.5
#   @test altitude(a, 500.0, 0.0) ≈ 1.0 atol=1e-6
#   @test altitude(a, 1000.0, 0.0) ≈ -1.0 atol=1e-6

#   @test ReflectionCoef(0.5 + 0.3im) isa ReflectionModel
#   @test reflectioncoef(ReflectionCoef(0.5 + 0.3im), 1000.0, 0.0) == 0.5 + 0.3im
#   @test RayleighReflectionCoef(1.0, 1.0) isa ReflectionModel
#   @test RayleighReflectionCoef(1.0, 1.0, 0.0) isa ReflectionModel
#   @test reflectioncoef(RayleighReflectionCoef(1.0, 1.0, 0.0), 1000.0, 0.0) ≈ 0.0
#   @test reflectioncoef(RayleighReflectionCoef(0.0, 1.0, 0.0), 1000.0, 0.0) ≈ -1.0
#   @test Rock isa RayleighReflectionCoef
#   @test Pebbles isa RayleighReflectionCoef
#   @test SandyGravel isa RayleighReflectionCoef
#   @test CoarseSand isa RayleighReflectionCoef
#   @test MediumSand isa RayleighReflectionCoef
#   @test FineSand isa RayleighReflectionCoef
#   @test VeryFineSand isa RayleighReflectionCoef
#   @test ClayeySand isa RayleighReflectionCoef
#   @test CoarseSilt isa RayleighReflectionCoef
#   @test SandySilt isa RayleighReflectionCoef
#   @test Silt isa RayleighReflectionCoef
#   @test FineSilt isa RayleighReflectionCoef
#   @test SandyClay isa RayleighReflectionCoef
#   @test SiltyClay isa RayleighReflectionCoef
#   @test Clay isa RayleighReflectionCoef
#   @test Vacuum isa ReflectionModel
#   @test reflectioncoef(Vacuum, 1000.0, 0.0) ≈ -1.0
#   @test 0.0 < abs(reflectioncoef(Rock, 1000.0, 0.0)) < 1.0
#   @test abs(reflectioncoef(Silt, 1000.0, 0.0)) < abs(reflectioncoef(SandyGravel, 1000.0, 0.0))
#   @test abs(reflectioncoef(CoarseSilt, 1000.0, 0.0)) < abs(reflectioncoef(Rock, 1000.0, 0.0))
#   @test SurfaceLoss(5.0) isa ReflectionModel
#   @test SeaState0 isa SurfaceLoss
#   @test SeaState1 isa SurfaceLoss
#   @test SeaState2 isa SurfaceLoss
#   @test SeaState3 isa SurfaceLoss
#   @test SeaState4 isa SurfaceLoss
#   @test SeaState5 isa SurfaceLoss
#   @test SeaState6 isa SurfaceLoss
#   @test SeaState7 isa SurfaceLoss
#   @test SeaState8 isa SurfaceLoss
#   @test SeaState9 isa SurfaceLoss
#   @test 0.0 < abs(reflectioncoef(SeaState2, 1000.0, 0.0)) < 1.0
#   @test abs(reflectioncoef(SeaState2, 1000.0, 0.0)) > abs(reflectioncoef(SeaState5, 1000.0, 0.0))
#   @test abs(reflectioncoef(SeaState5, 1000.0, 0.0)) > abs(reflectioncoef(SeaState5, 2000.0, 0.0))

#   rx = AcousticReceiverGrid2D(100.0, 2.0, 100, -100.0, 5.0, 10)
#   @test rx isa AcousticReceiverGrid2D
#   @test rx isa AbstractMatrix
#   @test size(rx) == (100, 10)
#   @test rx[1,1] == AcousticReceiver(100.0, -100.0)
#   @test rx[end,end] == AcousticReceiver(298.0, -55.0)
#   rx = AcousticReceiverGrid3D(100.0, 2.0, 100, 0.0, 1.0, 100, -100.0, 5.0, 10)
#   @test rx isa AcousticReceiverGrid3D
#   @test rx isa AbstractArray
#   @test size(rx) == (100, 100, 10)
#   @test rx[1,1,1] == AcousticReceiver(100.0, 0.0, -100.0)
#   @test rx[end,end,end] == AcousticReceiver(298.0, 99.0, -55.0)

# end

# @testset "pm-pekeris" begin

#   @test PekerisRayModel in models()

#   env = UnderwaterEnvironment()
#   pm = PekerisRayModel(env, 7)
#   @test pm isa PekerisRayModel

#   arr = arrivals(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
#   @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
#   @test length(arr) == 7
#   @test arr[1].time ≈ 0.0650 atol=0.0001
#   @test arr[2].time ≈ 0.0657 atol=0.0001
#   @test arr[3].time ≈ 0.0670 atol=0.0001
#   @test all([arr[j].time > arr[j-1].time for j ∈ 2:7])
#   @test abs(arr[1].phasor) ≈ 0.01 atol=0.001
#   @test real(arr[2].phasor) < 0.0
#   @test imag(arr[2].phasor) ≈ 0.0
#   @test all([abs(arr[j].phasor) < abs(arr[j-1].phasor) for j ∈ 2:7])
#   @test [(arr[j].surface, arr[j].bottom) for j ∈ 1:7] == [(0,0), (1,0), (0,1), (1,1), (1,1), (2,1), (1,2)]
#   @test all([abs(arr[j].arrivalangle) == abs(arr[j].launchangle) for j ∈ 1:7])

#   r = eigenrays(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
#   @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
#   @test length(r) == 7
#   @test all([abs(r[j].arrivalangle) == abs(r[j].launchangle) for j ∈ 1:7])
#   @test all([r[j].raypath[1] == (0.0, 0.0, -5.0) for j ∈ 1:7])
#   @test all([r[j].raypath[end] == (100.0, 0.0, -10.0) for j ∈ 1:7])
#   @test all([length(r[j].raypath) == r[j].surface + r[j].bottom + 2 for j ∈ 1:7])

#   r = rays(pm, AcousticSource(0.0, -5.0, 1000.0), -60°:15°:60°, 100.0)
#   @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
#   @test length(r) == 9
#   @test all([r[j].launchangle for j ∈ 1:9] .≈ -60°:15°:60°)
#   @test all([abs(r[j].arrivalangle) == abs(r[j].launchangle) for j ∈ 1:9])
#   @test all([r[j].raypath[1] == (0.0, 0.0, -5.0) for j ∈ 1:9])
#   @test all([r[j].raypath[end][1] ≥ 100.0 for j ∈ 1:9])
#   @test all([length(r[j].raypath) == r[j].surface + r[j].bottom + 2 for j ∈ 1:9])

#   ir1 = impulseresponse(arr, 10000.0; reltime=true, approx=true)
#   ir2 = impulseresponse(arr, 10000.0; reltime=false, approx=true)
#   @test length(ir2) ≈ length(ir1) + round(Int, 10000.0 * arr[1].time) atol=1
#   @test length(ir2) == round(Int, 10000.0 * arr[end].time) + 1
#   @test sum(ir1 .!= 0.0) == 7
#   @test sum(ir2 .!= 0.0) == 7
#   ndx = findall(abs.(ir1) .> 0)
#   @test (ndx .- 1) ./ 10000.0 ≈ [arr[j].time - arr[1].time for j ∈ 1:7] atol=1e-4
#   ndx = findall(abs.(ir2) .> 0)
#   @test (ndx .- 1) ./ 10000.0 ≈ [arr[j].time for j ∈ 1:7] atol=1e-4

#   ir1a = impulseresponse(arr, 10000.0; reltime=true)
#   ir2a = impulseresponse(arr, 10000.0; reltime=false)
#   @test length(ir2a) ≈ length(ir1a) + round(Int, 10000.0 * arr[1].time) atol=1
#   @test length(ir2a) ≥ length(ir2)
#   @test sum(abs2.(ir1a))/sum(abs2.(ir1)) ≈ 1.0 atol=0.05
#   @test sum(abs2.(ir2a))/sum(abs2.(ir2)) ≈ 1.0 atol=0.05

#   @test length(impulseresponse(arr, 10000.0, 256; reltime=true, approx=true)) == 256
#   @test length(impulseresponse(arr, 10000.0, 64; reltime=true, approx=true)) == 64
#   @test length(impulseresponse(arr, 10000.0, 256; reltime=true, approx=false)) == 256
#   @test length(impulseresponse(arr, 10000.0, 64; reltime=true, approx=false)) == 64
#   @test length(impulseresponse(arr, 10000.0, 1024; reltime=false, approx=true)) == 1024
#   @test length(impulseresponse(arr, 10000.0, 700; reltime=false, approx=true)) == 700
#   @test length(impulseresponse(arr, 10000.0, 1024; reltime=false, approx=false)) == 1024
#   @test length(impulseresponse(arr, 10000.0, 700; reltime=false, approx=false)) == 700

#   env = UnderwaterEnvironment(ssp=IsoSSP(1500.0))
#   pm = PekerisRayModel(env, 2)
#   d = (√1209.0)/4.0
#   x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
#   @test x isa Complex
#   @test abs(x) ≈ 0.0 atol=0.0002
#   x′ = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
#   @test x′ isa Complex
#   @test imag(x′) == 0.0
#   @test abs(x′) > 1/100.0
#   d = (√2409.0)/8.0
#   x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
#   @test abs(x) > abs(x′)
#   y = transmissionloss(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
#   @test -10 * log10(abs2(x)) ≈ y atol=0.1
#   x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
#   @test abs(x) ≈ abs(x′) atol=0.0001
#   y = transmissionloss(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
#   @test -10 * log10(abs2(x)) ≈ y atol=0.1
#   x1 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
#   x2 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
#   x3 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
#   x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
#   @test x isa AbstractVector
#   @test [x1, x2, x3] == x
#   x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
#   @test x isa AbstractMatrix
#   @test size(x) == (1, 3)
#   @test [x1 x2 x3] == x
#   x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
#   @test x isa AbstractMatrix
#   @test size(x) == (3, 3)
#   @test [x1, x2, x3] == x[1,:]
#   x1 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
#   x2 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
#   x3 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
#   x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
#   @test x isa AbstractVector
#   @test [x1, x2, x3] == x
#   x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
#   @test x isa AbstractMatrix
#   @test size(x) == (1, 3)
#   @test [x1 x2 x3] == x
#   x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
#   @test x isa AbstractMatrix
#   @test size(x) == (3, 3)
#   @test [x1, x2, x3] == x[1,:]

#   env = UnderwaterEnvironment()
#   pm = PekerisRayModel(env, 7)
#   tx = AcousticSource(0.0, -5.0, 1000.0)
#   rx = AcousticReceiver(100.0, -10.0)
#   sig = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig) == (44100,)
#   sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
#   @test size(sig) == (44100,)
#   sig = recorder(pm, tx, rx)(1.0, 44100.0)
#   @test size(sig) == (44100,)
#   sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
#   @test size(sig) == (44100,)
#   tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
#   rx = AcousticReceiver(100.0, -10.0)
#   sig = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig) == (44100,)
#   sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
#   @test size(sig) == (44100,)
#   sig = recorder(pm, tx, rx)(1.0, 44100.0)
#   @test size(sig) == (44100,)
#   sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
#   @test size(sig) == (44100,)
#   tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
#   rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
#   sig = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig) == (44100,2)
#   sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
#   @test size(sig) == (44100,2)
#   sig = recorder(pm, tx, rx)(1.0, 44100.0)
#   @test size(sig) == (44100,2)
#   sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
#   @test size(sig) == (44100,2)
#   tx = AcousticSource(0.0, -5.0, 1000.0)
#   rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
#   sig = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig) == (44100,2)
#   sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
#   @test size(sig) == (44100,2)
#   sig = recorder(pm, tx, rx)(1.0, 44100.0)
#   @test size(sig) == (44100,2)
#   sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
#   @test size(sig) == (44100,2)

#   env = UnderwaterEnvironment(noise=missing)
#   pm = PekerisRayModel(env, 7)
#   tx = Pinger(0.0, -5.0, 1000.0; interval=0.3)
#   rx = AcousticReceiver(100.0, -10.0)
#   sig1 = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig1) == (44100,)
#   sig2 = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
#   @test size(sig2) == (44100,)
#   @test sig1[22051:end] ≈ sig2[1:22050]
#   rx = AcousticReceiver(100.0, -11.0)
#   sig3 = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig3) == (44100,)
#   @test !(sig1 ≈ sig3)
#   rx = AcousticReceiver(100.0/√2, 100.0/√2, -10.0)
#   sig3 = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig3) == (44100,)
#   @test sig1 ≈ sig3
#   tx = [Pinger(0.0, -5.0, 1000.0; interval=0.3), Pinger(1.0, -5.0, 2000.0; interval=0.5)]
#   rx = AcousticReceiver(100.0, 0.0, -10.0)
#   sig1 = record(pm, tx, rx, 1.0, 44100.0)
#   rx = AcousticReceiver(100.0/√2, 100.0/√2, -10.0)
#   sig2 = record(pm, tx, rx, 1.0, 44100.0)
#   @test !(sig1 ≈ sig2)

#   env = UnderwaterEnvironment(noise=WhiteGaussianNoise(db2amp(120.0)))
#   pm = PekerisRayModel(env, 7)
#   tx = SampledAcousticSource(0.0, -5.0, sin.(10000π .* (1:44100) ./ 44100.0); fs=44100.0)
#   @test tx.frequency ≈ 5000.0 atol=10.0
#   rx = AcousticReceiver(100.0, -10.0)
#   sig = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig) == (44100,)
#   @test eltype(sig) === ComplexF64
#   tx = SampledAcousticSource(0.0, -5.0, cis.(10000π .* (1:44100) ./ 44100.0); fs=44100.0)
#   @test tx.frequency ≈ 5000.0 atol=10.0
#   sig = record(pm, tx, rx, 1.0, 44100.0)
#   @test size(sig) == (44100,)
#   @test eltype(sig) === ComplexF64
#   tx = SampledAcousticSource(0.0, -5.0, zeros(44100); fs=44100.0, frequency=5000.0)
#   sig = record(pm, tx, [rx, rx], 1.0, 44100.0)
#   @test pow2db(mean(abs2, sig[:,1])) ≈ 120.0 atol=0.5
#   @test pow2db(mean(abs2, sig[:,2])) ≈ 120.0 atol=0.5

#   env = UnderwaterEnvironment()
#   pm = PekerisRayModel(env)
#   tx = AcousticSource(0.0, -5.0, 10000.0)
#   rx = AcousticReceiver(100.0, -10.0)
#   rec = recorder(pm, tx, rx)
#   sig = rec(zeros(44100); fs=44100.0)
#   @test sig isa AbstractVector
#   @test eltype(sig) === Float64
#   @test 44100 < length(sig) < 44700
#   sig = rec(zeros(44100); fs=44100.0, reltime=false)
#   @test sig isa AbstractVector
#   @test eltype(sig) === Float64
#   @test 47000 < length(sig) < 47600
#   sig = rec(zeros(ComplexF64, 44100); fs=44100.0)
#   @test sig isa AbstractVector
#   @test eltype(sig) === ComplexF64
#   sig = rec(zeros(ComplexF32, 44100); fs=44100.0)
#   @test sig isa AbstractVector
#   @test eltype(sig) === ComplexF32
#   sig = rec(zeros(Float32, 44100); fs=44100.0)
#   @test sig isa AbstractVector
#   @test eltype(sig) === Float32
#   rec = recorder(pm, [tx, tx], rx)
#   sig = rec(zeros(44100, 2); fs=44100.0)
#   @test sig isa AbstractVector
#   rec = recorder(pm, tx, [rx, rx])
#   sig = rec(zeros(44100); fs=44100.0)
#   @test sig isa AbstractMatrix
#   @test size(sig, 2) == 2
#   rec = recorder(pm, [tx, tx], [rx, rx])
#   sig = rec(zeros(44100, 2); fs=44100.0)
#   @test sig isa AbstractMatrix
#   @test size(sig, 2) == 2

#   env = UnderwaterEnvironment(noise=nothing)
#   pm = PekerisRayModel(env)
#   tx = AcousticSource(0.0, -5.0, 10000.0)
#   rx = AcousticReceiver(100.0, -10.0)
#   rec = recorder(pm, tx, rx)
#   ir = impulseresponse(arrivals(pm, tx, rx), 44100.0, 500; approx=true, reltime=true)
#   X1 = channelmatrix(rec, 44100.0, 500; approx=true)
#   X2 = channelmatrix(arrivals(pm, tx, rx), 44100.0, 500; approx=true)
#   @test X1 == X2
#   @test X1 * vcat([1.0], zeros(size(X1,1)-1)) ≈ ir

# end

# function ∂(f, x, i, ϵ)
#   x1 = copy(x)
#   x1[i] = x[i] + ϵ
#   f1 = f(x1)
#   x1[i] = x[i] - ϵ
#   (f1 - f(x1)) / (2ϵ)
# end

# @testset "pm-∂pekeris" begin

#   function ℳ(x)
#     D, R, d1, d2, f, c = x
#     env = UnderwaterEnvironment(ssp=IsoSSP(c), bathymetry=ConstantDepth(D))
#     pm = PekerisRayModel(env, 7)
#     transmissionloss(pm, AcousticSource(0.0, -d1, f), AcousticReceiver(R, -d2))
#   end

#   x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
#   ∇ℳ = ForwardDiff.gradient(ℳ, x)
#   for i ∈ 1:length(x)
#     @test ∇ℳ[i] ≈ ∂(ℳ, x, i, 0.0001) atol=0.1
#   end

#   x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
#   ∇ℳ = ForwardDiff.gradient(ℳ, x)
#   for i ∈ 1:length(x)
#     @test ∇ℳ[i] ≈ ∂(ℳ, x, i, 0.0001) atol=0.1
#   end

#   x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
#   ∇ℳ = ReverseDiff.gradient(ℳ, x)
#   for i ∈ 1:length(x)
#     @test ∇ℳ[i] ≈ ∂(ℳ, x, i, 0.0001) atol=0.1
#   end

#   x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
#   ∇ℳ = ReverseDiff.gradient(ℳ, x)
#   for i ∈ 1:length(x)
#     @test ∇ℳ[i] ≈ ∂(ℳ, x, i, 0.0001) atol=0.1
#   end

#   x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
#   ∇ℳ = Zygote.gradient(ℳ, x) |> first
#   for i ∈ 1:length(x)
#     @test ∇ℳ[i] ≈ ∂(ℳ, x, i, 0.0001) atol=0.1
#   end

#   x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
#   ∇ℳ = Zygote.gradient(ℳ, x) |> first
#   for i ∈ 1:length(x)
#     @test ∇ℳ[i] ≈ ∂(ℳ, x, i, 0.0001) atol=0.1
#   end

#   function ℳ2(x)
#     D, R, d1, d2, f, c = x
#     env = UnderwaterEnvironment(ssp=IsoSSP(c), bathymetry=ConstantDepth(D))
#     pm = PekerisRayModel(env, 7)
#     # FIXME: approx=true needed because FFT is not yet differentiable
#     X = channelmatrix(arrivals(pm, AcousticSource(0.0, -d1, f), AcousticReceiver(R, -d2)), 44100.0, 500; approx=true)
#     sum(abs2, X)
#   end

#   x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
#   ∇ℳ = ForwardDiff.gradient(ℳ2, x)
#   for i ∈ 1:length(x)
#     @test ∇ℳ[i] ≈ ∂(ℳ2, x, i, 0.0001) atol=0.1
#   end
#   ∇ℳ = ReverseDiff.gradient(ℳ2, x)
#   for i ∈ 1:length(x)
#     @test ∇ℳ[i] ≈ ∂(ℳ2, x, i, 0.0001) atol=0.1
#   end
#   # FIXME: Zygote.gradient fails due to mutation in impulseresponse() and no adjoint for TriangularToeplitz

# end
