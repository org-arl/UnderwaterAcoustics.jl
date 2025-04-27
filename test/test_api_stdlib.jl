using TestItems

@testsnippet PekerisSetup begin
  env = UnderwaterEnvironment(
    bathymetry = 20.0u"m",
    temperature = 27.0u"°C",
    salinity = 35.0u"ppt"
  )
end

@testitem "models" setup=[PekerisSetup] begin
  m = @inferred models()
  @test m isa Vector{Type{<:UnderwaterAcoustics.AbstractPropagationModel}}
  @test length(m) == 2
  @test PekerisRayTracer ∈ m
  @test PekerisModeSolver ∈ m
  m = @inferred models(env)
  @test m isa Vector{Type{<:UnderwaterAcoustics.AbstractPropagationModel}}
  @test length(m) == 2
  @test PekerisRayTracer ∈ m
  @test PekerisModeSolver ∈ m
  env = UnderwaterEnvironment(soundspeed = SampledField([1500, 1490, 1520]; z=0:-10:-20, interp=:cubic))
  m = @inferred models(env)
  @test m isa Vector{Type{<:UnderwaterAcoustics.AbstractPropagationModel}}
  @test length(m) == 0
end

@testitem "env" setup=[PekerisSetup] begin
  for e ∈ (env, UnderwaterEnvironment(bathymetry=20.0u"m"))
    @test e isa UnderwaterEnvironment
    @test ! @inferred is_range_dependent(e)
    @test @inferred is_isovelocity(e)
    @test e.bathymetry == 20.0
    @test e.altimetry == 0.0
    @test e.temperature == 27.0
    @test e.salinity == 35.0
    @test e.soundspeed ≈ 1539.0 atol=0.1
    @test e.density ≈ 1022.7 atol=0.1
    @test e.surface == PressureReleaseBoundary
    @test e.seabed == RigidBoundary
    @test UnderwaterAcoustics.env_type(env) == Float64
  end
end

@testitem "src" begin
  import UnderwaterAcoustics: distance
  src = @inferred AcousticSource(-5.0, 1000.0; spl=180.0)
  @test src isa UnderwaterAcoustics.AbstractAcousticSource
  @test @inferred(location(src)) == (x=0.0, y=0.0, z=-5.0)
  @test @inferred(distance(location(src), (x=0.0, y=0.0, z=-5.0))) == 0.0
  @test @inferred(distance(location(src), (x=1.0, y=1.0, z=-4.0))) == sqrt(3)
  @test @inferred(frequency(src)) == 1000.0
  @test @inferred(spl(src)) == 180.0
  src = @inferred AcousticSource((10.0, -5.0), 1000.0)
  @test @inferred(location(src)) == (x=10.0, y=0.0, z=-5.0)
  @test @inferred(frequency(src)) == 1000.0
  @test @inferred(spl(src)) == 0.0
  src = @inferred AcousticSource((10.0, 20.0, -5.0), 1000.0)
  @test @inferred(location(src)) == (x=10.0, y=20.0, z=-5.0)
  @test @inferred(frequency(src)) == 1000.0
  @test @inferred(spl(src)) == 0.0
  src = @inferred AcousticSource(-5u"m", 1u"kHz"; spl=180u"dB")
  @test @inferred(location(src)) == (x=0.0, y=0.0, z=-5.0)
  @test @inferred(frequency(src)) == 1000.0
  @test @inferred(spl(src)) == 180.0
  src = @inferred AcousticSource((10u"m", -5u"m"), 1000.0)
  @test @inferred(location(src)) == (x=10.0, y=0.0, z=-5.0)
  @test @inferred(frequency(src)) == 1000.0
  @test @inferred(spl(src)) == 0.0
  src = @inferred AcousticSource((10u"m", 20u"m", -5u"m"), 1000.0)
  @test @inferred(location(src)) == (x=10.0, y=20.0, z=-5.0)
  @test @inferred(frequency(src)) == 1000.0
  src = @inferred AcousticSource((x=10u"m", y=20u"m", z=-5u"m"), 1u"kHz")
  @test @inferred(location(src)) == (x=10.0, y=20.0, z=-5.0)
  @test @inferred(frequency(src)) == 1000.0
  src = @inferred AcousticSource(nothing, 1000.0)
  @test @inferred(location(src)) === nothing
  @test @inferred(frequency(src)) == 1000.0
  src = @inferred AcousticSource(missing, 1000.0)
  @test @inferred(location(src)) === missing
  @test @inferred(frequency(src)) == 1000.0
end

@testitem "rcv" begin
  rcv = @inferred AcousticReceiver(10.0, 20.0, 0.0)
  @test rcv isa UnderwaterAcoustics.AbstractAcousticReceiver
  @test @inferred(location(rcv)) == (x=10.0, y=20.0, z=0.0)
  rcv = @inferred AcousticReceiver((10.0, 20.0, 0.0))
  @test @inferred(location(rcv)) == (x=10.0, y=20.0, z=0.0)
  rcv = @inferred AcousticReceiver(10.0, -20.0)
  @test @inferred(location(rcv)) == (x=10.0, y=0.0, z=-20.0)
  rcv = @inferred AcousticReceiver((10.0, -20.0))
  @test @inferred(location(rcv)) == (x=10.0, y=0.0, z=-20.0)
  rcv = @inferred AcousticReceiver(-20.0)
  @test @inferred(location(rcv)) == (x=0.0, y=0.0, z=-20.0)
  rcv = @inferred AcousticReceiver(10u"m", 20u"m", 0u"m")
  @test @inferred(location(rcv)) == (x=10.0, y=20.0, z=0.0)
  rcv = @inferred AcousticReceiver((10u"m", 20u"m", 0u"m"))
  @test @inferred(location(rcv)) == (x=10.0, y=20.0, z=0.0)
  rcv = @inferred AcousticReceiver(10u"m", -20u"m")
  @test @inferred(location(rcv)) == (x=10.0, y=0.0, z=-20.0)
  rcv = @inferred AcousticReceiver((10u"m", -20u"m"))
  @test @inferred(location(rcv)) == (x=10.0, y=0.0, z=-20.0)
  rcv = @inferred AcousticReceiver(-20u"m")
  @test @inferred(location(rcv)) == (x=0.0, y=0.0, z=-20.0)
end

@testitem "rcv-grid" begin
  rcv = @inferred AcousticReceiverGrid2D(0.0:10.0:1000.0, -20.0:0.0)
  @test rcv isa AbstractMatrix{AcousticReceiver}
  @test size(rcv) == (101, 21)
  @test @inferred(location(rcv[1, 1])) == (x=0.0, y=0.0, z=-20.0)
  @test @inferred(location(rcv[101, 21])) == (x=1000.0, y=0.0, z=0.0)
  rcv = @inferred AcousticReceiverGrid2D(0u"m":10u"m":1u"km", -20u"m":1u"m":0u"m")
  @test size(rcv) == (101, 21)
  @test @inferred(location(rcv[1, 1])) == (x=0.0, y=0.0, z=-20.0)
  @test @inferred(location(rcv[101, 21])) == (x=1000.0, y=0.0, z=0.0)
  rcv = @inferred AcousticReceiverGrid3D(0.0:10.0:1000.0, -100:100, -20.0:0.0)
  @test rcv isa AbstractArray{AcousticReceiver,3}
  @test size(rcv) == (101, 201, 21)
  @test @inferred(location(rcv[1, 1, 1])) == (x=0.0, y=-100.0, z=-20.0)
  @test @inferred(location(rcv[101, 201, 21])) == (x=1000.0, y=100.0, z=0.0)
  rcv = @inferred AcousticReceiverGrid3D(0u"m":10u"m":1000u"m", -100u"m":1u"m":100u"m", -20u"m":1u"m":0u"m")
  @test size(rcv) == (101, 201, 21)
  @test @inferred(location(rcv[1, 1, 1])) == (x=0.0, y=-100.0, z=-20.0)
  @test @inferred(location(rcv[101, 201, 21])) == (x=1000.0, y=100.0, z=0.0)
end

@testitem "boundary" begin
  @test FluidBoundary <: UnderwaterAcoustics.AbstractAcousticBoundary
  @test PressureReleaseBoundary isa FluidBoundary
  @test RigidBoundary isa FluidBoundary
  for sb ∈ (Rock, Pebbles, SandyGravel, CoarseSand, MediumSand, FineSand, VeryFineSand, ClayeySand, CoarseSilt, SandySilt, Silt, FineSilt, SandyClay, SiltyClay, Clay)
    @test sb isa FluidBoundary
    rc1 = @inferred reflection_coef(sb, 1000.0, 0.1, 1023.0, 1540.0)
    rc2 = @inferred reflection_coef(0.1, sb.ρ/1023.0, sb.c/1540.0, sb.δ)
    @test rc1 == rc2
  end
  @test WindySurface <: UnderwaterAcoustics.AbstractAcousticBoundary
  for ss ∈ (SeaState0, SeaState1, SeaState2, SeaState3, SeaState4, SeaState5, SeaState6, SeaState7, SeaState8, SeaState9)
    @test ss isa WindySurface
    rc1 = @inferred reflection_coef(ss, 1000.0, 0.1, 1023.0, 1540.0)
    rc2 = @inferred surface_reflection_coef(ss.windspeed, 1000.0, 0.1)
    @test rc1 == rc2
  end
end

@testitem "field" begin
  fld = @inferred SampledField([0.0, 10.0, 0.0]; x=[0.0, 10.0, 20.0])
  @test @inferred is_range_dependent(fld)
  @test @inferred(fld(0.0)) == 0.0
  @test @inferred(fld(5.0)) == 5.0
  @test @inferred(fld(10.0)) == 10.0
  @test @inferred(fld(15.0)) == 5.0
  @test @inferred(fld(20.0)) == 0.0
  @test @inferred(fld((x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(fld((x=5.0, y=0.0, z=0.0))) == 5.0
  @test @inferred(fld((x=10.0, y=0.0, z=0.0))) == 10.0
  @test @inferred(fld((x=15.0, y=0.0, z=0.0))) == 5.0
  @test @inferred(fld((x=20.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(fld((x=10.0, y=0.0, z=5.0))) == 10.0
  @test @inferred(fld((x=10.0, y=5.0, z=0.0))) == 10.0
  fld = @inferred SampledField([0.0, 10.0, 0.0]; z=[0.0, 10.0, 20.0])
  @test @inferred !is_range_dependent(fld)
  @test @inferred(fld(0.0)) == 0.0
  @test @inferred(fld(5.0)) == 5.0
  @test @inferred(fld(10.0)) == 10.0
  @test @inferred(fld(15.0)) == 5.0
  @test @inferred(fld(20.0)) == 0.0
  @test @inferred(fld((x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(fld((x=0.0, y=0.0, z=5.0))) == 5.0
  @test @inferred(fld((x=0.0, y=0.0, z=10.0))) == 10.0
  @test @inferred(fld((x=0.0, y=0.0, z=15.0))) == 5.0
  @test @inferred(fld((x=0.0, y=0.0, z=20.0))) == 0.0
  @test @inferred(fld((x=5.0, y=0.0, z=10.0))) == 10.0
  @test @inferred(fld((x=0.0, y=5.0, z=10.0))) == 10.0
  fld = @inferred SampledField([0.0 1.0; 1.0 2.0]; x=[0.0, 1.0], y=[0.0, 1.0])
  @test @inferred is_range_dependent(fld)
  @test @inferred(fld(0.0, 0.0)) == 0.0
  @test @inferred(fld(0.5, 0.5)) == 1.0
  @test @inferred(fld(1.0, 1.0)) == 2.0
  @test @inferred(fld((x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(fld((x=0.5, y=0.5, z=0.0))) == 1.0
  @test @inferred(fld((x=1.0, y=1.0, z=0.0))) == 2.0
  @test @inferred(fld((x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(fld((x=1.0, y=0.0, z=1.0))) == 1.0
  fld = SampledField([0.0; 1.0;; 1.0; 2.0;;; 0.0; 1.0;; 1.0; 2.0]; x=[0.0, 1.0], y=[0.0, 1.0], z=[0.0, 1.0])
  @test @inferred is_range_dependent(fld)
  @test @inferred(fld(0.0, 0.0, 0.0)) == 0.0
  @test @inferred(fld(0.5, 0.5, 0.5)) == 1.0
  @test @inferred(fld(1.0, 1.0, 1.0)) == 2.0
  @test @inferred(fld((x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(fld((x=0.5, y=0.5, z=0.5))) == 1.0
  @test @inferred(fld((x=1.0, y=1.0, z=1.0))) == 2.0
end

@testitem "noise" begin
  @test WhiteGaussianNoise <: UnderwaterAcoustics.AbstractNoiseModel
  @test RedGaussianNoise <: UnderwaterAcoustics.AbstractNoiseModel
  @test WhiteGaussianNoise(1.0) == WhiteGaussianNoise(0.0, 2.0)
  @test WhiteGaussianNoise(1f0) == WhiteGaussianNoise(0f0, 2f0)
  x = @inferred rand(WhiteGaussianNoise(1.0), 100; fs=1)
  @test x isa AbstractVector{Float64}
  @test length(x) == 100
  x = @inferred rand(WhiteGaussianNoise(1f0), 100; fs=1)
  @test x isa AbstractVector{Float32}
  @test length(x) == 100
  x = @inferred rand(WhiteGaussianNoise(0.0, 2.0), 100; fs=1)
  @test x isa AbstractVector{Float64}
  @test length(x) == 100
  x = @inferred rand(WhiteGaussianNoise(0f0, 2f0), 100; fs=1)
  @test x isa AbstractVector{Float32}
  @test length(x) == 100
  x = @inferred rand(RedGaussianNoise(1.0), 100; fs=1)
  @test x isa AbstractVector{Float64}
  @test length(x) == 100
  x = @inferred rand(RedGaussianNoise(1f0), 100; fs=1)
  @test x isa AbstractVector{Float32}
  @test length(x) == 100
  x = @inferred rand(WhiteGaussianNoise(1.0), 100, 2; fs=1)
  @test x isa AbstractMatrix{Float64}
  @test size(x) == (100, 2)
  x = @inferred rand(WhiteGaussianNoise(1f0), 100, 2; fs=1)
  @test x isa AbstractMatrix{Float32}
  @test size(x) == (100, 2)
  x = @inferred rand(WhiteGaussianNoise(0.0, 2.0), 100, 2; fs=1)
  @test x isa AbstractMatrix{Float64}
  @test size(x) == (100, 2)
  x = @inferred rand(WhiteGaussianNoise(0f0, 2f0), 100, 2; fs=1)
  @test x isa AbstractMatrix{Float32}
  @test size(x) == (100, 2)
  x = @inferred rand(RedGaussianNoise(1.0), 100, 2; fs=1)
  @test x isa AbstractMatrix{Float64}
  @test size(x) == (100, 2)
  x = @inferred rand(RedGaussianNoise(1f0), 100, 2; fs=1)
  @test x isa AbstractMatrix{Float32}
  @test size(x) == (100, 2)
  x = @inferred rand(WhiteGaussianNoise(1.0), (100, 2); fs=1)
  @test x isa AbstractMatrix{Float64}
  @test size(x) == (100, 2)
  x = @inferred rand(WhiteGaussianNoise(1f0), (100, 2); fs=1)
  @test x isa AbstractMatrix{Float32}
  @test size(x) == (100, 2)
  x = @inferred rand(WhiteGaussianNoise(0.0, 2.0), (100, 2); fs=1)
  @test x isa AbstractMatrix{Float64}
  @test size(x) == (100, 2)
  x = @inferred rand(WhiteGaussianNoise(0f0, 2f0), (100, 2); fs=1)
  @test x isa AbstractMatrix{Float32}
  @test size(x) == (100, 2)
  x = @inferred rand(RedGaussianNoise(1.0), (100, 2); fs=1)
  @test x isa AbstractMatrix{Float64}
  @test size(x) == (100, 2)
  x = @inferred rand(RedGaussianNoise(1f0), (100, 2); fs=1)
  @test x isa AbstractMatrix{Float32}
  @test size(x) == (100, 2)
end
