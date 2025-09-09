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
  @test length(m) == 3
  @test PekerisRayTracer ∈ m
  @test PekerisModeSolver ∈ m
  @test AdiabaticExt ∈ m
  m = @inferred models(UnderwaterAcoustics.AbstractModePropagationModel)
  @test m isa Vector{Type{<:UnderwaterAcoustics.AbstractPropagationModel}}
  @test length(m) == 2
  @test PekerisModeSolver ∈ m
  @test AdiabaticExt ∈ m
end

@testitem "env" setup=[PekerisSetup] begin
  for e ∈ (env, UnderwaterEnvironment(bathymetry=20.0u"m"))
    @test e isa UnderwaterEnvironment
    @test startswith(sprint(show, e), "UnderwaterEnvironment")
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

@testitem "channel" begin
  using Statistics
  env = UnderwaterEnvironment(seabed=FluidBoundary(water_density(), 1500.0, 0.0), soundspeed=1500.0)
  pm = @inferred PekerisRayTracer(env)
  tx = @inferred AcousticSource((x=0u"m", z=-5u"m"), 10u"kHz"; spl=170)
  rx = @inferred AcousticReceiver(100u"m", -10u"m")
  ch = @inferred channel(pm, tx, rx, 192000)
  @test ch isa UnderwaterAcoustics.SampledPassbandChannel
  @test startswith(sprint(show, ch), "SampledPassbandChannel")
  @test ch.fs == 192000
  @test ch.noise === nothing
  x = vcat(1.0, 0.5, zeros(98))
  y1 = @inferred transmit(ch, x)
  @test length(y1) ≥ length(x)
  @test all(y1[1:3] .> 1)
  @test all(y1[129:131] .< -1)
  y1[1:3] .= 0
  y1[129:131] .= 0
  @test all(abs.(y1) .< 1)
  ch = @inferred channel(pm, tx, rx, 192000; noise=RedGaussianNoise(0.5e6))
  @test ch isa UnderwaterAcoustics.SampledPassbandChannel
  @test ch.fs == 192000
  @test ch.noise isa RedGaussianNoise
  y2 = @inferred transmit(ch, x)
  @test length(y2) ≥ length(x)
  @test std(y1) < std(y2)
end

@testitem "src" begin
  import UnderwaterAcoustics: distance
  src = @inferred AcousticSource(-5.0, 1000.0; spl=180.0)
  @test src isa UnderwaterAcoustics.AbstractAcousticSource
  @test startswith(sprint(show, src), "TX")
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
  @test startswith(sprint(show, rcv), "RX")
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
  @test rcv isa AbstractMatrix{<:AcousticReceiver}
  @test size(rcv) == (101, 21)
  @test @inferred(location(rcv[1, 1])) == (x=0.0, y=0.0, z=-20.0)
  @test @inferred(location(rcv[101, 21])) == (x=1000.0, y=0.0, z=0.0)
  rcv = @inferred AcousticReceiverGrid2D(0u"m":10u"m":1u"km", -20u"m":1u"m":0u"m")
  @test size(rcv) == (101, 21)
  @test @inferred(location(rcv[1, 1])) == (x=0.0, y=0.0, z=-20.0)
  @test @inferred(location(rcv[101, 21])) == (x=1000.0, y=0.0, z=0.0)
  rcv = @inferred AcousticReceiverGrid3D(0.0:10.0:1000.0, -100:100, -20.0:0.0)
  @test rcv isa AbstractArray{<:AcousticReceiver,3}
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
  bc = FluidBoundary(1200, 1500)
  @test startswith(sprint(show, bc), "FluidBoundary")
  @test bc isa FluidBoundary
  @test bc == FluidBoundary(ρ=1200, c=1500)
  @test bc == FluidBoundary(1.2u"g/cm^3", 1500u"m/s")
  @test FluidBoundary(1200, 1500, 0.1) == FluidBoundary(ρ=1200, c=1500, δ=0.1)
  @test FluidBoundary(1200, 1500, 0.1, 0.2) == FluidBoundary(ρ=1200, c=1500, δ=0.1, σ=0.2)
  @test PressureReleaseBoundary isa FluidBoundary
  @test RigidBoundary isa FluidBoundary
  @test sprint(show, PressureReleaseBoundary) == "PressureReleaseBoundary"
  @test sprint(show, RigidBoundary) == "RigidBoundary"
  for sb ∈ (Rock, Pebbles, SandyGravel, VeryCoarseSand, MuddySandyGravel, CoarseSand, GravellyMuddySand, MediumSand, MuddyGravel, FineSand,
            MuddySand, VeryFineSand, ClayeySand, CoarseSilt, SandySilt, MediumSilt, SandyMud, FineSilt, SandyClay, VeryFineSilt, SiltyClay, Clay)
    @test sb isa FluidBoundary
    rc1 = @inferred reflection_coef(sb, 1000.0, 0.1, 1023.0, 1540.0)
    rc2 = @inferred reflection_coef(0.1, sb.ρ/1023.0, sb.c/1540.0, sb.δ)
    @test rc1 == rc2
  end
  @test ElasticBoundary <: UnderwaterAcoustics.AbstractAcousticBoundary
  bc = ElasticBoundary(1200, 1500, 500)
  @test startswith(sprint(show, bc), "ElasticBoundary")
  @test bc isa ElasticBoundary
  @test bc == ElasticBoundary(ρ=1200, cₚ=1500, cₛ=500)
  @test bc == ElasticBoundary(1.2u"g/cm^3", 1500u"m/s", 500u"m/s")
  @test ElasticBoundary(1200, 1500, 500, 0.1, 0.11) == ElasticBoundary(ρ=1200, cₚ=1500, cₛ=500, δₚ=0.1, δₛ=0.11)
  @test ElasticBoundary(1200, 1500, 500, 0.1, 0, 0.2) == ElasticBoundary(ρ=1200, cₚ=1500, cₛ=500, δₚ=0.1, σ=0.2)
  @test ElasticBoundary(FluidBoundary(1200, 1560, 0.1, 0.2)) == ElasticBoundary(ρ=1200, cₚ=1560, cₛ=288.72, δₚ=0.1, σ=0.2)
  @test ElasticBoundary(FluidBoundary(1200, 1560, 0.1, 0.2), 0.3) == ElasticBoundary(ρ=1200, cₚ=1560, cₛ=288.72, δₚ=0.1, δₛ=0.3, σ=0.2)
  @test ElasticBoundary(FluidBoundary(1200, 1560, 0.1, 0.2), 500, 0.3) == ElasticBoundary(ρ=1200, cₚ=1560, cₛ=500, δₚ=0.1, δₛ=0.3, σ=0.2)
  rc1a = @inferred reflection_coef(bc, 1000.0, 0.1, 1023.0, 1540.0)
  rc2a = @inferred reflection_coef(0.1, bc.ρ/1023.0, bc.cₚ/1540.0, bc.δₚ)
  @test rc1a == rc2a
  @test MultilayerElasticBoundary <: UnderwaterAcoustics.AbstractAcousticBoundary
  bc = MultilayerElasticBoundary([
    (5.2, 1300, 1700, 100, 0.1, 0.2, 0.3),
    (5.2, 1300, 1700, 100, 0.1, 0.2),
    (5.2, 1300, 1700, 100, 0.1),
    (Inf, 2000, 2500, 500)
  ])
  @test startswith(sprint(show, bc), "MultilayerElasticBoundary")
  @test startswith(sprint(show, MIME"text/plain"(), bc), "MultilayerElasticBoundary")
  @test length(bc.layers) == 4
  @test bc.layers == MultilayerElasticBoundary([
    (h = 5.2, ρ = 1300, cₚ = 1700, cₛ = 100, δₚ = 0.1, δₛ = 0.2, σ = 0.3),
    (h = 5.2, ρ = 1300, cₚ = 1700, cₛ = 100, δₚ = 0.1, δₛ = 0.2),
    (h = 5.2, ρ = 1300, cₚ = 1700, cₛ = 100, δₚ = 0.1),
    (h = Inf, ρ = 2000, cₚ = 2500, cₛ = 500)
  ]).layers
  @test bc.layers == MultilayerElasticBoundary([
    (5.2u"m", 1.3u"g/cm^3", 1700u"m/s", 100u"m/s", 0.1, 0.2, 0.3),
    (5.2u"m", 1.3u"g/cm^3", 1700u"m/s", 100u"m/s", 0.1, 0.2),
    (5.2u"m", 1.3u"g/cm^3", 1700u"m/s", 100u"m/s", 0.1),
    (Inf, 2u"g/cm^3", 2500u"m/s", 500u"m/s")
  ]).layers
  erock = @inferred ElasticBoundary(Rock, 0.1)
  bc1 = MultilayerElasticBoundary([(40, FineSand), (Inf, erock)])
  @test length(bc1.layers) == 2
  @test bc1.layers[1] == (h=40, ρ=FineSand.ρ, cₚ=FineSand.c, cₛ=0.0, δₚ=FineSand.δ, δₛ=0.0, σ=FineSand.σ)
  @test bc1.layers[2] == (h=Inf, ρ=erock.ρ, cₚ=erock.cₚ, cₛ=erock.cₛ, δₚ=erock.δₚ, δₛ=erock.δₛ, σ=erock.σ)
  rc1a = @inferred reflection_coef(bc, 1000.0, 0.1, 1023.0, 1540.0)
  rc2a = @inferred reflection_coef(0.1, bc.layers[1].ρ/1023.0, bc.layers[1].cₚ/1540.0, bc.layers[1].δₚ)
  @test rc1a == rc2a
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
  @test fld isa UnderwaterAcoustics.SampledFieldX
  @test startswith(sprint(show, fld), "SampledField")
  @test @inferred is_range_dependent(fld)
  @test @inferred(fld(0.0)) == 0.0
  @test @inferred(fld(5.0)) == 5.0
  @test @inferred(fld(10.0)) == 10.0
  @test @inferred(fld(15.0)) == 5.0
  @test @inferred(fld(20.0)) == 0.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(value(fld, (x=5.0, y=0.0, z=0.0))) == 5.0
  @test @inferred(value(fld, (x=10.0, y=0.0, z=0.0))) == 10.0
  @test @inferred(value(fld, (x=15.0, y=0.0, z=0.0))) == 5.0
  @test @inferred(value(fld, (x=20.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(value(fld, (x=10.0, y=0.0, z=5.0))) == 10.0
  @test @inferred(value(fld, (x=10.0, y=5.0, z=0.0))) == 10.0
  fld = @inferred SampledField([0.0, 10.0, 0.0]; z=[0.0, 10.0, 20.0])
  @test fld isa UnderwaterAcoustics.SampledFieldZ
  @test startswith(sprint(show, fld), "SampledField")
  @test @inferred !is_range_dependent(fld)
  @test @inferred(fld(0.0)) == 0.0
  @test @inferred(fld(5.0)) == 5.0
  @test @inferred(fld(10.0)) == 10.0
  @test @inferred(fld(15.0)) == 5.0
  @test @inferred(fld(20.0)) == 0.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=5.0))) == 5.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=10.0))) == 10.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=15.0))) == 5.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=20.0))) == 0.0
  @test @inferred(value(fld, (x=5.0, y=0.0, z=10.0))) == 10.0
  @test @inferred(value(fld, (x=0.0, y=5.0, z=10.0))) == 10.0
  fld = @inferred SampledField([0.0 1.0; 1.0 2.0]; x=[0.0, 1.0], y=[0.0, 1.0])
  @test fld isa UnderwaterAcoustics.SampledFieldXY
  @test startswith(sprint(show, fld), "SampledField")
  @test @inferred is_range_dependent(fld)
  @test @inferred(fld(0.0, 0.0)) == 0.0
  @test @inferred(fld(0.5, 0.5)) == 1.0
  @test @inferred(fld(1.0, 1.0)) == 2.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(value(fld, (x=0.5, y=0.5, z=0.0))) == 1.0
  @test @inferred(value(fld, (x=1.0, y=1.0, z=0.0))) == 2.0
  @test @inferred(value(fld, (x=1.0, y=0.0, z=1.0))) == 1.0
  fld = @inferred SampledField([0.0 1.0; 1.0 2.0]; x=[0.0, 1.0], z=[0.0, 1.0])
  @test fld isa UnderwaterAcoustics.SampledFieldXZ
  @test startswith(sprint(show, fld), "SampledField")
  @test @inferred is_range_dependent(fld)
  @test @inferred(fld(0.0, 0.0)) == 0.0
  @test @inferred(fld(0.5, 0.5)) == 1.0
  @test @inferred(fld(1.0, 1.0)) == 2.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(value(fld, (x=0.5, y=0.0, z=0.5))) == 1.0
  @test @inferred(value(fld, (x=1.0, y=0.0, z=1.0))) == 2.0
  @test @inferred(value(fld, (x=1.0, y=1.0, z=0.0))) == 1.0
  fld = SampledField([0.0; 1.0;; 1.0; 2.0;;; 0.0; 1.0;; 1.0; 2.0]; x=[0.0, 1.0], y=[0.0, 1.0], z=[0.0, 1.0])
  @test fld isa UnderwaterAcoustics.SampledFieldXYZ
  @test startswith(sprint(show, fld), "SampledField")
  @test @inferred is_range_dependent(fld)
  @test @inferred(fld(0.0, 0.0, 0.0)) == 0.0
  @test @inferred(fld(0.5, 0.5, 0.5)) == 1.0
  @test @inferred(fld(1.0, 1.0, 1.0)) == 2.0
  @test @inferred(value(fld, (x=0.0, y=0.0, z=0.0))) == 0.0
  @test @inferred(value(fld, (x=0.5, y=0.5, z=0.5))) == 1.0
  @test @inferred(value(fld, (x=1.0, y=1.0, z=1.0))) == 2.0
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
