using TestItems

@testitem "pekeris-rays+ir" begin
  env = UnderwaterEnvironment(bathymetry = 20.0u"m", temperature = 27.0u"°C", salinity = 35.0u"ppt", seabed = SandySilt)
  pm = PekerisRayTracer(env, 3)
  tx = AcousticSource((x=0.0, z=-5.0), 1000.0)
  rx = AcousticReceiver((x=100.0, z=-10.0))
  arr = arrivals(pm, tx, rx)
  @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(arr) == 7
  @test arr[1].t ≈ 0.0650 atol=0.0001
  @test arr[2].t ≈ 0.0657 atol=0.0001
  @test arr[3].t ≈ 0.0670 atol=0.0001
  @test all([arr[j].t > arr[j-1].t for j ∈ 2:7])
  @test abs(arr[1].ϕ) ≈ 0.01 atol=0.001
  @test real(arr[2].ϕ) < 0.0
  @test imag(arr[2].ϕ) ≈ 0.0
  @test all([abs(arr[j].ϕ) < abs(arr[j-1].ϕ) for j ∈ 2:7])
  @test [(arr[j].ns, arr[j].nb) for j ∈ 1:7] == [(0,0), (1,0), (0,1), (1,1), (1,1), (2,1), (1,2)]
  @test all([abs(arr[j].θᵣ) == abs(arr[j].θₛ) for j ∈ 1:7])
  @test all([arr[j].path[1] == (x=0.0, y=0.0, z=-5.0) for j ∈ 1:7])
  @test all([arr[j].path[end] == (x=100.0, y=0.0, z=-10.0) for j ∈ 1:7])
  @test all([length(arr[j].path) == arr[j].ns + arr[j].nb + 2 for j ∈ 1:7])
  ir1 = impulse_response(pm, tx, rx, 10000.0; abstime=false)
  ir2 = impulse_response(pm, tx, rx, 10000.0; abstime=true)
  @test length(ir2) ≈ length(ir1) + round(Int, 10000.0 * arr[1].t) atol=1
  @test length(ir2) ≈ round(Int, 10000.0 * arr[end].t) + 1 atol=1
  @test all(abs.(ir1[[round(Int, 10000.0 * (arr[j].t - arr[1].t)) + 1 for j ∈ 1:7]]) .> 0)
  @test all(abs.(ir2[[round(Int, 10000.0 * arr[j].t) + 1 for j ∈ 1:7]]) .> 0)
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=256, abstime=true)) == 256
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=64, abstime=true)) == 64
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=1024, abstime=false)) == 1024
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=700, abstime=false)) == 700
end

@testitem "pekeris-field+tl" begin
  env = UnderwaterEnvironment(bathymetry = 20.0u"m", soundspeed = 1500.0u"m/s",
    seabed = FluidBoundary(water_density(), 1500.0u"m/s", 0.0))  # non-reflecting
  pm = PekerisRayTracer(env, 1)
  d = (√1209.0)/4.0
  x = acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test x isa Complex
  @test abs(x) ≈ 0.0 atol=0.0002
  x′ = acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test x′ isa Complex
  @test imag(x′) == 0.0
  @test abs(x′) > 1/100.0
  d = (√2409.0)/8.0
  x = acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test abs(x) > abs(x′)
  y = transmission_loss(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x = acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test abs(x) ≈ abs(x′) atol=0.0001
  y = transmission_loss(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  tx = AcousticSource((x=0.0, z=-5.0), 1000.0)
  x1 = acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-5.0)))
  x2 = acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-10.0)))
  x3 = acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-15.0)))
  x = acoustic_field(pm, tx, [AcousticReceiver((x=100.0, z=-d)) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = acoustic_field(pm, tx, AcousticReceiverGrid2D(100.0:100.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = acoustic_field(pm, tx, AcousticReceiverGrid2D(100.0:10.0:120.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]
  x1 = transmission_loss(pm, tx, AcousticReceiver(100.0, -5.0))
  x2 = transmission_loss(pm, tx, AcousticReceiver(100.0, -10.0))
  x3 = transmission_loss(pm, tx, AcousticReceiver(100.0, -15.0))
  x = transmission_loss(pm, tx, [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = transmission_loss(pm, tx, AcousticReceiverGrid2D(100.0:100.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = transmission_loss(pm, tx, AcousticReceiverGrid2D(100.0:10.0:120.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]
end

@testitem "pekeris-∂" begin
  using DifferentiationInterface
  import ForwardDiff, FiniteDifferences, Zygote, Mooncake, Enzyme
  fd = AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(5, 1))
  function ℳ₁((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = SandySilt)
    pm = PekerisRayTracer(env, 3)
    transmission_loss(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)))
  end
  function ℳ₂((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = SandySilt)
    pm = PekerisRayTracer(env, 3)
    ir = impulse_response(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)), 10000.0; abstime=true)
    sum(abs2, samples(ir))
  end
  x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
  ∇ℳ₁ = gradient(ℳ₁, fd, x)
  ∇ℳ₂ = gradient(ℳ₂, fd, x)
  @test gradient(ℳ₁, AutoForwardDiff(), x) ≈ ∇ℳ₁
  @test gradient(ℳ₁, AutoZygote(), x) ≈ ∇ℳ₁
  @test gradient(ℳ₁, AutoMooncake(config=nothing), x) ≈ ∇ℳ₁
  @test gradient(ℳ₁, AutoEnzyme(), x) ≈ ∇ℳ₁
  @test gradient(ℳ₂, AutoForwardDiff(), x) ≈ ∇ℳ₂
  # FIXME: Zygote fails due to mutation in impulse_response()
  # @test gradient(ℳ₂, AutoZygote(), x) ≈ ∇ℳ₂
  @test gradient(ℳ₂, AutoMooncake(config=nothing), x) ≈ ∇ℳ₂
  @test gradient(ℳ₂, AutoEnzyme(), x) ≈ ∇ℳ₂
  x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
  ∇ℳ₁ = gradient(ℳ₁, fd, x)
  ∇ℳ₂ = gradient(ℳ₂, fd, x)
  @test gradient(ℳ₁, AutoForwardDiff(), x) ≈ ∇ℳ₁
  @test gradient(ℳ₁, AutoZygote(), x) ≈ ∇ℳ₁
  @test gradient(ℳ₁, AutoMooncake(config=nothing), x) ≈ ∇ℳ₁
  @test gradient(ℳ₁, AutoEnzyme(), x) ≈ ∇ℳ₁
  @test gradient(ℳ₂, AutoForwardDiff(), x) ≈ ∇ℳ₂
  # FIXME: Zygote fails due to mutation in impulse_response()
  # @test gradient(ℳ₂, AutoZygote(), x) ≈ ∇ℳ₂
  @test gradient(ℳ₂, AutoMooncake(config=nothing), x) ≈ ∇ℳ₂
  @test gradient(ℳ₂, AutoEnzyme(), x) ≈ ∇ℳ₂
end
