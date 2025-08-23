using TestItems

@testitem "pekeris-rays-arrivals+ir" begin
  env = @inferred UnderwaterEnvironment(bathymetry = 20.0u"m", temperature = 27.0u"°C", salinity = 35.0u"ppt", seabed = SandySilt)
  pm = @inferred PekerisRayTracer(env)
  @test startswith(sprint(show, pm), "PekerisRayTracer")
  tx = @inferred AcousticSource((x=0.0, z=-5.0), 1000.0)
  rx = @inferred AcousticReceiver((x=100.0, z=-10.0))
  arr = @inferred arrivals(pm, tx, rx)
  @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test startswith(sprint(show, arr[1]), "∠")
  @test length(arr) == 7
  @test arr[1].t == arr[1].time
  @test arr[1].ϕ == arr[1].phasor
  @test arr[1].ns == arr[1].surface_bounces
  @test arr[1].nb == arr[1].bottom_bounces
  @test arr[1].θₛ == arr[1].launch_angle
  @test arr[1].θᵣ == arr[1].arrival_angle
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
  ir1 = @inferred impulse_response(pm, tx, rx, 10000.0; abstime=false)
  ir2 = @inferred impulse_response(pm, tx, rx, 10000.0; abstime=true)
  @test ir1 isa AbstractVector{<:Complex}
  @test ir2 isa AbstractVector{<:Complex}
  @test length(ir2) ≈ length(ir1) + round(Int, 10000.0 * arr[1].t) - 10 atol=1
  @test length(ir2) ≈ round(Int, 10000.0 * arr[end].t) + 11 atol=1
  @test all(abs.(ir1[[round(Int, 10000.0 * (arr[j].t - arr[1].t + 0.001)) + 1 for j ∈ 1:7]]) .> 0)
  @test all(abs.(ir2[[round(Int, 10000.0 * arr[j].t) + 1 for j ∈ 1:7]]) .> 0)
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=256, abstime=true)) == 256
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=64, abstime=true)) == 64
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=1024, abstime=false)) == 1024
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=700, abstime=false)) == 700
end

@testitem "pekeris-rays-field+tl" begin
  env = @inferred UnderwaterEnvironment(bathymetry = 20.0u"m", soundspeed = 1500.0u"m/s",
    seabed = FluidBoundary(water_density(), 1500.0u"m/s", 0.0))  # non-reflecting
  pm = @inferred PekerisRayTracer(env; max_bounces=1)
  d = (√1209.0)/4.0
  x = @inferred acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test x isa Complex
  @test abs(x) ≈ 0.0 atol=0.0002
  x′ = @inferred acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test x′ isa Complex
  @test imag(x′) == 0.0
  @test abs(x′) > 1/100.0
  d = (√2409.0)/8.0
  x = @inferred acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test abs(x) > abs(x′)
  y = @inferred transmission_loss(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x = @inferred acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test abs(x) ≈ abs(x′) atol=0.0001
  y = @inferred transmission_loss(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  tx = @inferred AcousticSource((x=0.0, z=-5.0), 1000.0)
  x1 = @inferred acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-5.0)))
  x2 = @inferred acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-10.0)))
  x3 = @inferred acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-15.0)))
  x = @inferred acoustic_field(pm, tx, [AcousticReceiver((x=100.0, z=-d)) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = @inferred acoustic_field(pm, tx, AcousticReceiverGrid2D(100.0:100.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = @inferred acoustic_field(pm, tx, AcousticReceiverGrid2D(100.0:10.0:120.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]
  x1 = @inferred transmission_loss(pm, tx, AcousticReceiver(100.0, -5.0))
  x2 = @inferred transmission_loss(pm, tx, AcousticReceiver(100.0, -10.0))
  x3 = @inferred transmission_loss(pm, tx, AcousticReceiver(100.0, -15.0))
  x = @inferred transmission_loss(pm, tx, [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = @inferred transmission_loss(pm, tx, AcousticReceiverGrid2D(100.0:100.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = @inferred transmission_loss(pm, tx, AcousticReceiverGrid2D(100.0:10.0:120.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]
end

@testitem "pekeris-rays-∂" begin
  using DifferentiationInterface
  import ForwardDiff, FiniteDifferences, Zygote, Mooncake, Enzyme
  fd = AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(5, 1))
  function ℳ₁((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = SandySilt)
    pm = PekerisRayTracer(env)
    transmission_loss(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)))
  end
  function ℳ₂((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = SandySilt)
    pm = PekerisRayTracer(env)
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
  @test gradient(ℳ₂, AutoZygote(), x) ≈ ∇ℳ₂
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
  @test gradient(ℳ₂, AutoZygote(), x) ≈ ∇ℳ₂
  @test gradient(ℳ₂, AutoMooncake(config=nothing), x) ≈ ∇ℳ₂
  @test gradient(ℳ₂, AutoEnzyme(), x) ≈ ∇ℳ₂
end

@testitem "pekeris-modes-arrivals+tl" begin
  env = @inferred UnderwaterEnvironment(
    bathymetry = 5000,
    soundspeed = 1500,
    density = 1000,
    seabed = FluidBoundary(2000, 2000)
  )
  pm = @inferred PekerisModeSolver(env)
  @test startswith(sprint(show, pm), "PekerisModeSolver")
  tx = @inferred AcousticSource(0, -500, 10)
  rx = @inferred AcousticReceiver(200000, -2500)
  m = @inferred arrivals(pm, tx, rx)
  @test m isa Vector{<:UnderwaterAcoustics.ModeArrival}
  @test contains(sprint(show, m[1]), "mode")
  @test length(m) == 44
  @test m[1].m == m[1].mode
  @test m[1].kᵣ == m[1].hwavenumber
  @test m[1].v == m[1].group_velocity
  @test m[1].ψ == m[1].mode_function
  k = [0.04188332253060325,  0.04186958032482773,  0.04184666447436872,  0.041814556737327056,
       0.041773231610265325, 0.041722656355626477, 0.041662791026605805, 0.041593588482909706,
       0.041514994390511865, 0.04142694719862171,  0.04132937808752444,  0.04122221088162445,
       0.04110536192280678,  0.04097873990001411,  0.04084224563163022,  0.04069577179780036,
       0.040539202620165,    0.04037241348662042,  0.04019527051864795,  0.040007630078488614,
       0.03980933821300124,  0.03960023003045652,  0.039380129005813426, 0.03914884620922843,
       0.03890617945169683,  0.03865191234086701,  0.03838581323927118,  0.03810763411659509,
       0.037817109287352026, 0.03751395402576234,  0.037197863051343545, 0.03686850888271845,
       0.036525540065337404, 0.03616857929474754,  0.03579722148769253,  0.0354110319129428,
       0.035009544613159714, 0.03459226159966202,  0.03415865386378722,  0.03370816662868136,
       0.033240235114164754, 0.0327543299843113,   0.03225010949713772,  0.03172822169111466]
  @test [real(m1.kᵣ) for m1 ∈ m] ≈ k
  v = log(absorption(frequency(tx), 1.0, pm.S, pm.T, pm.h / 2, pm.pH))
  @test all([imag(m1.kᵣ) ≈ v for m1 ∈ m])
  v = [1499.8406693989984, 1499.3624730362858, 1498.5648007263292, 1497.4466468441854,
       1496.006626357596,  1494.2429954061495, 1492.153674826432,  1489.736274891665,
       1486.9881195724402, 1483.9062688059496, 1480.487537543422,  1476.7285106821262,
       1472.6255533344163, 1468.1748162053236, 1463.3722361172465, 1458.2135319228933,
       1452.6941961837783, 1446.8094830677962, 1440.5543929471391, 1433.9236541715757,
       1426.9117024670375, 1419.5126583808806, 1411.720303178339,  1403.5280536058474,
       1394.9289359954637, 1385.9155603166237, 1376.480095025709,  1366.614243982029,
       1356.30922739289,   1345.555769895444,  1334.3441007877284, 1322.6639746483686,
       1310.5047261965228, 1297.8553833330877, 1284.7048811936058, 1271.0424571833448,
       1256.8583844525185, 1242.1453751290705, 1226.9014118153054, 1211.1359455191348,
       1194.8852077517781, 1178.257775701711,  1161.6206969257646, 1147.1591430052147]
  @test [m1.v for m1 ∈ m] ≈ v
  rxs = @inferred AcousticReceiverGrid2D(200000:1000:220000, -2500)
  x = @inferred transmission_loss(pm, tx, rxs)
  y = [93.12294364780713, 86.50268953343257, 86.44886663114649,  84.0362459489373,
       88.77243070715846, 86.52199965747599, 97.59552929416957,  89.697653630194,
       92.2097226818825,  95.47903086283519, 88.48249620410911,  87.33789340605915,
      108.19712606519451, 93.39736419463111, 91.1187296811673,   91.65297876114658,
      87.12267265661998,  87.92951345976743, 90.81157040842854, 102.8825858777283,
      95.08527791225038]
  @test vec(x) ≈ y atol=0.1
  env = @inferred UnderwaterEnvironment(
    bathymetry = 5000,
    soundspeed = 1500,
    density = 1000,
    seabed = RigidBoundary
  )
  pm = @inferred PekerisModeSolver(env)
  m = @inferred arrivals(pm, tx, rx)
  @test length(m) == 67
  k = [0.041887, 0.041877, 0.041858, 0.04183, 0.041792, 0.041745, 0.041688, 0.041622, 0.041546, 0.04146]
  v = [1499.96, 1499.62, 1498.94, 1497.93, 1496.58, 1494.89, 1492.85, 1490.48, 1487.76, 1484.69]
  @test [m1.kᵣ for m1 ∈ m[1:10]] ≈ k atol=0.01
  @test [m1.v for m1 ∈ m[1:10]] ≈ v atol=0.01
  x = @inferred transmission_loss(pm, tx, rxs)
  y = [94.3836, 86.7388, 90.8998, 83.0766, 85.6271, 91.2697, 85.4745, 80.0225, 81.1577, 85.27,
       84.0323, 101.713, 88.4222, 88.8432, 82.2533, 87.7657, 90.7141, 91.4239, 91.3287, 87.339, 86.9366]
  @test vec(x) ≈ y atol=0.1
  env = @inferred UnderwaterEnvironment(
    bathymetry = 5000,
    soundspeed = 1500,
    density = 1000,
    seabed = PressureReleaseBoundary
  )
  pm = @inferred PekerisModeSolver(env)
  m = @inferred arrivals(pm, tx, rx)
  @test length(m) == 66
  k = [0.04188318, 0.0418690, 0.0418454, 0.0418124, 0.0417699, 0.04171791, 0.0416563, 0.0415852, 0.0415044, 0.041413]
  v = [1499.831, 1499.324, 1498.480, 1497.297, 1495.775, 1493.912, 1491.708, 1489.160, 1486.268, 1483.028]
  @test [m1.kᵣ for m1 ∈ m[1:10]] ≈ k atol=0.01
  @test [m1.v for m1 ∈ m[1:10]] ≈ v atol=0.01
  x = @inferred transmission_loss(pm, tx, rxs)
  y = [92.4121, 87.5373, 109.99, 84.9219, 91.5823, 94.6702, 89.9848, 90.8624, 92.803, 85.8212, 94.3161, 82.8287, 89.2737, 82.7979, 89.2751, 92.9164, 85.2361, 82.5535, 82.8422, 88.6535, 94.7227]
  @test vec(x) ≈ y atol=0.1
end

@testitem "pekeris-modes-ir" begin
  env = @inferred UnderwaterEnvironment(bathymetry=100, seabed=SandyGravel)
  pm = @inferred PekerisModeSolver(env)
  tx = @inferred AcousticSource(0, -30, 300)
  rx = @inferred AcousticReceiver(1000, -40)
  ir1 = @inferred impulse_response(pm, tx, rx, 8000; abstime=false)
  ir2 = @inferred impulse_response(pm, tx, rx, 8000; abstime=true)
  @test ir1 isa AbstractVector{<:Complex}
  @test ir2 isa AbstractVector{<:Complex}
  x = abs.(ir2)
  let i = argmax(x)
    @test abs((i - 1) / 8000 - hypot(1000, 10) / env.soundspeed) < 0.002
  end
  function get_Δt_first5(x)
    ndx = Int[]
    θ = maximum(x) / 10
    i = argmax(x)
    for j ∈ 1:5
      while i < length(x) && x[i+1] > x[i]
        i += 1
      end
      push!(ndx, i)
      while i < length(x) && x[i] > θ
        i += 1
      end
      while i < length(x) && x[i] < θ
        i += 1
      end
    end
    diff(ndx)
  end
  d51 = get_Δt_first5(abs.(ir1))
  d52 = get_Δt_first5(abs.(ir2))
  @test d51 == [179, 38, 54, 38]
  @test d52 == [179, 38, 54, 38]
  @test length(@inferred(impulse_response(pm, tx, rx, 8000; ntaps=4096))) == 4096
  @test length(@inferred(impulse_response(pm, tx, rx, 8000; ntaps=700))) == 700
end

@testitem "pekeris-modes-∂" begin
  using DifferentiationInterface
  import ForwardDiff, FiniteDifferences, Zygote
  fd = AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(5, 1))
  function ℳ₁((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = RigidBoundary)
    pm = PekerisModeSolver(env)
    transmission_loss(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)))
  end
  function ℳ₂((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = PressureReleaseBoundary)
    pm = PekerisModeSolver(env)
    transmission_loss(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)))
  end
  function ℳ₃((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = Rock)
    pm = PekerisModeSolver(env)
    transmission_loss(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)))
  end
  x = [20.0, 100.0, 5.0, 10.0, 500.0, 1500.0]
  @test gradient(ℳ₁, AutoForwardDiff(), x) ≈ gradient(ℳ₁, fd, x)
  @test gradient(ℳ₂, AutoForwardDiff(), x) ≈ gradient(ℳ₂, fd, x)
  @test gradient(ℳ₃, AutoForwardDiff(), x) ≈ gradient(ℳ₃, fd, x)
  @test gradient(ℳ₁, AutoZygote(), x) ≈ gradient(ℳ₁, fd, x)
  @test gradient(ℳ₂, AutoZygote(), x) ≈ gradient(ℳ₂, fd, x)
  # FIXME broken (originally due to https://github.com/SciML/NonlinearSolve.jl/issues/581)
  # @test gradient(ℳ₃, AutoZygote(), x) ≈ gradient(ℳ₃, fd, x)
  x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
  @test gradient(ℳ₁, AutoForwardDiff(), x) ≈ gradient(ℳ₁, fd, x)
  @test gradient(ℳ₂, AutoForwardDiff(), x) ≈ gradient(ℳ₂, fd, x)
  @test gradient(ℳ₃, AutoForwardDiff(), x) ≈ gradient(ℳ₃, fd, x)
  @test gradient(ℳ₁, AutoZygote(), x) ≈ gradient(ℳ₁, fd, x)
  @test gradient(ℳ₂, AutoZygote(), x) ≈ gradient(ℳ₂, fd, x)
  # FIXME broken (originally due to https://github.com/SciML/NonlinearSolve.jl/issues/581)
  # @test gradient(ℳ₃, AutoZygote(), x) ≈ gradient(ℳ₃, fd, x)
end

@testitem "pekeris-modes-adiabatic" begin
  env = @inferred UnderwaterEnvironment(
    bathymetry = SampledField([200, 150, 200]; x=[0, 2000, 5000], interp=:linear),
    soundspeed = 1500,
    density = 1000,
    seabed = FluidBoundary(2000, 2000)
  )
  pm = @inferred AdiabaticExt(PekerisModeSolver, env)
  @test startswith(sprint(show, pm), "AdiabaticExt")
  tx = @inferred AcousticSource(0, -20, 300)
  rxs = @inferred AcousticReceiverGrid2D(100:100:5000, -150)
  x = @inferred transmission_loss(pm, tx, rxs)
  y = [62.233, 46.392, 46.456, 51.631, 45.669, 54.690, 48.586,
    50.710, 49.661, 48.844, 52.003, 50.907, 59.314, 55.976, 55.497, 53.547,
    61.082, 55.273, 60.577, 53.672, 50.658, 51.736, 52.712, 53.007, 53.595,
    53.587, 53.936, 59.569, 60.200, 61.482, 57.863, 58.901, 56.467, 70.807,
    53.391, 57.480, 58.761, 55.793, 69.782, 57.063, 67.308, 61.774, 64.339,
    56.078, 65.519, 62.946, 67.441, 68.531, 54.778, 79.428]
  @test vec(x) ≈ y atol=0.1
  x = @inferred transmission_loss(pm, tx, rxs[end])
  @test x ≈ y[end] atol=0.1
  arr = @inferred arrivals(pm, tx, AcousticReceiver(5000, -50))
  @test length(arr) == 40
  @test all([a.v < 1500 for a ∈ arr])
  @test all([a.vₚ > 1500 for a ∈ arr])
  env = @inferred UnderwaterEnvironment(
    bathymetry = SampledField([200, 150, 200]; x=[-5000, -2000, 0], interp=:linear),
    soundspeed = 1500,
    density = 1000,
    seabed = FluidBoundary(2000, 2000)
  )
  pm = @inferred AdiabaticExt(PekerisModeSolver, env)
  rxs = @inferred AcousticReceiverGrid2D(-100:-100:-5000, -150)
  x = @inferred transmission_loss(pm, tx, rxs)
  @test vec(x) ≈ y atol=0.1
  x = @inferred transmission_loss(pm, tx, rxs[end])
  @test x ≈ y[end] atol=0.1
  arr = @inferred arrivals(pm, tx, AcousticReceiver(-5000, -50))
  @test length(arr) == 40
  @test all([a.v < 1500 for a ∈ arr])
  @test all([a.vₚ > 1500 for a ∈ arr])
  env = @inferred UnderwaterEnvironment(
    bathymetry = 200,
    soundspeed = 1500,
    density = 1000,
    seabed = FluidBoundary(2000, 2000)
  )
  pm1 = @inferred AdiabaticExt(PekerisModeSolver, env)
  pm2 = @inferred PekerisModeSolver(env)
  x1 = @inferred transmission_loss(pm1, tx, rxs)
  x2 = @inferred transmission_loss(pm2, tx, rxs)
  @test x1 == x2
  arr1 = @inferred arrivals(pm1, tx, AcousticReceiver(5000, -50))
  arr2 = @inferred arrivals(pm2, tx, AcousticReceiver(5000, -50))
  @test length(arr1) == length(arr2)
  @test [a.kᵣ for a ∈ arr1] == [a.kᵣ for a ∈ arr2]
  @test [a.v for a ∈ arr1] == [a.v for a ∈ arr2]
  @test [a.vₚ for a ∈ arr1] == [a.vₚ for a ∈ arr2]
end
