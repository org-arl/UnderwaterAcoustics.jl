using TestItems

@testsnippet ReframeSetup begin
  # rigid transform: rotate by θ about z and translate by (x₀, y₀)
  function world(p; θ=deg2rad(37.0), x₀=120.0, y₀=-80.0)
    (x = x₀ + p.x * cos(θ) - p.y * sin(θ), y = y₀ + p.x * sin(θ) + p.y * cos(θ), z = p.z)
  end
end

@testitem "reframe-types" begin
  env = @inferred UnderwaterEnvironment(bathymetry=20.0u"m")
  pm1 = @inferred Reframe2D(PekerisRayTracer, env)
  pm2 = @inferred Reframe2D(PekerisModeSolver, env)
  @test pm1 isa UnderwaterAcoustics.AbstractRayPropagationModel
  @test pm2 isa UnderwaterAcoustics.AbstractModePropagationModel
  @test UnderwaterAcoustics.Reframe2DRay ∈ models(UnderwaterAcoustics.AbstractRayPropagationModel)
  @test UnderwaterAcoustics.Reframe2DMode ∈ models(UnderwaterAcoustics.AbstractModePropagationModel)
  @test startswith(sprint(show, pm1), "Reframe2D(")
  @test startswith(sprint(show, pm2), "Reframe2D(")
  @test UnderwaterAcoustics._model(pm1) === PekerisRayTracer
  @test UnderwaterAcoustics._model(pm2) === PekerisModeSolver
end

@testitem "reframe-xform" begin
  using UnderwaterAcoustics: RigidXform, _forward, _inverse, _isidentity, xyz
  t = RigidXform(120.0, -80.0, cosd(37.0), sind(37.0))
  for p ∈ (xyz(0, 0, -5), xyz(100, 50, -10), xyz(-30, -70, 0))
    @test collect(_inverse(t, _forward(t, p))) ≈ collect(p) atol=1e-12
    @test collect(_forward(t, _inverse(t, p))) ≈ collect(p) atol=1e-12
    @test _forward(t, p).z == p.z
  end
  @test _isidentity(RigidXform(0.0, 0.0, 1.0, 0.0))
  @test !_isidentity(RigidXform(1.0, 0.0, 1.0, 0.0))
  @test !_isidentity(t)
  # frames with non-float number types are never treated as identity (AD safety)
  @test !_isidentity(RigidXform(big(0.0), big(0.0), big(1.0), big(0.0)))
end

@testitem "reframe-passthrough" setup=[ReframeSetup] begin
  env = UnderwaterEnvironment(bathymetry=20.0u"m", seabed=SandySilt)
  pm0 = PekerisRayTracer(env)
  pm = Reframe2D(PekerisRayTracer, env)
  tx = AcousticSource((x=0.0, z=-5.0), 1000.0)
  rx = AcousticReceiver((x=100.0, z=-10.0))
  @test acoustic_field(pm, tx, rx) == acoustic_field(pm0, tx, rx)
  @test transmission_loss(pm, tx, rx) == transmission_loss(pm0, tx, rx)
  arr = arrivals(pm, tx, rx)
  arr0 = arrivals(pm0, tx, rx)
  @test [a.t for a ∈ arr] == [a.t for a ∈ arr0]
  @test [a.path for a ∈ arr] == [a.path for a ∈ arr0]
  @test samples(impulse_response(pm, tx, rx, 10000.0)) == samples(impulse_response(pm0, tx, rx, 10000.0))
  # canonical scenario yields identity transform and untouched environment
  t = UnderwaterAcoustics._derive_xform(pm, tx, rx)
  @test UnderwaterAcoustics._isidentity(t)
  @test UnderwaterAcoustics._reframe(t, env) === env
end

@testitem "reframe-rotation-invariance" setup=[ReframeSetup] begin
  # PekerisRayTracer is natively frame-invariant with true 3D ray paths, so it
  # provides a gold standard to validate the reframing and back-transformation
  env = UnderwaterEnvironment(bathymetry=20.0u"m", seabed=SandySilt)
  pm0 = PekerisRayTracer(env)
  pm = Reframe2D(PekerisRayTracer, env)
  tx0 = AcousticSource((x=0.0, z=-5.0), 1000.0)
  rx0 = AcousticReceiver((x=800.0, z=-10.0))
  tx = AcousticSource(world((x=0.0, y=0.0, z=-5.0)), 1000.0)
  rx = AcousticReceiver(world((x=800.0, y=0.0, z=-10.0)))
  @test acoustic_field(pm, tx, rx) ≈ acoustic_field(pm0, tx0, rx0) rtol=1e-10
  @test acoustic_field(pm, tx, rx) ≈ acoustic_field(pm0, tx, rx) rtol=1e-10
  arr = arrivals(pm, tx, rx)
  arr0 = arrivals(pm0, tx0, rx0)
  @test [a.t for a ∈ arr] ≈ [a.t for a ∈ arr0] rtol=1e-12
  @test [a.θₛ for a ∈ arr] ≈ [a.θₛ for a ∈ arr0] rtol=1e-12
  for (a, a0) ∈ zip(arr, arr0)
    # eigenray paths must be returned in world coordinates
    @test all(collect(p) ≈ collect(world(p0)) for (p, p0) ∈ zip(a.path, a0.path))
  end
  @test collect(arr[1].path[1]) ≈ collect(location(tx)) atol=1e-9
  @test collect(arr[1].path[end]) ≈ collect(location(rx)) atol=1e-9
end

@testitem "reframe-modes-offaxis" setup=[ReframeSetup] begin
  env = UnderwaterEnvironment(bathymetry=20.0u"m", soundspeed=1500.0, seabed=Rock)
  pm0 = PekerisModeSolver(env)
  pm = Reframe2D(PekerisModeSolver, env)
  tx0 = AcousticSource((x=0.0, z=-5.0), 1000.0)
  rx0 = AcousticReceiver((x=500.0, z=-10.0))
  x0 = acoustic_field(pm0, tx0, rx0)
  # off-axis scenario, same range
  tx = AcousticSource(world((x=0.0, y=0.0, z=-5.0)), 1000.0)
  rx = AcousticReceiver(world((x=500.0, y=0.0, z=-10.0)))
  @test acoustic_field(pm, tx, rx) ≈ x0 rtol=1e-10
  # the mode solver itself now computes horizontal range from both x and y
  @test acoustic_field(pm0, tx, rx) ≈ x0 rtol=1e-10
  @test acoustic_field(pm0, tx0, AcousticReceiver((x=300.0, y=400.0, z=-10.0))) ≈ x0 rtol=1e-10
end

@testitem "reframe-validation" setup=[ReframeSetup] begin
  env = UnderwaterEnvironment(bathymetry=20.0u"m", soundspeed=1500.0, seabed=Rock)
  pm = Reframe2D(PekerisModeSolver, env)
  tx = AcousticSource((x=100.0, y=50.0, z=-5.0), 1000.0)
  # single receiver anywhere is always a valid 2D scenario
  @test acoustic_field(pm, tx, AcousticReceiver((x=700.0, y=-250.0, z=-10.0))) isa Complex
  # non-colinear receivers are rejected
  rxs = [AcousticReceiver((x=500.0, y=50.0, z=-10.0)), AcousticReceiver((x=500.0, y=300.0, z=-10.0))]
  @test_throws r"not 2D" acoustic_field(pm, tx, rxs)
  # colinearity tolerance is honored
  rxs = [AcousticReceiver((x=500.0, y=50.0, z=-10.0)), AcousticReceiver((x=300.0, y=50.05, z=-10.0))]
  @test acoustic_field(pm, tx, rxs) isa AbstractVector
  pm1 = Reframe2D(PekerisModeSolver, env; atol=1e-3)
  @test_throws r"not 2D" acoustic_field(pm1, tx, rxs)
  # receiver directly below the source
  @test acoustic_field(pm, tx, AcousticReceiver((x=100.0, y=50.0, z=-10.0))) isa Complex
  # receivers on both sides of the source (supported by the mode solver)
  rxs = [AcousticReceiver((x=600.0, y=50.0, z=-10.0)), AcousticReceiver((x=-400.0, y=50.0, z=-10.0))]
  x = acoustic_field(pm, tx, rxs)
  @test x[1] ≈ x[2] rtol=1e-10
  # source location is required
  @test_throws r"known location" acoustic_field(pm, AcousticSource(nothing, 1000.0), AcousticReceiver((x=500.0, z=-10.0)))
end

@testitem "reframe-grid" setup=[ReframeSetup] begin
  env = UnderwaterEnvironment(bathymetry=20.0u"m", soundspeed=1500.0, seabed=Rock)
  pm0 = PekerisModeSolver(env)
  pm = Reframe2D(PekerisModeSolver, env)
  # source offset along x: grid remains a grid (single model run)
  tx = AcousticSource((x=500.0, z=-5.0), 1000.0)
  grid = AcousticReceiverGrid2D(1000.0:10.0:1100.0, -15.0:5.0:-5.0)
  grid0 = AcousticReceiverGrid2D(500.0:10.0:600.0, -15.0:5.0:-5.0)
  tx0 = AcousticSource((x=0.0, z=-5.0), 1000.0)
  x = transmission_loss(pm, tx, grid)
  @test size(x) == size(grid)
  @test x ≈ transmission_loss(pm0, tx0, grid0) rtol=1e-10
  t = UnderwaterAcoustics._derive_xform(pm, tx, grid)
  @test UnderwaterAcoustics._reframe(t, grid) isa AcousticReceiverGrid2D
  # grid entirely behind the source (180° rotation)
  gridb = AcousticReceiverGrid2D(-1000.0:-10.0:-1100.0, -15.0:5.0:-5.0)
  xb = transmission_loss(pm, tx0, gridb)
  @test xb ≈ transmission_loss(pm0, tx0, AcousticReceiverGrid2D(1000.0:10.0:1100.0, -15.0:5.0:-5.0)) rtol=1e-10
  # single vertical column of receivers works with a source anywhere
  col = AcousticReceiverGrid2D(1000.0, -15.0:5.0:-5.0)
  txo = AcousticSource((x=400.0, y=-800.0, z=-5.0), 1000.0)
  r = hypot(1000.0 - 400.0, 800.0)
  xc = transmission_loss(pm, txo, col)
  @test vec(xc) ≈ vec(transmission_loss(pm0, tx0, AcousticReceiverGrid2D(r, -15.0:5.0:-5.0))) rtol=1e-10
  # source off the plane of a proper grid is rejected
  txy = AcousticSource((x=0.0, y=100.0, z=-5.0), 1000.0)
  @test_throws r"not 2D" transmission_loss(pm, txy, grid)
end

@testitem "reframe-env-fields" begin
  using UnderwaterAcoustics: RigidXform, _reframe_field, _reframe, TransformedField
  struct _TestField <: UnderwaterAcoustics.PositionDependent end
  Base.minimum(::_TestField) = 0.0
  Base.maximum(::_TestField) = 1000.0
  (::_TestField)(pos::UnderwaterAcoustics.XYZ) = 100.0 + 0.01 * pos.x - 0.02 * pos.y
  θ = deg2rad(30.0)
  t = RigidXform(100.0, -50.0, cos(θ), sin(θ))
  wx(x′) = 100.0 + x′ * cos(θ)                    # world x on track at range x′
  wxy(x′) = (100.0 + x′ * cos(θ), -50.0 + x′ * sin(θ))
  # constants and depth-dependent fields are untouched
  @test _reframe_field(t, 1500.0) === 1500.0
  @test _reframe_field(t, RigidBoundary) === RigidBoundary
  qz = SampledField([1500.0, 1520.0]; z=[-100.0, 0.0])
  @test _reframe_field(t, qz) === qz
  # x-varying fields are remapped exactly, preserving type
  q = SampledField([200.0, 150.0, 210.0]; x=[0.0, 2000.0, 5000.0])
  q2 = _reframe_field(t, q)
  @test q2 isa UnderwaterAcoustics.SampledFieldX
  @test all(q2((x=x′, y=0.0, z=0.0)) ≈ q((x=wx(x′), y=0.0, z=0.0)) for x′ ∈ -3000.0:97.0:8000.0)
  # 180° - θ rotation (cosθ < 0) keeps the sample grid ascending
  t2 = RigidXform(100.0, -50.0, cos(deg2rad(150.0)), sin(deg2rad(150.0)))
  q3 = _reframe_field(t2, q)
  @test issorted(q3.xrange)
  @test all(q3((x=x′, y=0.0, z=0.0)) ≈ q((x=100.0 + x′ * cos(deg2rad(150.0)), y=0.0, z=0.0)) for x′ ∈ 0.0:50.0:5000.0)
  # track perpendicular to the direction of variation: field becomes constant
  t3 = RigidXform(1000.0, 0.0, 0.0, 1.0)
  @test _reframe_field(t3, q) ≈ q((x=1000.0, y=0.0, z=0.0))
  # xz-varying fields
  qxz = SampledField([1500.0 1520.0; 1490.0 1500.0; 1480.0 1490.0]; x=[0.0, 500.0, 1000.0], z=[-100.0, 0.0])
  qxz2 = _reframe_field(t, qxz)
  @test qxz2 isa UnderwaterAcoustics.SampledFieldXZ
  @test all(qxz2((x=x′, y=0.0, z=z)) ≈ qxz((x=wx(x′), y=0.0, z=z)) for x′ ∈ -500.0:37.0:2000.0, z ∈ -90.0:13.0:0.0)
  qxz3 = _reframe_field(t3, qxz)
  @test qxz3 isa UnderwaterAcoustics.SampledFieldZ
  @test all(qxz3((x=0.0, y=0.0, z=z)) ≈ qxz((x=1000.0, y=0.0, z=z)) for z ∈ -100.0:10.0:0.0)
  # xy-varying fields are sliced along the track
  qxy = SampledField([100.0 110.0; 120.0 130.0]; x=[0.0, 1000.0], y=[-500.0, 500.0])
  qxy2 = _reframe_field(t, qxy)
  @test qxy2 isa UnderwaterAcoustics.SampledFieldX
  @test all(qxy2((x=x′, y=0.0, z=0.0)) ≈ qxy((x=wxy(x′)[1], y=wxy(x′)[2], z=0.0)) for x′ ∈ 0.0:23.0:800.0)
  # xyz-varying fields
  vals = reshape(1.0:24.0, 3, 4, 2)
  qxyz = SampledField(vals; x=[0.0, 500.0, 1000.0], y=[-800.0, -300.0, 200.0, 700.0], z=[-50.0, 0.0])
  qxyz2 = _reframe_field(t, qxyz)
  @test qxyz2 isa UnderwaterAcoustics.SampledFieldXZ
  @test all(qxyz2((x=x′, y=0.0, z=z)) ≈ qxyz((x=wxy(x′)[1], y=wxy(x′)[2], z=z)) for x′ ∈ 0.0:31.0:600.0, z ∈ -45.0:11.0:0.0)
  # arbitrary position-dependent callables are wrapped lazily
  f = _TestField()
  f2 = _reframe_field(t, f)
  @test f2 isa TransformedField
  @test startswith(sprint(show, f2), "TransformedField")
  @test is_range_dependent(f2)
  @test !is_constant(f2)
  @test minimum(f2) == 0.0 && maximum(f2) == 1000.0
  @test all(value(f2, (x′, 0.0, 0.0)) ≈ f((x=wxy(x′)[1], y=wxy(x′)[2], z=0.0)) for x′ ∈ 0.0:100.0:2000.0)
  # environment reframing preserves untouched environments
  env = UnderwaterEnvironment(bathymetry=20.0u"m", soundspeed=1500.0)
  @test _reframe(t, env) === env
  env2 = UnderwaterEnvironment(bathymetry=q, soundspeed=1500.0)
  env2t = _reframe(t, env2)
  @test env2t !== env2
  @test env2t.soundspeed == env2.soundspeed
  @test value(env2t.bathymetry, (500.0, 0.0, 0.0)) ≈ value(env2.bathymetry, (wx(500.0), 0.0, 0.0))
end

@testitem "reframe-scatterers" setup=[ReframeSetup] begin
  using UnderwaterAcoustics: RigidXform, _reframe, _forward
  # translation along x: Ellipse remaps analytically
  t = RigidXform(300.0, 0.0, 1.0, 0.0)
  env = UnderwaterEnvironment(bathymetry=50.0, soundspeed=1500.0,
    scatterers=Scatterer(Ellipse(500.0, -20.0, 10.0, 4.0; θ=0.3)))
  env1 = _reframe(t, env)
  e1 = env1.scatterers[1].shape
  @test e1 isa Ellipse
  @test e1.x == 200.0 && e1.z == -20.0 && e1.a == 10.0 && e1.b == 4.0 && e1.θ == 0.3
  @test env1.scatterers[1].boundary === RigidBoundary
  # 180° rotation: mirrored ellipse, boundary traced consistently
  t2 = RigidXform(300.0, 0.0, -1.0, 0.0)
  e2 = _reframe(t2, env).scatterers[1].shape
  @test e2.x == -200.0 && e2.θ == -0.3
  @test is_inside(e2, (x=-200.0, y=0.0, z=-20.0))
  bb = bounding_box(e2)
  @test bb.min.x ≈ -(bounding_box(Ellipse(200.0, -20.0, 10.0, 4.0; θ=0.3)).max.x)
  # generic 2-D shape: boundary points map through the transform, normals stay outward
  circ = ParametricCurve(u -> (500.0 + 10 * cospi(2u), -20.0 + 10 * sinpi(2u)))
  envc = UnderwaterEnvironment(bathymetry=50.0, soundspeed=1500.0, scatterers=Scatterer(circ))
  for tt ∈ (t, t2)
    c1 = _reframe(tt, envc).scatterers[1].shape
    @test ndims(c1) == 2
    c = _forward(tt, (x=500.0, y=0.0, z=-20.0))     # transformed center
    for u ∈ 0.0:0.1:0.9
      p = boundary_point(c1, u)
      @test hypot(p.x - c.x, p.z - c.z) ≈ 10.0 atol=1e-9
      n = normal(c1, u)
      @test n.x * (p.x - c.x) + n.z * (p.z - c.z) > 0   # outward
    end
    @test is_inside(c1, c)
  end
  # 3-D shape: transforms rigidly under any rotation
  sph = ParametricSurface((u, v) -> (500.0 + 10 * cospi(2u) * sinpi(v),
    200.0 + 10 * sinpi(2u) * sinpi(v), -20.0 - 10 * cospi(v)))
  t3 = RigidXform(120.0, -80.0, cosd(37.0), sind(37.0))
  s1 = UnderwaterAcoustics._reframe_shape(t3, sph)
  @test ndims(s1) == 3
  c = _forward(t3, (x=500.0, y=200.0, z=-20.0))
  for u ∈ 0.05:0.2:0.85, v ∈ 0.1:0.2:0.9
    p = boundary_point(s1, u, v)
    @test hypot(p.x - c.x, p.y - c.y, p.z - c.z) ≈ 10.0 atol=1e-9
    n = normal(s1, u, v)
    @test n.x * (p.x - c.x) + n.y * (p.y - c.y) + n.z * (p.z - c.z) > 0
  end
  # 2-D scatterers + rotated propagation plane is not a 2D scenario
  pm = Reframe2D(PekerisRayTracer, env)
  tx = AcousticSource((x=0.0, z=-5.0), 1000.0)
  @test_throws r"not 2D" acoustic_field(pm, tx, AcousticReceiver((x=700.0, y=400.0, z=-10.0)))
  # in-plane scenarios reach the wrapped model, which decides on scatterer support
  @test_throws r"scatterers" acoustic_field(pm, tx, AcousticReceiver((x=700.0, y=0.0, z=-10.0)))
end

@testitem "reframe-∂" setup=[ReframeSetup] begin
  using DifferentiationInterface
  import ForwardDiff, FiniteDifferences, Zygote
  fd = AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(5, 1))
  function ℳ((txx, txy, txz, rxx, rxy, rxz))
    env = UnderwaterEnvironment(bathymetry=20.0, soundspeed=1500.0, seabed=SandySilt)
    pm = Reframe2D(PekerisRayTracer, env)
    transmission_loss(pm, AcousticSource((x=txx, y=txy, z=txz), 1000.0), AcousticReceiver((x=rxx, y=rxy, z=rxz)))
  end
  # rotated scenario
  x = [values(world((x=0.0, y=0.0, z=-5.0)))..., values(world((x=800.0, y=0.0, z=-10.0)))...]
  @test gradient(ℳ, AutoForwardDiff(), x) ≈ gradient(ℳ, fd, x) rtol=1e-4
  @test gradient(ℳ, AutoZygote(), x) ≈ gradient(ℳ, fd, x) rtol=1e-4
  # canonical scenario with y exactly 0: partials w.r.t. y must still be computed
  xc = [0.0, 0.0, -5.0, 800.0, 0.0, -10.0]
  @test gradient(ℳ, AutoForwardDiff(), xc) ≈ gradient(ℳ, fd, xc) rtol=1e-4
  @test gradient(ℳ, AutoZygote(), xc) ≈ gradient(ℳ, fd, xc) rtol=1e-4
end

@testitem "reframe-ir+channel" setup=[ReframeSetup] begin
  env = UnderwaterEnvironment(bathymetry=20.0u"m", seabed=SandySilt)
  pm0 = PekerisRayTracer(env)
  pm = Reframe2D(PekerisRayTracer, env)
  tx0 = AcousticSource((x=0.0, z=-5.0), 1000.0)
  rx0 = AcousticReceiver((x=800.0, z=-10.0))
  tx = AcousticSource(world((x=0.0, y=0.0, z=-5.0)), 1000.0)
  rx = AcousticReceiver(world((x=800.0, y=0.0, z=-10.0)))
  ir0 = impulse_response(pm0, tx0, rx0, 10000.0)
  ir = impulse_response(pm, tx, rx, 10000.0)
  @test samples(ir) ≈ samples(ir0) rtol=1e-10
  ch0 = channel(pm0, tx0, rx0, 10000.0)
  ch = channel(pm, tx, rx, 10000.0)
  x = randn(500)
  @test samples(transmit(ch, x)) ≈ samples(transmit(ch0, x)) rtol=1e-8
  # mode models
  envm = UnderwaterEnvironment(bathymetry=20.0u"m", soundspeed=1500.0, seabed=Rock)
  pmm0 = PekerisModeSolver(envm)
  pmm = Reframe2D(PekerisModeSolver, envm)
  irm0 = impulse_response(pmm0, tx0, rx0, 8000.0)
  irm = impulse_response(pmm, tx, rx, 8000.0)
  @test samples(irm) ≈ samples(irm0) rtol=1e-8
end
