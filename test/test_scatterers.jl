using TestItems

@testitem "scatterers-construction" begin
  @test Ellipse <: AbstractShape
  e = Ellipse(500.0, -40.0, 20.0, 8.0)
  @test e isa Ellipse
  @test e == Ellipse(x=500, z=-40, a=20, b=8)
  @test e == Ellipse(500u"m", -40.0u"m", 20u"m", 8000u"mm")
  @test e == @inferred Ellipse(500.0, -40.0, 20.0, 8.0)
  @test e isa Ellipse{Float64}
  @test startswith(sprint(show, e), "Ellipse")
  e = Ellipse(500, -40, 20, 8; θ=30u"°")
  @test e.θ ≈ π/6
  @test e == Ellipse(x=500, z=-40, a=20, b=8, θ=π/6)
  @test_throws ErrorException Ellipse(0, 0, -1, 1)
  @test_throws ErrorException Ellipse(0, 0, 1, 0)
  s = Scatterer(e)
  @test s isa Scatterer
  @test Scatterer <: UnderwaterAcoustics.AbstractScatterer
  @test s.shape === e
  @test s.boundary === RigidBoundary
  @test s == Scatterer(e, RigidBoundary)
  @test Scatterer(e; boundary=PressureReleaseBoundary).boundary === PressureReleaseBoundary
  @test startswith(sprint(show, s), "Scatterer")
  @test location(s) == (x=500.0, y=0.0, z=-40.0)
  @test location(e) == (x=500.0, y=0.0, z=-40.0)
end

@testitem "scatterers-geometry" begin
  e = Ellipse(100.0, -50.0, 20.0, 8.0; θ=30u"°")
  @test ndims(e) == 2
  # points on the boundary satisfy the ellipse equation in the ellipse frame
  onellipse(e, p) = begin
    s, c = sincos(e.θ)
    px = c * (p.x - e.x) + s * (p.z - e.z)
    pz = -s * (p.x - e.x) + c * (p.z - e.z)
    abs2(px / e.a) + abs2(pz / e.b) ≈ 1
  end
  for u ∈ 0:0.05:1
    p = @inferred boundary_point(e, u)
    @test p.y == 0
    @test onellipse(e, p)
  end
  p0 = boundary_point(e, 0)
  p1 = boundary_point(e, 1)
  @test all(isapprox.(values(p0), values(p1); atol=1e-12))
  # analytic normals are unit length and outward
  for u ∈ 0:0.05:1
    n = @inferred normal(e, u)
    @test hypot(n.x, n.y, n.z) ≈ 1
    p = boundary_point(e, u)
    @test (p.x - e.x) * n.x + (p.z - e.z) * n.z > 0
  end
  # curvature of a circle is 1/r
  circ = Ellipse(0.0, 0.0, 5.0, 5.0)
  @test all(curvature(circ, u) ≈ 0.2 for u ∈ 0:0.1:1)
  # curvature at the ends of the major/minor axes
  e0 = Ellipse(0.0, 0.0, 4.0, 2.0)
  @test curvature(e0, 0) ≈ 4 / 2^2
  @test curvature(e0, 0.25) ≈ 2 / 4^2
  # AD-derived normal/curvature from boundary_point() agree with analytic ones
  pc = ParametricCurve(u -> boundary_point(e, u))
  @test ndims(pc) == 2
  for u ∈ 0:0.05:1
    na = normal(e, u)
    nd = normal(pc, u)
    @test all(isapprox.(values(na), values(nd); atol=1e-10))
    @test curvature(pc, u) ≈ curvature(e, u)
  end
  # sampled boundary is closed and on the ellipse
  pts = boundary_points(e; n=64)
  @test length(pts) == 65
  @test all(isapprox.(values(pts[1]), values(pts[end]); atol=1e-12))
  @test all(onellipse(e, p) for p ∈ pts)
end

@testitem "scatterers-queries" begin
  e = Ellipse(100.0, -50.0, 20.0, 8.0; θ=30u"°")
  pc = ParametricCurve(u -> boundary_point(e, u))
  # is_inside
  @test is_inside(e, (x=100.0, y=0.0, z=-50.0))
  @test !is_inside(e, (x=150.0, y=0.0, z=-50.0))
  # just outside / just inside the boundary (exactly on the boundary is float-ambiguous)
  let p = boundary_point(e, 0.1), n = normal(e, 0.1), δ = 1e-6
    @test !is_inside(e, (x=p.x + δ * n.x, y=0.0, z=p.z + δ * n.z))
    @test is_inside(e, (x=p.x - δ * n.x, y=0.0, z=p.z - δ * n.z))
  end
  for x ∈ 70:5:130, z ∈ -65:2.5:-35
    @test is_inside(pc, (x=Float64(x), y=0.0, z=Float64(z))) == is_inside(e, (x=Float64(x), y=0.0, z=Float64(z)))
  end
  # distance is signed, zero on the boundary, and matches offsets along the normal
  for u ∈ 0:0.1:0.9
    p = boundary_point(e, u)
    n = normal(e, u)
    @test abs(distance(e, p)) < 1e-8
    for δ ∈ (0.1, 2.0)
      pout = (x=p.x + δ * n.x, y=0.0, z=p.z + δ * n.z)
      pin = (x=p.x - δ * n.x, y=0.0, z=p.z - δ * n.z)
      @test distance(e, pout) ≈ δ
      @test distance(e, pin) ≈ -δ
      @test distance(pc, pout) ≈ δ atol=1e-6
      @test distance(pc, pin) ≈ -δ atol=1e-6
    end
  end
  # boundary_projection returns the nearest boundary point
  p = boundary_point(e, 0.15)
  n = normal(e, 0.15)
  q = (x=p.x + 3n.x, y=0.0, z=p.z + 3n.z)
  for shape ∈ (e, pc)
    bp = boundary_projection(shape, q)
    @test all(isapprox.(values(bp.point), values(p); atol=1e-6))
    @test all(isapprox.(values(bp.normal), values(n); atol=1e-6))
    @test bp.u ≈ 0.15 atol=1e-6
    @test bp.distance ≈ 3
    @test bp.curvature ≈ curvature(e, 0.15) atol=1e-6
  end
  # intersect_ray on an axis-aligned ellipse with known answers
  e0 = Ellipse(0.0, 0.0, 4.0, 2.0)
  pc0 = ParametricCurve(u -> boundary_point(e0, u))
  for shape ∈ (e0, pc0)
    hit = intersect_ray(shape, (x=-10.0, y=0.0, z=0.0), (x=1.0, y=0.0, z=0.0))
    @test hit !== nothing
    @test hit.t ≈ 6
    @test hit.u ≈ 0.5
    @test all(isapprox.(values(hit.point), (-4.0, 0.0, 0.0); atol=1e-9))
    @test all(isapprox.(values(hit.normal), (-1.0, 0.0, 0.0); atol=1e-9))
    @test hit.curvature ≈ 1
    # from inside, the first boundary crossing along the ray is returned
    hit = intersect_ray(shape, (x=0.0, y=0.0, z=0.0), (x=1.0, y=0.0, z=0.0))
    @test hit.t ≈ 4
    # misses
    @test intersect_ray(shape, (x=-10.0, y=0.0, z=3.0), (x=1.0, y=0.0, z=0.0)) === nothing
    @test intersect_ray(shape, (x=-10.0, y=0.0, z=0.0), (x=-1.0, y=0.0, z=0.0)) === nothing
  end
  # bounding box
  bb = bounding_box(e0)
  @test all(isapprox.(values(bb.min), (-4.0, 0.0, -2.0); atol=1e-9))
  @test all(isapprox.(values(bb.max), (4.0, 0.0, 2.0); atol=1e-9))
  bb = bounding_box(Ellipse(10.0, -20.0, 4.0, 2.0; θ=90u"°"))
  @test all(isapprox.(values(bb.min), (8.0, 0.0, -24.0); atol=1e-9))
  @test all(isapprox.(values(bb.max), (12.0, 0.0, -16.0); atol=1e-9))
end

@testitem "scatterers-parametric-surface" begin
  # unit sphere as a parametric surface
  sph = ParametricSurface() do u, v
    s, c = sincospi(2u)
    sv, cv = sincospi(v)
    (x=10 + sv * c, y=sv * s, z=-40 - cv)
  end
  @test ndims(sph) == 3
  p = boundary_point(sph, 0.25, 0.5)
  @test all(isapprox.(values(p), (10.0, 1.0, -40.0); atol=1e-12))
  n = normal(sph, 0.1, 0.3)
  @test hypot(n.x, n.y, n.z) ≈ 1
  # outward normal of a sphere is radial
  p = boundary_point(sph, 0.1, 0.3)
  @test all(isapprox.(values(n), (p.x - 10, p.y, p.z + 40); atol=1e-9))
  @test curvature(sph, 0.1, 0.3) === missing
  pts = boundary_points(sph; n=16)
  @test length(pts) == 17^2
  @test all(hypot(p.x - 10, p.y, p.z + 40) ≈ 1 for p ∈ pts)
  bb = bounding_box(sph)
  @test all(isapprox.(values(bb.min), (9.0, -1.0, -41.0); atol=1e-3))
  @test all(isapprox.(values(bb.max), (11.0, 1.0, -39.0); atol=1e-3))
end

@testitem "scatterers-env" begin
  s = Scatterer(Ellipse(500.0, -40.0, 20.0, 8.0))
  env = UnderwaterEnvironment(scatterers=s)
  @test env.scatterers == (s,)
  @test UnderwaterEnvironment(scatterers=[s]).scatterers == (s,)
  @test UnderwaterEnvironment(scatterers=(s,)).scatterers == (s,)
  s2 = Scatterer(Ellipse(700.0, -20.0, 5.0, 5.0), PressureReleaseBoundary)
  @test UnderwaterEnvironment(scatterers=(s, s2)).scatterers == (s, s2)
  @test UnderwaterEnvironment().scatterers == ()
  @test has_scatterers(env)
  @test !has_scatterers(UnderwaterEnvironment())
  @test occursin("scatterers", sprint(show, env))
  @test !occursin("scatterers", sprint(show, UnderwaterEnvironment()))
  @test (@inferred UnderwaterEnvironment(scatterers=(s,))) isa UnderwaterEnvironment
  @test UnderwaterAcoustics.env_type(UnderwaterEnvironment()) === Float64
  @test UnderwaterAcoustics.env_type(env) === Float64
  import ForwardDiff
  d = ForwardDiff.Dual(20.0, 1.0)
  envd = UnderwaterEnvironment(scatterers=Scatterer(Ellipse(500.0, -40.0, d, 8.0)))
  @test UnderwaterAcoustics.env_type(envd) <: ForwardDiff.Dual
  # custom shapes work too
  envp = UnderwaterEnvironment(scatterers=Scatterer(ParametricCurve(u -> (500 + 10cospi(2u), -40 + 10sinpi(2u)))))
  @test UnderwaterAcoustics.env_type(envp) === Float64
end

@testitem "scatterers-models" begin
  env = UnderwaterEnvironment(bathymetry=100.0, scatterers=Scatterer(Ellipse(500.0, -40.0, 20.0, 8.0)))
  @test_throws ErrorException PekerisRayTracer(env)
  @test_throws ErrorException PekerisModeSolver(env)
  @test_throws ErrorException AdiabaticExt(PekerisModeSolver, env)
end

@testitem "scatterers-∂" begin
  import ForwardDiff
  # differentiable geometry wrt shape parameters
  q = (x=150.0, y=0.0, z=-45.0)
  f = p -> distance(Ellipse(p[1], p[2], p[3], p[4]; θ=p[5]), q)
  p₀ = [100.0, -50.0, 20.0, 8.0, 0.5]
  g = ForwardDiff.gradient(f, p₀)
  @test all(isfinite, g)
  # compare with central finite differences
  for i ∈ eachindex(p₀)
    h = 1e-6
    p⁺ = copy(p₀); p⁺[i] += h
    p⁻ = copy(p₀); p⁻[i] -= h
    @test g[i] ≈ (f(p⁺) - f(p⁻)) / 2h atol=1e-4
  end
  # gradient of the signed distance wrt position is the unit outward direction,
  # both through the analytic path and the generic Newton fallback
  e = Ellipse(100.0, -50.0, 20.0, 8.0; θ=0.5)
  pc = ParametricCurve(u -> boundary_point(e, u))
  for shape ∈ (e, pc)
    gp = ForwardDiff.gradient(x -> distance(shape, (x=x[1], y=0.0, z=x[2])), [150.0, -45.0])
    @test hypot(gp...) ≈ 1
  end
  # boundary_point is differentiable wrt shape parameters
  gb = ForwardDiff.derivative(a -> boundary_point(Ellipse(0.0, 0.0, a, 2.0), 0.0).x, 4.0)
  @test gb ≈ 1
end
