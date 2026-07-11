import ForwardDiff

export AbstractShape, Ellipse, ParametricCurve, ParametricSurface, Scatterer
export boundary_point, normal, curvature, distance, boundary_projection
export intersect_ray, is_inside, boundary_points, bounding_box, has_scatterers

###############################################################################
### shapes API

"""
Superclass for all scatterer shapes.

A shape describes the geometry of a closed object through a parametric map of
its boundary over a unit parameter domain. The boundary of a 2-D shape is a
curve parametrized by `u ∈ [0,1]`, and the boundary of a 3-D shape is a surface
parametrized by `(u, v) ∈ [0,1]²`. A concrete shape must define:

- `Base.ndims(shape)` — spatial dimension of the shape (2 or 3)
- `boundary_point(shape, u...)` — parametric boundary position map

Optionally, a shape may provide analytic `normal()` and `curvature()` methods;
if absent, they are automatically derived from `boundary_point()` using
automatic differentiation. A shape may also override the derived queries
`distance()`, `boundary_projection()`, `intersect_ray()`, `is_inside()`,
`boundary_points()`, `bounding_box()` and `location()` with analytic
implementations for speed and accuracy.
"""
abstract type AbstractShape end

"""
    boundary_point(shape, u)
    boundary_point(shape, u, v)

Get the position on the boundary of `shape` at the given boundary parameters.
The boundary of a 2-D shape is parametrized by a single parameter `u ∈ [0,1]`,
and the boundary of a 3-D shape by two parameters `(u, v) ∈ [0,1]²`. The
parametric map must trace the full closed boundary as the parameters sweep the
unit domain, with `u = 0` and `u = 1` mapping to the same point. The position
is returned as an `XYZ` named tuple with coordinates in meters.

2-D shapes lie in the x–z plane (`y = 0`), with `z` negative downward. The
boundary of a 2-D shape must be traversed counter-clockwise (with `x` rightward
and `z` upward) as `u` increases, so that automatically derived normals point
outward. Similarly, the parametrization of a 3-D shape must be right-handed
(`∂p/∂u × ∂p/∂v` pointing outward).
"""
function boundary_point end

"""
    normal(shape, u)
    normal(shape, u, v)

Get the unit outward normal to the boundary of `shape` at the given boundary
parameters, as an `XYZ` named tuple. If a shape does not provide an analytic
normal, it is derived from `boundary_point()` using automatic differentiation.
"""
function normal(shape::AbstractShape, u)
  dx, dz = _dbdu(shape, u)
  m = hypot(dz, dx)
  xyz(dz / m, -dx / m)
end

function normal(shape::AbstractShape, u, v)
  pu = _dbduv(shape, u, v, 1)
  pv = _dbduv(shape, u, v, 2)
  nx = pu[2] * pv[3] - pu[3] * pv[2]
  ny = pu[3] * pv[1] - pu[1] * pv[3]
  nz = pu[1] * pv[2] - pu[2] * pv[1]
  m = hypot(nx, ny, nz)
  xyz(nx / m, ny / m, nz / m)
end

"""
    curvature(shape, u)
    curvature(shape, u, v)

Get the signed curvature (in 1/m) of the boundary of a 2-D `shape` at boundary
parameter `u`. The curvature is positive where the boundary curves toward the
interior of the shape (convex), and negative where it curves away (concave).
If a shape does not provide an analytic curvature, it is derived from
`boundary_point()` using automatic differentiation.

The return convention for the curvature of 3-D shapes is not yet settled, and
`missing` is returned for 3-D shapes that do not provide a curvature.
"""
function curvature(shape::AbstractShape, u)
  fx = u -> boundary_point(shape, u).x
  fz = u -> boundary_point(shape, u).z
  x′ = ForwardDiff.derivative(fx, u)
  z′ = ForwardDiff.derivative(fz, u)
  x″ = ForwardDiff.derivative(u -> ForwardDiff.derivative(fx, u), u)
  z″ = ForwardDiff.derivative(u -> ForwardDiff.derivative(fz, u), u)
  (x′ * z″ - z′ * x″) / hypot(x′, z′)^3
end

curvature(shape::AbstractShape, u, v) = missing

"""
    distance(shape, pos)

Get the signed distance (in meters) from position `pos` to the boundary of
`shape`. The distance is negative if `pos` is inside the shape, positive if
outside, and zero on the boundary. The magnitude is the Euclidean distance to
the nearest boundary point. A shape may return a conservative lower bound on
the magnitude far away from the boundary, but the distance must be exact (with
the correct sign and zero level set) near the boundary.

For 2-D shapes, the distance is computed in the x–z plane, ignoring the `y`
coordinate of `pos`.
"""
distance(shape::AbstractShape, pos::XYZ) = boundary_projection(shape, pos).distance
distance(shape::AbstractShape, pos) = distance(shape, xyz(pos))

"""
    boundary_projection(shape, pos)

Project position `pos` onto the boundary of `shape`, i.e., find the boundary
point nearest to `pos`. Returns a named tuple with fields:

- `u`: boundary parameter(s) of the nearest boundary point
- `point`: position of the nearest boundary point
- `normal`: unit outward normal at that point
- `curvature`: curvature at that point
- `distance`: signed distance from `pos` to the boundary (negative inside)

For 2-D shapes, the projection is computed in the x–z plane, ignoring the `y`
coordinate of `pos`.
"""
function boundary_projection(shape::AbstractShape, pos::XYZ)
  ndims(shape) == 2 || error("Generic boundary_projection is currently only available for 2-D shapes")
  n = 64
  # coarse global search over the parameter domain, followed by Newton refinement
  u = 0 / n
  d² = _dist²(shape, u, pos)
  for i ∈ 1:n-1
    d²ᵢ = _dist²(shape, i / n, pos)
    d²ᵢ < d² && ((u, d²) = (i / n, d²ᵢ))
  end
  u = _newton_projection(shape, float(u), pos, 1 / 2n)
  _projection_bundle(shape, u, pos)
end

boundary_projection(shape::AbstractShape, pos) = boundary_projection(shape, xyz(pos))

"""
    intersect_ray(shape, origin, dir)

Find the nearest intersection of the ray `origin + t * dir` (for `t ≥ 0`) with
the boundary of `shape`. Returns `nothing` if the ray does not intersect the
boundary, and otherwise a named tuple with fields:

- `t`: ray parameter at the intersection (distance along the ray in units of
  the length of `dir`)
- `u`: boundary parameter(s) of the intersection point
- `point`: position of the intersection point
- `normal`: unit outward normal at the intersection point
- `curvature`: curvature at the intersection point

For 2-D shapes, the intersection is computed in the x–z plane, ignoring the
`y` coordinates of `origin` and `dir`.
"""
function intersect_ray(shape::AbstractShape, origin::XYZ, dir::XYZ)
  ndims(shape) == 2 || error("Generic intersect_ray is currently only available for 2-D shapes")
  n = 256
  # seed candidates from a sampled polygon, then refine on the parametric system
  pts = boundary_points(shape; n)
  best = nothing
  for i ∈ 1:n
    p1 = pts[i]
    p2 = pts[i+1]
    # solve origin + t dir = p1 + s (p2 - p1) for (t, s)
    ex = p2.x - p1.x
    ez = p2.z - p1.z
    det = dir.z * ex - dir.x * ez
    det == 0 && continue
    bx = p1.x - origin.x
    bz = p1.z - origin.z
    t = (bz * ex - bx * ez) / det
    s = (dir.x * bz - dir.z * bx) / det
    (-0.01 ≤ s ≤ 1.01 && t ≥ -1e-6) || continue
    u, t = _newton_ray(shape, float((i - 1 + s) / n), float(t), origin, dir)
    t ≥ 0 || continue
    (best === nothing || t < best[2]) && (best = (u, t))
  end
  best === nothing && return nothing
  u, t = best
  u = mod(u, 1)
  (t=t, u=u, point=boundary_point(shape, u), normal=normal(shape, u), curvature=curvature(shape, u))
end

intersect_ray(shape::AbstractShape, origin, dir) = intersect_ray(shape, xyz(origin), xyz(dir))

"""
    is_inside(shape, pos)

Return `true` if position `pos` is strictly inside `shape`, and `false`
otherwise. For 2-D shapes, the test is performed in the x–z plane, ignoring
the `y` coordinate of `pos`.
"""
is_inside(shape::AbstractShape, pos::XYZ) = distance(shape, pos) < 0
is_inside(shape::AbstractShape, pos) = is_inside(shape, xyz(pos))

"""
    boundary_points(shape; n=128)

Sample points on the boundary of `shape`, e.g. for plotting. For 2-D shapes,
`n+1` points tracing the closed boundary are returned (with the first and last
point coinciding). For 3-D shapes, points sampled on an `(n+1) × (n+1)` grid
over the parameter domain are returned.
"""
function boundary_points(shape::AbstractShape; n=128)
  if ndims(shape) == 2
    [boundary_point(shape, u) for u ∈ range(0, 1; length=n+1)]
  else
    vec([boundary_point(shape, u, v) for u ∈ range(0, 1; length=n+1), v ∈ range(0, 1; length=n+1)])
  end
end

"""
    bounding_box(shape)

Get an axis-aligned bounding box of `shape`, as a named tuple `(min, max)` of
positions. The generic implementation estimates the box from sampled boundary
points, and may slightly underestimate the extent of the shape.
"""
function bounding_box(shape::AbstractShape)
  pts = boundary_points(shape; n=256)
  (min=xyz(minimum(p.x for p ∈ pts), minimum(p.y for p ∈ pts), minimum(p.z for p ∈ pts)),
   max=xyz(maximum(p.x for p ∈ pts), maximum(p.y for p ∈ pts), maximum(p.z for p ∈ pts)))
end

"""
    location(shape::AbstractShape)
    location(s::Scatterer)

Get the location of the center of a shape or scatterer.
"""
function location(shape::AbstractShape)
  bb = bounding_box(shape)
  xyz((bb.min.x + bb.max.x) / 2, (bb.min.y + bb.max.y) / 2, (bb.min.z + bb.max.z) / 2)
end

### AD helpers for shapes that only define boundary_point()

# ∂/∂u of the x and z coordinates of the boundary of a 2-D shape
function _dbdu(shape::AbstractShape, u)
  dx = ForwardDiff.derivative(u -> boundary_point(shape, u).x, u)
  dz = ForwardDiff.derivative(u -> boundary_point(shape, u).z, u)
  (dx, dz)
end

# partial derivative of the boundary of a 3-D shape wrt parameter i
function _dbduv(shape::AbstractShape, u, v, i)
  f = i == 1 ? (u -> boundary_point(shape, u, v)) : (v -> boundary_point(shape, u, v))
  w = i == 1 ? u : v
  (ForwardDiff.derivative(w -> f(w).x, w),
   ForwardDiff.derivative(w -> f(w).y, w),
   ForwardDiff.derivative(w -> f(w).z, w))
end

### numerical machinery for the generic derived queries (2-D)

# squared distance from pos to the boundary point at parameter u, in the x-z plane
function _dist²(shape::AbstractShape, u, pos::XYZ)
  p = boundary_point(shape, u)
  abs2(p.x - pos.x) + abs2(p.z - pos.z)
end

# Newton refinement of the nearest-point parameter, seeded near the global minimum;
# converges to a stationary point so that AD duals propagate correctly (envelope theorem)
function _newton_projection(shape::AbstractShape, u, pos::XYZ, maxstep)
  for _ ∈ 1:32
    p = boundary_point(shape, u)
    x′, z′ = _dbdu(shape, u)
    x″ = ForwardDiff.derivative(u -> _dbdu(shape, u)[1], u)
    z″ = ForwardDiff.derivative(u -> _dbdu(shape, u)[2], u)
    g = (p.x - pos.x) * x′ + (p.z - pos.z) * z′
    g′ = x′^2 + z′^2 + (p.x - pos.x) * x″ + (p.z - pos.z) * z″
    g′ == 0 && break
    Δ = clamp(g / g′, -maxstep, maxstep)
    u -= Δ
    abs(Δ) < 1e-13 && break
  end
  u
end

# assemble the boundary_projection() result bundle at parameter u
function _projection_bundle(shape::AbstractShape, u, pos::XYZ)
  u = mod(u, 1)
  p = boundary_point(shape, u)
  n̂ = normal(shape, u)
  κ = curvature(shape, u)
  d = hypot(p.x - pos.x, p.z - pos.z)
  (pos.x - p.x) * n̂.x + (pos.z - p.z) * n̂.z < 0 && (d = -d)
  (u=u, point=p, normal=n̂, curvature=κ, distance=d)
end

# Newton refinement of a ray-boundary intersection on the parametric system
# boundary_point(u) == origin + t * dir, in the x-z plane
function _newton_ray(shape::AbstractShape, u, t, origin::XYZ, dir::XYZ)
  for _ ∈ 1:32
    p = boundary_point(shape, u)
    x′, z′ = _dbdu(shape, u)
    fx = p.x - origin.x - t * dir.x
    fz = p.z - origin.z - t * dir.z
    det = -x′ * dir.z + z′ * dir.x
    det == 0 && break
    Δu = (-fx * dir.z + fz * dir.x) / det
    Δt = (x′ * fz - z′ * fx) / det
    u -= Δu
    t -= Δt
    abs(Δu) < 1e-13 && abs(Δt) < 1e-13 && break
  end
  (u, t)
end

###############################################################################
### shapes

"""
    Ellipse(x, z, a, b; θ=0)
    Ellipse(; x, z, a, b, θ=0)

Create a 2-D elliptical shape in the x–z plane, centered at `(x, z)` (in
meters), with semi-major axis `a` and semi-minor axis `b` (in meters), and
with the major axis rotated by angle `θ` (in radians) counter-clockwise from
the x-axis.

# Examples
```julia-repl
julia> Ellipse(500.0, -40.0, 20.0, 8.0)
Ellipse(x=500.0, z=-40.0, a=20.0, b=8.0)

julia> Ellipse(x=500, z=-40, a=20, b=8, θ=30u"°")
Ellipse(x=500.0, z=-40.0, a=20.0, b=8.0, θ=0.5235987755982988)
```
"""
struct Ellipse{T} <: AbstractShape
  x::T
  z::T
  a::T
  b::T
  θ::T
  function Ellipse(x, z, a, b, θ)
    x = in_units(u"m", x)
    z = in_units(u"m", z)
    a = in_units(u"m", a)
    b = in_units(u"m", b)
    θ = in_units(u"rad", θ)
    a > 0 || error("Semi-major axis must be positive")
    b > 0 || error("Semi-minor axis must be positive")
    x, z, a, b, θ = float.(promote(x, z, a, b, θ))
    new{typeof(x)}(x, z, a, b, θ)
  end
end

Ellipse(x, z, a, b; θ=0) = Ellipse(x, z, a, b, θ)
Ellipse(; x, z, a, b, θ=0) = Ellipse(x, z, a, b, θ)

function Base.show(io::IO, e::Ellipse)
  print(io, "Ellipse(x=$(e.x), z=$(e.z), a=$(e.a), b=$(e.b)")
  e.θ == 0 || print(io, ", θ=$(e.θ)")
  print(io, ")")
end

Base.ndims(::Ellipse) = 2

# world position of a point given in the ellipse frame
function _eworld(e::Ellipse, px, pz)
  sθ, cθ = sincos(e.θ)
  xyz(e.x + cθ * px - sθ * pz, 0, e.z + sθ * px + cθ * pz)
end

# ellipse frame coordinates of a world position
function _eframe(e::Ellipse, pos::XYZ)
  sθ, cθ = sincos(e.θ)
  dx = pos.x - e.x
  dz = pos.z - e.z
  (cθ * dx + sθ * dz, -sθ * dx + cθ * dz)
end

function boundary_point(e::Ellipse, u)
  s, c = sincospi(2u)
  _eworld(e, e.a * c, e.b * s)
end

function normal(e::Ellipse, u)
  s, c = sincospi(2u)
  nx = e.b * c
  nz = e.a * s
  m = hypot(nx, nz)
  sθ, cθ = sincos(e.θ)
  xyz((cθ * nx - sθ * nz) / m, 0, (sθ * nx + cθ * nz) / m)
end

function curvature(e::Ellipse, u)
  s, c = sincospi(2u)
  e.a * e.b / hypot(e.a * s, e.b * c)^3
end

function is_inside(e::Ellipse, pos::XYZ)
  px, pz = _eframe(e, pos)
  abs2(px / e.a) + abs2(pz / e.b) < 1
end

function bounding_box(e::Ellipse)
  sθ, cθ = sincos(e.θ)
  hx = hypot(e.a * cθ, e.b * sθ)
  hz = hypot(e.a * sθ, e.b * cθ)
  (min=xyz(e.x - hx, 0, e.z - hz), max=xyz(e.x + hx, 0, e.z + hz))
end

location(e::Ellipse) = xyz(e.x, 0, e.z)

function boundary_projection(e::Ellipse, pos::XYZ)
  # nearest point on the ellipse p(t) = (a cos t, b sin t) to q, in the ellipse frame:
  # Newton on g(t) = (p - q)⋅p′ = (b² - a²) sin t cos t + a qx sin t - b qz cos t
  qx, qz = _eframe(e, pos)
  a = e.a
  b = e.b
  t = atan(a * qz, b * qx)
  # coarse candidates guard against convergence to a non-nearest stationary point
  f = t -> abs2(a * cos(t) - qx) + abs2(b * sin(t) - qz)
  for t₀ ∈ (0:15) .* (π / 8)
    f(t₀) < f(t) && (t = t₀ + 0 * t)
  end
  for _ ∈ 1:32
    s, c = sincos(t)
    g = (b^2 - a^2) * s * c + a * qx * s - b * qz * c
    g′ = (b^2 - a^2) * (c^2 - s^2) + a * qx * c + b * qz * s
    g′ == 0 && break
    Δ = clamp(g / g′, -0.5, 0.5)
    t -= Δ
    abs(Δ) < 1e-13 && break
  end
  u = mod(t / 2π, 1)
  p = boundary_point(e, u)
  n̂ = normal(e, u)
  d = hypot(p.x - pos.x, p.z - pos.z)
  is_inside(e, pos) && (d = -d)
  (u=u, point=p, normal=n̂, curvature=curvature(e, u), distance=d)
end

function intersect_ray(e::Ellipse, origin::XYZ, dir::XYZ)
  # transform the ray to a frame where the ellipse is a unit circle and solve the quadratic
  ox, oz = _eframe(e, origin)
  sθ, cθ = sincos(e.θ)
  dx = cθ * dir.x + sθ * dir.z
  dz = -sθ * dir.x + cθ * dir.z
  ox /= e.a
  dx /= e.a
  oz /= e.b
  dz /= e.b
  A = dx^2 + dz^2
  A == 0 && return nothing
  B = ox * dx + oz * dz
  C = ox^2 + oz^2 - 1
  Δ = B^2 - A * C
  Δ < 0 && return nothing
  t = (-B - √Δ) / A
  t < 0 && (t = (-B + √Δ) / A)
  t < 0 && return nothing
  u = mod(atan(oz + t * dz, ox + t * dx) / 2π, 1)
  (t=t, u=u, point=boundary_point(e, u), normal=normal(e, u), curvature=curvature(e, u))
end

"""
    ParametricCurve(f)

Wrap a function `f(u)` as a 2-D shape. The function must map a boundary
parameter `u ∈ [0,1]` to a position on the boundary of the shape, tracing a
closed curve counter-clockwise in the x–z plane (with `x` rightward and `z`
upward) as `u` increases, so that automatically derived normals point outward.
The position may be given in any form accepted by `xyz()` (e.g. an `(x, z)`
tuple or an `XYZ` named tuple). The function must be differentiable if
normals, curvatures or the generic derived queries are used.

# Examples
```julia-repl
julia> circle = ParametricCurve(u -> (10 * cospi(2u), -40 + 10 * sinpi(2u)));

julia> boundary_point(circle, 0.25)
(x = 0.0, y = 0.0, z = -30.0)
```
"""
struct ParametricCurve{F} <: AbstractShape
  f::F
end

Base.ndims(::ParametricCurve) = 2
boundary_point(s::ParametricCurve, u) = xyz(s.f(u))
Base.show(io::IO, s::ParametricCurve) = print(io, "ParametricCurve($(s.f))")

"""
    ParametricSurface(f)

Wrap a function `f(u, v)` as a 3-D shape. The function must map boundary
parameters `(u, v) ∈ [0,1]²` to a position on the boundary of the shape,
covering the full closed surface with a right-handed parametrization
(`∂p/∂u × ∂p/∂v` pointing outward), so that automatically derived normals
point outward. The position may be given in any form accepted by `xyz()`.
The function must be differentiable if normals are used.
"""
struct ParametricSurface{F} <: AbstractShape
  f::F
end

Base.ndims(::ParametricSurface) = 3
boundary_point(s::ParametricSurface, u, v) = xyz(s.f(u, v))
Base.show(io::IO, s::ParametricSurface) = print(io, "ParametricSurface($(s.f))")

###############################################################################
### scatterers

"""
Superclass for all scatterers. A scatterer is a discrete object in the water
column that scatters sound.
"""
abstract type AbstractScatterer end

"""
    Scatterer(shape; boundary=RigidBoundary)
    Scatterer(shape, boundary)

Create a scatterer with the given `shape` and acoustic `boundary` condition.
A scatterer is a non-penetrable object in the water column, described by its
geometry (an [`AbstractShape`](@ref)) and the acoustic property of its surface
(an `AbstractAcousticBoundary`, e.g. `RigidBoundary` for a sound-hard object
or `PressureReleaseBoundary` for a sound-soft object such as a bubble cloud).

Propagation models that support scatterers are expected to apply the boundary
condition under a local tangent-plane (Kirchhoff) approximation: the reflection
coefficient of the boundary is applied using the local surface normal at each
interaction point. This approximation is valid when the local radius of
curvature of the scatterer is much larger than the acoustic wavelength.

# Examples
```julia-repl
julia> Scatterer(Ellipse(500.0, -40.0, 20.0, 8.0))
Scatterer(Ellipse(x=500.0, z=-40.0, a=20.0, b=8.0), RigidBoundary)

julia> Scatterer(Ellipse(500.0, -40.0, 20.0, 8.0), PressureReleaseBoundary)
Scatterer(Ellipse(x=500.0, z=-40.0, a=20.0, b=8.0), PressureReleaseBoundary)
```
"""
struct Scatterer{G<:AbstractShape,B<:AbstractAcousticBoundary} <: AbstractScatterer
  shape::G
  boundary::B
end

Scatterer(shape::AbstractShape; boundary=RigidBoundary) = Scatterer(shape, boundary)

Base.show(io::IO, s::Scatterer) = print(io, "Scatterer(", s.shape, ", ", s.boundary, ")")

location(s::Scatterer) = location(s.shape)

"""
    has_scatterers(env)

Return `true` if the environment `env` has any scatterers in the water column,
and `false` otherwise.
"""
has_scatterers(env::UnderwaterEnvironment) = !isempty(env.scatterers)

# base number type of the geometry of a scatterer, for env_type()
_geom_type(s::Scatterer) = _geom_type(s.shape)
_geom_type(::Ellipse{T}) where {T} = T
function _geom_type(shape::AbstractShape)
  p = ndims(shape) == 2 ? boundary_point(shape, 0.0) : boundary_point(shape, 0.0, 0.0)
  typeof(p.x)
end
