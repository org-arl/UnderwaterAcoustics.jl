export Reframe2D

################################################################################
# Reframe2D — automatic reframing of 2D scenarios into a model's canonical
# coordinate system
#
# Many 2D propagation models require scenarios to be specified in a canonical
# coordinate system with the source at (0, 0, z) and receivers in the x-z plane
# with x ≥ 0. These are implementation details rather than fundamental
# limitations. When a scenario is genuinely 2D (source and receivers lie in a
# single vertical plane), it can be rigidly transformed (translation + rotation
# about the z-axis) into the canonical frame, the model run there, and the
# results transformed back to world coordinates. The environment is guaranteed
# to be correct only on the vertical plane containing the source and receivers,
# as that is the only region a 2D model queries.

"""
    Reframe2D(model, env; atol=0.1, kwargs...)

Wrap a 2D propagation `model` so that scenarios may be specified in world
coordinates rather than the model's canonical coordinate system. Any `kwargs`
passed in are transferred to the underlying `model` when it is constructed.

Many 2D propagation models require the source to be located at `(0, 0, z)` and
all receivers to lie in the x-z plane with `x ≥ 0`. The wrapper automatically
translates and rotates (about the z-axis) the source, receivers, and
environment (bathymetry, altimetry, sound speed profile, etc) into the
coordinate system required by the model, and transforms results (e.g. eigenray
paths) back to the original coordinates.

A scenario is considered 2D if all receivers lie within a horizontal distance
`atol` (in meters) of the vertical plane containing the source and the farthest
receiver. Receivers within `atol` of the plane are projected onto it. A single
receiver anywhere is always a valid 2D scenario. If the scenario is already in
canonical coordinates, the wrapper passes it through unchanged.

The transformed environment is only guaranteed to match the original
environment on the vertical plane containing the source and receivers.

Scatterers in the water column are transformed along with the environment.
Since 2-D scatterer shapes lie in the x-z plane, scenarios with 2-D scatterers
are considered 2D only if the propagation plane coincides with the x-z plane
(i.e. the scenario is at most translated along the x-axis); 3-D scatterer
shapes support any scenario orientation.

# Examples
```julia-repl
julia> pm = Reframe2D(PekerisRayTracer, UnderwaterEnvironment(bathymetry=20.0))
Reframe2D(PekerisRayTracer)

julia> tx = AcousticSource((x=100.0, y=-50.0, z=-5.0), 1000.0);

julia> rx = AcousticReceiver((x=500.0, y=250.0, z=-10.0));

julia> transmission_loss(pm, tx, rx)
```
"""
function Reframe2D end

struct Reframe2DRay{T1,T2,M} <: AbstractRayPropagationModel
  env::T1
  kwargs::T2
  atol::Float64
  function Reframe2DRay(model::Type{<:AbstractRayPropagationModel}, env, atol, kwargs)
    new{typeof(env),typeof(kwargs),model}(env, kwargs, atol)
  end
end

struct Reframe2DMode{T1,T2,M} <: AbstractModePropagationModel
  env::T1
  kwargs::T2
  atol::Float64
  function Reframe2DMode(model::Type{<:AbstractModePropagationModel}, env, atol, kwargs)
    new{typeof(env),typeof(kwargs),model}(env, kwargs, atol)
  end
end

const AbstractReframe2D = Union{Reframe2DRay,Reframe2DMode}

function Reframe2D(model::Type{<:AbstractRayPropagationModel}, env; atol=0.1, kwargs...)
  Reframe2DRay(model, env, Float64(in_units(u"m", atol)), kwargs)
end

function Reframe2D(model::Type{<:AbstractModePropagationModel}, env; atol=0.1, kwargs...)
  Reframe2DMode(model, env, Float64(in_units(u"m", atol)), kwargs)
end

_model(pm::Reframe2DRay{T1,T2,M}) where {T1,T2,M} = M
_model(pm::Reframe2DMode{T1,T2,M}) where {T1,T2,M} = M

Base.show(io::IO, pm::AbstractReframe2D) = print(io, "Reframe2D($(_model(pm)))")

## interface methods
#
# methods are defined on the concrete wrapper types rather than on
# AbstractReframe2D, as methods on the Union would be ambiguous against
# methods on AbstractRayPropagationModel or AbstractModePropagationModel
# (e.g. the impulse_response fallbacks)

function _arrivals2d(pm, tx, rx; kwargs...)
  t, pm1, tx1, rx1 = _setup2d(pm, tx, rx)
  _unreframe(t, arrivals(pm1, tx1, rx1; kwargs...))
end

function _acoustic_field2d(pm, tx, rxs; kwargs...)
  _, pm1, tx1, rxs1 = _setup2d(pm, tx, rxs)
  acoustic_field(pm1, tx1, rxs1; kwargs...)
end

function _impulse_response2d(pm, tx, rx, fs; kwargs...)
  _, pm1, tx1, rx1 = _setup2d(pm, tx, rx)
  impulse_response(pm1, tx1, rx1, fs; kwargs...)
end

arrivals(pm::Reframe2DRay, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; kwargs...) = _arrivals2d(pm, tx, rx; kwargs...)
arrivals(pm::Reframe2DMode, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; kwargs...) = _arrivals2d(pm, tx, rx; kwargs...)
acoustic_field(pm::Reframe2DRay, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; kwargs...) = _acoustic_field2d(pm, tx, rx; kwargs...)
acoustic_field(pm::Reframe2DMode, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; kwargs...) = _acoustic_field2d(pm, tx, rx; kwargs...)
acoustic_field(pm::Reframe2DRay, tx::AbstractAcousticSource, rxs::AbstractArray{<:AbstractAcousticReceiver}; kwargs...) = _acoustic_field2d(pm, tx, rxs; kwargs...)
acoustic_field(pm::Reframe2DMode, tx::AbstractAcousticSource, rxs::AbstractArray{<:AbstractAcousticReceiver}; kwargs...) = _acoustic_field2d(pm, tx, rxs; kwargs...)
impulse_response(pm::Reframe2DRay, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver, fs; kwargs...) = _impulse_response2d(pm, tx, rx, fs; kwargs...)
impulse_response(pm::Reframe2DMode, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver, fs; kwargs...) = _impulse_response2d(pm, tx, rx, fs; kwargs...)

## rigid transform

# rigid horizontal transform: world → canonical is a translation by (-x₀, -y₀)
# followed by a rotation such that the world direction (cosθ, sinθ) maps to
# the canonical +x axis; z is unaffected
struct RigidXform{T}
  x₀::T
  y₀::T
  cosθ::T
  sinθ::T
end

function _forward(t::RigidXform, p::XYZ)
  dx = p.x - t.x₀
  dy = p.y - t.y₀
  xyz(dx * t.cosθ + dy * t.sinθ, dy * t.cosθ - dx * t.sinθ, p.z)
end

function _inverse(t::RigidXform, p::XYZ)
  xyz(t.x₀ + p.x * t.cosθ - p.y * t.sinθ, t.y₀ + p.x * t.sinθ + p.y * t.cosθ, p.z)
end

# only exact floating point frames may be treated as identity, so that
# transforms carrying AD number types (e.g. ForwardDiff duals) always take the
# full transform path and propagate their partial derivatives
_isidentity(t::RigidXform{<:Union{Float32,Float64}}) = t.x₀ == 0 && t.y₀ == 0 && t.sinθ == 0 && t.cosθ == 1
_isidentity(::RigidXform) = false

## transform derivation and validation

function _txpos(tx)
  p₀ = location(tx)
  (p₀ === nothing || p₀ === missing) && error("Reframe2D requires a source with a known location")
  p₀
end

# derive transform from horizontal offsets Δ of all receivers relative to the
# source at p₀, validating that the scenario is 2D
function _xform(pm::AbstractReframe2D, p₀, Δ)
  r = map(d -> hypot(d[1], d[2]), Δ)
  i = argmax(r)
  # bearing towards the farthest receiver; if all receivers are (numerically)
  # directly above/below the source, the bearing is undefined and we default
  # to the identity rotation
  cosθ, sinθ = r[i] ≤ pm.atol ? (one(p₀.x), zero(p₀.x)) : (Δ[i][1] / r[i], Δ[i][2] / r[i])
  dev = maximum(d -> abs(d[2] * cosθ - d[1] * sinθ), Δ)
  dev ≤ pm.atol || error("Scenario is not 2D: receivers deviate up to " *
    "$(round(dev; sigdigits=3)) m from the vertical plane containing the source " *
    "(atol = $(pm.atol) m); use a 3D propagation model, compute the field one " *
    "receiver at a time, or increase `atol` if the deviation is physically negligible")
  RigidXform(promote(p₀.x, p₀.y, cosθ, sinθ)...)
end

function _derive_xform(pm::AbstractReframe2D, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver)
  p₀ = _txpos(tx)
  p = location(rx)
  _xform(pm, p₀, ((p.x - p₀.x, p.y - p₀.y),))
end

function _derive_xform(pm::AbstractReframe2D, tx::AbstractAcousticSource, rxs::AbstractArray{<:AbstractAcousticReceiver})
  p₀ = _txpos(tx)
  _xform(pm, p₀, vec(map(rx -> (location(rx).x - p₀.x, location(rx).y - p₀.y), rxs)))
end

# a receiver grid spans the y = 0 plane, so the only rigid transforms that
# keep it a grid are a translation along x, possibly with a 180° rotation;
# the source must therefore lie in the y = 0 plane (unless the grid is a
# single vertical column of receivers)
function _derive_xform(pm::AbstractReframe2D, tx::AbstractAcousticSource, rxs::AcousticReceiverGrid2D)
  p₀ = _txpos(tx)
  rxs.xrange.len == 1 && return _xform(pm, p₀, ((rxs.xrange[1] - p₀.x, -p₀.y),))
  abs(p₀.y) ≤ pm.atol || error("Scenario is not 2D: source (y = $(p₀.y)) does not " *
    "lie in the y = 0 plane of the receiver grid (atol = $(pm.atol) m)")
  x1, x2 = extrema(rxs.xrange)
  far = abs(x2 - p₀.x) ≥ abs(x1 - p₀.x) ? x2 : x1
  s = far ≥ p₀.x ? one(p₀.x) : -one(p₀.x)
  RigidXform(promote(p₀.x, p₀.y, s, zero(s))...)
end

## reframing of scenario components

function _setup2d(pm::AbstractReframe2D, tx, rxs)
  t = _derive_xform(pm, tx, rxs)
  _check_scatterers(pm, t)
  pm1 = _model(pm)(_reframe(t, pm.env); pm.kwargs...)
  if _isidentity(t) && _inplane(rxs)
    (t, pm1, tx, rxs)
  else
    (t, pm1, _reframe(t, tx), _reframe(t, rxs))
  end
end

_inplane(rx::AbstractAcousticReceiver) = location(rx).y == 0
_inplane(rxs::AcousticReceiverGrid2D) = true
_inplane(rxs::AbstractArray) = all(rx -> location(rx).y == 0, rxs)

function _reframe(t::RigidXform, tx::AbstractAcousticSource)
  p = location(tx)
  AcousticSource(xyz(zero(p.x), zero(p.y), p.z), frequency(tx); spl=spl(tx))
end

function _reframe(t::RigidXform, rx::AbstractAcousticReceiver)
  p = _forward(t, location(rx))
  # project exactly onto the x-z plane, as many models require y == 0
  AcousticReceiver(xyz(p.x, zero(p.y), p.z))
end

_reframe(t::RigidXform, rxs::AbstractArray{<:AbstractAcousticReceiver}) = map(rx -> _reframe(t, rx), rxs)

function _reframe(t::RigidXform, rxs::AcousticReceiverGrid2D)
  if t.sinθ == 0 && abs(t.cosθ) == 1
    xr = t.cosθ > 0 ? rxs.xrange .- t.x₀ : t.x₀ .- rxs.xrange
    AcousticReceiverGrid2D(xr, rxs.zrange)
  else
    # single vertical column of receivers (guaranteed by _derive_xform)
    p = _forward(t, location(rxs[1,1]))
    AcousticReceiverGrid2D(p.x, rxs.zrange)
  end
end

## reframing of environment

function _reframe(t::RigidXform, env::UnderwaterEnvironment)
  _isidentity(t) && return env
  fields = (env.bathymetry, env.altimetry, env.temperature, env.salinity,
    env.pH, env.soundspeed, env.density, env.seabed, env.surface)
  fields1 = map(q -> _reframe_field(t, q), fields)
  scatterers1 = map(s -> _reframe(t, s), env.scatterers)
  (all(map(===, fields1, fields)) && isempty(scatterers1)) && return env
  UnderwaterEnvironment(fields1..., scatterers1)
end

# 2-D scatterer shapes lie in the x-z plane, so the propagation plane of a
# genuinely 2D scenario must coincide with it (3-D shapes have no such
# restriction, as they transform rigidly in 3D)
function _check_scatterers(pm::AbstractReframe2D, t::RigidXform)
  (_isidentity(t) || !has_scatterers(pm.env)) && return nothing
  any(_is2dscatterer, pm.env.scatterers) || return nothing
  (t.sinθ == 0 && abs(t.y₀) ≤ pm.atol) || error("Scenario is not 2D: 2-D scatterers " *
    "lie in the x-z plane, but the vertical plane through the source and receivers " *
    "does not coincide with it")
  nothing
end

_is2dscatterer(s::Scatterer) = ndims(s.shape) == 2
_is2dscatterer(s) = true    # conservative for unknown scatterer types

_reframe(t::RigidXform, s::Scatterer) = Scatterer(_reframe_shape(t, s.shape), s.boundary)

function _reframe_shape(t::RigidXform, shape::AbstractShape)
  if ndims(shape) == 3
    return ParametricSurface(let shape = shape, t = t
      (u, v) -> _forward(t, boundary_point(shape, u, v))
    end)
  end
  t.sinθ == 0 || error("Cannot reframe a 2-D scatterer onto a rotated propagation plane")
  if t.cosθ > 0
    ParametricCurve(let shape = shape, x₀ = t.x₀
      u -> (p = boundary_point(shape, u); (p.x - x₀, p.z))
    end)
  else
    # a 180° rotation mirrors the x-z plane, flipping the traversal direction;
    # reverse the parametrization to keep the boundary counter-clockwise so
    # that derived normals point outward
    ParametricCurve(let shape = shape, x₀ = t.x₀
      u -> (p = boundary_point(shape, 1 - u); (x₀ - p.x, p.z))
    end)
  end
end

function _reframe_shape(t::RigidXform, e::Ellipse)
  t.sinθ == 0 || error("Cannot reframe a 2-D scatterer onto a rotated propagation plane")
  t.cosθ > 0 ? Ellipse(e.x - t.x₀, e.z, e.a, e.b, e.θ) : Ellipse(t.x₀ - e.x, e.z, e.a, e.b, -e.θ)
end

"""
Position-dependent field re-expressed in a different coordinate frame. The
field is evaluated by transforming the position back to the original frame.
"""
struct TransformedField{T1,T2} <: PositionDependent
  q::T1
  t::RigidXform{T2}
end

(w::TransformedField)(pos::XYZ) = w.q(_inverse(w.t, pos))
Base.minimum(w::TransformedField) = minimum(w.q)   # rigid transform preserves extrema
Base.maximum(w::TransformedField) = maximum(w.q)
Base.show(io::IO, w::TransformedField) = print(io, "TransformedField($(w.q))")

# position-independent quantities and depth-only fields are unaffected by a
# horizontal rigid transform
_reframe_field(t::RigidXform, q) = q
_reframe_field(t::RigidXform, q::DepthDependent) = q
_reframe_field(t::RigidXform, q::PositionDependent) = TransformedField(q, t)

# tolerance below which the track is considered perpendicular to the x-axis
const _COSθ_MIN = 1e-8

# fields varying with x only remain x-varying along the canonical x-axis, with
# the sample grid remapped by x′ = (x - x₀) / cosθ; the wrapped function is
# exact everywhere on the canonical x-axis
function _reframe_field(t::RigidXform, q::SampledFieldX)
  abs(t.cosθ) ≤ _COSθ_MIN && return q.f(t.x₀)   # constant along track
  xr = (q.xrange .- t.x₀) ./ t.cosθ
  t.cosθ < 0 && (xr = reverse(xr))
  SampledFieldX(let f = q.f, x₀ = t.x₀, cosθ = t.cosθ
    x -> f(x₀ + x * cosθ)
  end, xr, q.interp)
end

function _reframe_field(t::RigidXform, q::SampledFieldXZ)
  if abs(t.cosθ) ≤ _COSθ_MIN
    return SampledField([q.f(t.x₀, z) for z ∈ q.zrange]; z=collect(q.zrange), interp=q.interp)
  end
  xr = collect((q.xrange .- t.x₀) ./ t.cosθ)
  t.cosθ < 0 && (xr = reverse(xr))
  SampledFieldXZ(let f = q.f, x₀ = t.x₀, cosθ = t.cosθ
    (x, z) -> f(x₀ + x * cosθ, z)
  end, xr, collect(q.zrange), q.interp)
end

# fields varying with x and y are sliced along the track; the wrapped function
# is exact on the canonical x-axis, and the sample grid is chosen to cover the
# projection of the field's sampled region onto the track
function _reframe_field(t::RigidXform, q::SampledFieldXY)
  xr = _track_nodes(t, q.xrange, q.yrange)
  SampledFieldX(let f = q.f, t = t
    x -> f(t.x₀ + x * t.cosθ, t.y₀ + x * t.sinθ)
  end, xr, q.interp)
end

function _reframe_field(t::RigidXform, q::SampledFieldXYZ)
  xr = collect(_track_nodes(t, q.xrange, q.yrange))
  SampledFieldXZ(let f = q.f, t = t
    (x, z) -> f(t.x₀ + x * t.cosθ, t.y₀ + x * t.sinθ, z)
  end, xr, collect(q.zrange), q.interp)
end

function _track_nodes(t::RigidXform, xrange, yrange)
  a, b = extrema((x - t.x₀) * t.cosθ + (y - t.y₀) * t.sinθ for x ∈ extrema(xrange), y ∈ extrema(yrange))
  Δ = min(_minspacing(xrange), _minspacing(yrange))
  Δ > 0 || (Δ = one(Δ))
  b - a ≤ Δ && ((a, b) = (a - Δ, b + Δ))
  range(a, b; length=min(ceil(Int, (b - a) / Δ) + 1, 10001))
end

_minspacing(r::AbstractRange) = abs(step(r))
_minspacing(v) = length(v) > 1 ? minimum(diff(v)) : zero(float(eltype(v)))

## back-transformation of results

_unreframe(t::RigidXform, arr::AbstractVector) = _isidentity(t) ? arr : map(a -> _unreframe(t, a), arr)
_unreframe(t::RigidXform, a::ModeArrival) = a

function _unreframe(t::RigidXform, a::RayArrival)
  (a.path === missing || length(a.path) == 0) && return a
  RayArrival(a.t, a.ϕ, a.ns, a.nb, a.θₛ, a.θᵣ, [_inverse(t, xyz(p)) for p ∈ a.path])
end
