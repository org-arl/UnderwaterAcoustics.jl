import SignalAnalysis: signal

export PekerisRayTracer

"""
    PekerisRayTracer(env; nbounces=3)

A fast differentiable ray tracer that only supports isovelocity constant depth
environments. `nbounces` is the number of surface/bottom bounces to consider
in the ray tracing.
"""
struct PekerisRayTracer{T1,T2,T3,T4,T5,T6,T7,T8} <: AbstractPropagationModel
  h::T1             # water depth
  c::T2             # sound speed
  ρ::T3             # density
  T::T4             # temperature
  S::T5             # salinity
  seabed::T6        # seabed properties
  surface::T7       # surface properties
  noise::T8         # noise model
  nbounces::Int
  function PekerisRayTracer(env; nbounces=16)
    nbounces ≥ 0 || error("Number of nbounces cannot be negative")
    isospeed(env) || error("Environment must be isovelocity")
    is_range_dependent(env) && error("Environment must be range independent")
    is_constant(env.temperature) || error("Temperature must be constant")
    is_constant(env.salinity) || error("Salinity must be constant")
    is_constant(env.density) || error("Density must be constant")
    h = value(env.bathymetry)
    c = value(env.soundspeed)
    ρ = value(env.density)
    T = value(env.temperature)
    S = value(env.salinity)
    ps = (h, c, ρ, T, S, env.seabed, env.surface, env.noise)
    new{typeof.(ps)...}(ps..., nbounces)
  end
end

function Base.show(io::IO, model::PekerisRayTracer)
  print(io, "PekerisRayTracer(h=$(model.h), nbounces=$(model.nbounces))")
end

"""
    arrivals(model, tx, rx; paths=true)

Compute the arrivals between a transmitter `tx` and a receiver `rx` in the
Pekeris waveguide.

If `paths=true`, the eigenray paths are computed and stored in the returned
arrivals. Setting `paths=false` avoids eigenray computation and is slightly
faster.
"""
function arrivals(model::PekerisRayTracer, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; paths=true)
  # based on Chitre (2007)
  f = frequency(tx)
  p1 = location(tx)
  p2 = location(rx)
  R² = abs2(p1.x - p2.x) + abs2(p1.y - p2.y)
  R = √R²
  d1 = -p1.z
  d2 = -p2.z
  if paths
    [_arrival(j, model, R, R², d1, d2, f, p1, p2) for j ∈ 1:1+2*model.nbounces]
  else
    [_arrival(j, model, R, R², d1, d2, f) for j ∈ 1:1+2*model.nbounces]
  end
end

"""
    acoustic_field(model, tx, rxs; mode=:coherent)

Compute the acoustic field at a receiver `rxs` due to a transmitter `tx` in the
Pekeris waveguide. The field can be computed incoherently if `mode=:incoherent`.
Otherwise, the field is computed coherently.
"""
function acoustic_field(model::PekerisRayTracer, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; mode=:coherent)
  arr = arrivals(model, tx, rx; paths=false)
  length(arr) == 0 && return zero(_phasortype(eltype(arr)))
  if mode === :incoherent
    Complex(√sum(a -> abs2(a.ϕ), arr), 0) * db2amp(spl(tx))
  else
    f = frequency(tx)
    sum(a -> a.ϕ * cispi(2f * a.t), arr) * db2amp(spl(tx))
  end
end

function acoustic_field(model::PekerisRayTracer, tx::AbstractAcousticSource, rxs::AbstractArray{<:AbstractAcousticReceiver}; mode=:coherent)
  tmap(rx -> acoustic_field(model, tx, rx; mode), rxs)
end

function impulse_response(model::PekerisRayTracer, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver, fs; abstime=false)
  arr = arrivals(model, tx, rx; paths=false)
  T = _phasortype(eltype(arr))
  length(arr) == 0 && return Vector{T}(undef, 0)
  t0, tmax = extrema(a -> a.t, arr)
  abstime && (t0 = zero(t0))
  n = ceil(Int, (tmax - t0) * fs) + 1
  x = zeros(T, n)
  for a ∈ arr
    # allocate arrival energy to 2 nearest samples
    t = (a.t - t0) * fs + 1
    t̄ = floor(Int, t)
    α = sqrt(t - t̄)
    β = sqrt(1 - t + t̄)
    x[t̄] += β * a.ϕ
    x[t̄+1] += α * a.ϕ
  end
  signal(x, fs)
end

### helpers

_phasortype(::Type{RayArrival{T1,T2,T3,T4,T5}}) where {T1,T2,T3,T4,T5} = Complex{T2}

# complex ForwardDiff friendly version of x^n
_ipow(x, n::Int) = prod(x for _ ∈ 1:n)

function _arrival(j, model, R, R², d1, d2, f, p1=missing, p2=missing)
  upward = iseven(j)
  s1 = 2 * upward - 1
  n = div(j, 2)
  s = div(n + upward, 2)
  b = div(n + (1 - upward), 2)
  s2 = 2 * iseven(n) - 1
  dz = 2 * b * model.h + s1 * d1 - s1 * s2 * d2
  D = √(R² + abs2(dz))
  θ = atan(R, dz)
  t = D / model.c
  A = Complex(1.0, 0.0) / D * absorption(f, D, model.S)
  s > 0 && (A *= _ipow(reflection_coef(model.surface, f, θ, model.c, model.ρ), s))
  b > 0 && (A *= _ipow(reflection_coef(model.seabed, f, θ, model.c, model.ρ), b))
  λ = π/2 - θ
  if p1 !== missing
    path = Array{typeof(p1)}(undef, 2 + s + b)
    path[1] = p1
    if s + b > 0
      dx = p2[1] - p1[1]
      dy = p2[2] - p1[2]
      z = (1 - upward) * model.h
      r = abs(z - d1) * tan(θ)
      path[2] = (p1[1] + r/R * dx, p1[2] + r/R * dy, -z)
      for i ∈ 3:length(path)-1
        r += model.h * tan(θ)
        z = model.h - z
        path[i] = (p1[1] + r/R * dx, p1[2] + r/R * dy, -z)
      end
    end
    path[end] = p2
    # conj(A) needed to match with Bellhop
    RayArrival(t, conj(A), s, b, s1 * λ, -s1 * s2 * λ, path)
  else
    RayArrival(t, conj(A), s, b, s1*λ, -s1*s2*λ, missing)
  end
end
