import SignalAnalysis: signal
import Roots: find_zero, Bisection

export PekerisRayTracer, PekerisModeSolver

################################################################################
# Pekeris ray tracer

"""
    PekerisRayTracer(env; nbounces=3)

A fast differentiable ray tracer that only supports isovelocity constant depth
environments. `nbounces` is the number of surface/bottom bounces to consider
in the ray tracing.
"""
struct PekerisRayTracer{T1,T2,T3,T4,T5,T6,T7} <: AbstractRayPropagationModel
  h::T1             # water depth
  c::T2             # sound speed
  ρ::T3             # density
  T::T4             # temperature
  S::T5             # salinity
  seabed::T6        # seabed properties
  surface::T7       # surface properties
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
    ps = (h, c, ρ, T, S, env.seabed, env.surface)
    new{typeof.(ps)...}(ps..., nbounces)
  end
end

function Base.show(io::IO, pm::PekerisRayTracer)
  print(io, "PekerisRayTracer(h=$(pm.h), nbounces=$(pm.nbounces))")
end

"""
    arrivals(pm, tx, rx; paths=true)

Compute the arrivals between a transmitter `tx` and a receiver `rx` in the
Pekeris waveguide described by propagation model `pm`.

If `paths=true`, the eigenray paths are computed and stored in the returned
arrivals. Setting `paths=false` avoids eigenray computation and is slightly
faster.
"""
function arrivals(pm::PekerisRayTracer, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; paths=true)
  # based on Chitre (2007)
  f = frequency(tx)
  p1 = location(tx)
  p2 = location(rx)
  R² = abs2(p1.x - p2.x) + abs2(p1.y - p2.y)
  R = √R²
  d1 = -p1.z
  d2 = -p2.z
  if paths
    [_arrival(j, pm, R, R², d1, d2, f, p1, p2) for j ∈ 1:1+2*pm.nbounces]
  else
    [_arrival(j, pm, R, R², d1, d2, f) for j ∈ 1:1+2*pm.nbounces]
  end
end

"""
    acoustic_field(pm, tx, rxs; mode=:coherent)

Compute the acoustic field at a receiver `rxs` due to a transmitter `tx` in the
Pekeris waveguide. The field is computed incoherently if `mode=:incoherent`.
Otherwise, the field is computed coherently.
"""
function acoustic_field(pm::PekerisRayTracer, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; mode=:coherent)
  arr = arrivals(pm, tx, rx; paths=false)
  length(arr) == 0 && return zero(_phasortype(eltype(arr)))
  if mode === :incoherent
    complex(√sum(a -> abs2(a.ϕ), arr)) * db2amp(spl(tx))
  elseif mode === :coherent
    f = frequency(tx)
    sum(a -> a.ϕ * cispi(2f * a.t), arr) * db2amp(spl(tx))
  else
    error("Unknown mode: $mode")
  end
end

function impulse_response(pm::PekerisRayTracer, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver, fs; abstime=false)
  arr = arrivals(pm, tx, rx; paths=false)
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

function _arrival(j, pm, R, R², d1, d2, f, p1=missing, p2=missing)
  upward = iseven(j)
  s1 = 2 * upward - 1
  n = div(j, 2)
  s = div(n + upward, 2)
  b = div(n + (1 - upward), 2)
  s2 = 2 * iseven(n) - 1
  dz = 2 * b * pm.h + s1 * d1 - s1 * s2 * d2
  D = √(R² + abs2(dz))
  θ = atan(R, dz)
  t = D / pm.c
  A = Complex(1.0, 0.0) / D * absorption(f, D, pm.S)
  s > 0 && (A *= _ipow(reflection_coef(pm.surface, f, θ, pm.ρ, pm.c), s))
  b > 0 && (A *= _ipow(reflection_coef(pm.seabed, f, θ, pm.ρ, pm.c), b))
  λ = π/2 - θ
  if p1 !== missing
    path = Array{typeof(p1)}(undef, 2 + s + b)
    path[1] = p1
    if s + b > 0
      dx = p2[1] - p1[1]
      dy = p2[2] - p1[2]
      z = (1 - upward) * pm.h
      r = abs(z - d1) * tan(θ)
      path[2] = (p1[1] + r/R * dx, p1[2] + r/R * dy, -z)
      for i ∈ 3:length(path)-1
        r += pm.h * tan(θ)
        z = pm.h - z
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

################################################################################
# Pekeris mode propagation model

"""
    PekerisModeSolver(env)

A fast differentiable mode propagation model that only supports isovelocity constant
depth environments.
"""
struct PekerisModeSolver{T1,T2,T3,T4,T5,T6,T7} <: AbstractModePropagationModel
  h::T1             # water depth
  c::T2             # sound speed
  ρ::T3             # density
  T::T4             # temperature
  S::T5             # salinity
  seabed::T6        # seabed properties
  surface::T7       # surface properties
  function PekerisModeSolver(env)
    isospeed(env) || error("Environment must be isovelocity")
    is_range_dependent(env) && error("Environment must be range independent")
    is_constant(env.temperature) || error("Temperature must be constant")
    is_constant(env.salinity) || error("Salinity must be constant")
    is_constant(env.density) || error("Density must be constant")
    env.seabed isa FluidBoundary || error("Seabed must be a fluid boundary")
    # TODO: handle non-pressure release surface
    env.surface.c == 0 || error("Surface must be a pressure release boundary")
    # TODO: handle seabed absorption
    env.seabed.δ == 0 || @warn "Seabed absorption will be ignored"
    h = value(env.bathymetry)
    c = value(env.soundspeed)
    ρ = value(env.density)
    T = value(env.temperature)
    S = value(env.salinity)
    ps = (h, c, ρ, T, S, env.seabed, env.surface)
    new{typeof.(ps)...}(ps...)
  end
end

function Base.show(io::IO, pm::PekerisModeSolver)
  print(io, "PekerisModeSolver(h=$(pm.h))")
end

"""
    arrivals(pm, tx, rx)

Compute the mode arrivals between a transmitter `tx` and a receiver `rx` in the
Pekeris waveguide.
"""
function arrivals(pm::PekerisModeSolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver)
  p1 = location(tx)
  p2 = location(rx)
  R = sqrt(abs2(p1.x - p2.x) + abs2(p1.y - p2.y))
  ω = 2π * frequency(tx)
  k₁ = ω / pm.c
  if pm.seabed.c == 0                            # pressure release boundary
    M = floor(Int, k₁ * pm.h / π)
    m = _mode.(1:M, (1:M) .* (π / pm.h), k₁, R, pm.c)
  elseif isinf(pm.seabed.c)                      # rigid boundary
    M = floor(Int, k₁ * pm.h / π + 0.5)
    m = _mode.(1:M, ((1:M) .- 0.5) .* (π / pm.h), k₁, R, pm.c)
  else                                           # acousto-elastic boundary
    k₂ = ω / pm.seabed.c
    J(γ) = sin(γ*pm.h) + γ / sqrt(k₁^2 - k₂^2 - γ^2) * cos(γ * pm.h)
    # TODO: add support for leaky modes
    k₂ < k₁ || return ModeArrival[]
    γgrid = range(0, sqrt(k₁^2 - k₂^2) - eps(); length=1000)
    ndx = findall(i -> sign(J(γgrid[i+1])) * sign(J(γgrid[i])) < 0, 1:length(γgrid)-1)
    γ = [find_zero(J, (γgrid[i], γgrid[i+1]), Bisection()) for i ∈ ndx]
    m = _mode.(1:length(γ), γ, k₁, R, pm.c)
  end
  m
end

"""
    acoustic_field(pm, tx, rxs; mode=:coherent)

Compute the acoustic field at a receiver `rxs` due to a transmitter `tx` in the
Pekeris waveguide. The field can be computed incoherently if `mode=:incoherent`.
Otherwise, the field is computed coherently.
"""
function acoustic_field(pm::PekerisModeSolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; mode=:coherent)
  only(acoustic_field(pm, tx, [rx]; mode))
end

function acoustic_field(pm::PekerisModeSolver, tx::AbstractAcousticSource, rxs::AbstractArray{<:AbstractAcousticReceiver}; mode=:coherent)
  mode ∈ (:coherent, :incoherent) || error("Unknown mode: $mode")
  p1 = location(tx)
  ρₛ = value(pm.ρ, p1)
  # modes don't depend on the receiver, so we can compute based on any receiver
  modes = arrivals(pm, tx, rxs[1])
  γ = [m.kz for m ∈ modes]
  kᵣ = [m.kr for m ∈ modes]
  ϕ = [m.ϕ for m ∈ modes]
  tmap(rxs) do rx
    p2 = location(rx)
    R = sqrt(abs2(p1.x - p2.x) + abs2(p1.y - p2.y))
    modal_terms = @. ϕ * sin(γ * -p1.z) * sin(γ * -p2.z) * hankelh1(0, kᵣ * R)
    if mode === :coherent
      sum(modal_terms) * im / (2 * pm.h * ρₛ) * db2amp(spl(tx))
    else
      complex(√sum(abs2, modal_terms) / (2 * pm.h * ρₛ)) * db2amp(spl(tx))
    end
  end
end

function impulse_response(pm::PekerisModeSolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver, fs; abstime=false)

end

### helpers

function _mode(m, γ, k, R, c)
  kᵣ = sqrt(k^2 - γ^2)
  cᵣ = kᵣ / k * c
  t = R / cᵣ                      # FIXME: this is not an accurate way to compute this
  θ = asin(γ / k)
  ϕ = complex(1.0, 0.0)
  ModeArrival(t, m, γ, kᵣ, θ, ϕ)
end
