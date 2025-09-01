import SignalAnalysis: signal, db2pow, tukey
import NonlinearSolve: IntervalNonlinearProblem, solve
import DSP: nextfastfft
import FFTW: ifft

export PekerisRayTracer, PekerisModeSolver

################################################################################
# Pekeris ray tracer

"""
    PekerisRayTracer(env; max_bounces=3)

A fast differentiable ray tracer that only supports iso-velocity constant depth
environments. `max_bounces` is the number of surface/bottom bounces to consider
in the ray tracing.
"""
struct PekerisRayTracer{T1,T2,T3,T4,T5,T6,T7,T8} <: AbstractRayPropagationModel
  h::T1             # water depth
  c::T2             # sound speed
  ρ::T3             # density
  T::T4             # temperature
  S::T5             # salinity
  pH::T6            # pH
  seabed::T7        # seabed properties
  surface::T8       # surface properties
  max_bounces::Int  # maximum number of bounces
  function PekerisRayTracer(env; max_bounces=3)
    max_bounces ≥ 0 || error("Maximum number of bounces cannot be negative")
    is_isovelocity(env) || error("Environment must be iso-velocity")
    is_range_dependent(env) && error("Environment must be range independent")
    is_constant(env.temperature) || error("Temperature must be constant")
    is_constant(env.salinity) || error("Salinity must be constant")
    is_constant(env.pH) || error("pH must be constant")
    is_constant(env.density) || error("Density must be constant")
    h = value(env.bathymetry)
    c = value(env.soundspeed)
    ρ = value(env.density)
    T = value(env.temperature)
    S = value(env.salinity)
    pH = value(env.pH)
    ps = (h, c, ρ, T, S, pH, env.seabed, env.surface)
    new{typeof.(ps)...}(ps..., max_bounces)
  end
end

function Base.show(io::IO, pm::PekerisRayTracer)
  print(io, "PekerisRayTracer(h=$(pm.h), max_bounces=$(pm.max_bounces))")
end

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
    [_arrival(j, pm, R, R², d1, d2, f, typeof(p1), p1, p2) for j ∈ 1:1+2*pm.max_bounces]
  else
    [_arrival(j, pm, R, R², d1, d2, f, typeof(p1)) for j ∈ 1:1+2*pm.max_bounces]
  end
end

"""
    acoustic_field(pm::PekerisRayTracer, tx, rxs; mode=:coherent)

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
    error("Unknown mode :$mode")
  end
end

function impulse_response(pm::AbstractRayPropagationModel, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver, fs; abstime=false, ntaps=nothing)
  arr = arrivals(pm, tx, rx; paths=false)
  T = _phasortype(eltype(arr))
  length(arr) == 0 && return signal(Vector{T}(undef, 0), fs)
  t0, tmax = extrema(a -> a.t, arr)
  # leave 1 ms guard period before and after the impulse response
  t0 = abstime ? zero(t0) : t0 - 0.001
  n = something(ntaps, ceil(Int, (tmax - t0 + 0.001) * fs) + 1)
  signal(_arr2ir([a.t for a ∈ arr], [a.ϕ for a ∈ arr]; T, t0, fs, n), fs)
end

### helpers

_phasortype(::Type{RayArrival{T1,T2,T3,T4,T5}}) where {T1,T2,T3,T4,T5} = Complex{T2}

# complex ForwardDiff friendly version of x^n
_ipow(x, n::Int) = prod(x for _ ∈ 1:n)

# Create an array of length n and type T with ϕs values at times ts, sampled at fs
# and time origin t0. Times ts need not correspond to integer samples.
function _arr2ir(ts, ϕs; T, t0, fs, n)
  x = zeros(T, n)
  for i ∈ eachindex(ts)
    # allocate arrival energy to 2 nearest samples
    t = (ts[i] - t0) * fs + 1
    t̄ = floor(Int, t)
    α, β = sincospi(0.5 * (t - t̄))
    t̄ ≤ n && (x[t̄] += β * ϕs[i])
    t̄ < n && (x[t̄+1] += α * ϕs[i])
  end
  x
end

function _arrival(j, pm, R, R², d1, d2, f, T, p1=missing, p2=missing)
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
  A = Complex(1.0, 0.0) / D * absorption(f, D, pm.S, pm.T, pm.h / 2, pm.pH)  # nominal absorption
  s > 0 && (A *= _ipow(reflection_coef(pm.surface, f, θ, pm.ρ, pm.c), s))
  b > 0 && (A *= _ipow(reflection_coef(pm.seabed, f, θ, pm.ρ, pm.c), b))
  λ = π/2 - θ
  if p1 !== missing
    path = Array{T}(undef, 2 + s + b)
    path[1] = p1
    if s + b > 0
      dx = p2[1] - p1[1]
      dy = p2[2] - p1[2]
      z = (1 - upward) * pm.h
      r = abs(z - d1) * tan(θ)
      path[2] = xyz(p1[1] + r/R * dx, p1[2] + r/R * dy, -z)
      for i ∈ 3:length(path)-1
        r += pm.h * tan(θ)
        z = pm.h - z
        path[i] = xyz(p1[1] + r/R * dx, p1[2] + r/R * dy, -z)
      end
    end
    path[end] = p2
    RayArrival(t, A, s, b, s1 * λ, -s1 * s2 * λ, path)
  else
    RayArrival(t, A, s, b, s1 * λ, -s1 * s2 * λ, T[])
  end
end

################################################################################
# Pekeris mode propagation model
#
# current limitations:
#   iso-velocity, range independent, no bottom absorption, pressure release surface,
#   fluid half-space seabed, no layers, no leaky modes

"""
    PekerisModeSolver(env; ngrid=0, nmodes=0)

A fast differentiable mode propagation model that only supports iso-velocity
constant depth environments.

`ngrid` is the number of grid points to use for modal root finding for fluid
bottom environments. If `ngrid` is too small, the mode solver may miss some
modes. If `ngrid` is too large, the mode solver may take a long time to
converge. The default value of `ngrid` of 0 will use a heuristic to automatically
determine the number of grid points to use.

`nmodes` controls the maximum number of modes computed. Setting `nmodes` to 0
causes all propagating modes to be computed.
"""
struct PekerisModeSolver{T1,T2,T3,T4,T5,T6,T7,T8} <: AbstractModePropagationModel
  h::T1             # water depth
  c::T2             # sound speed
  ρ::T3             # density
  T::T4             # temperature
  S::T5             # salinity
  pH::T6            # pH
  seabed::T7        # seabed properties
  surface::T8       # surface properties
  ngrid::Int        # number of grid points for mode computation
  nmodes::Int       # maximum number of modes (0 for no limit)
  function PekerisModeSolver(env; ngrid=0, nmodes=0)
    is_isovelocity(env) || error("Environment must be iso-velocity")
    is_range_dependent(env) && error("Environment must be range independent")
    is_constant(env.temperature) || error("Temperature must be constant")
    is_constant(env.salinity) || error("Salinity must be constant")
    is_constant(env.pH) || error("pH must be constant")
    is_constant(env.density) || error("Density must be constant")
    env.seabed isa FluidBoundary || error("Seabed must be a fluid boundary")
    env.surface.c == 0 || error("Surface must be a pressure release boundary")
    ngrid == 0 || ngrid > 2 || error("A minimum of 2 grid points are required")
    nmodes ≥ 0 || error("Number of modes cannot be negative")
    h = value(env.bathymetry)
    c = value(env.soundspeed)
    ρ = value(env.density)
    T = value(env.temperature)
    S = value(env.salinity)
    pH = value(env.pH)
    ps = (h, c, ρ, T, S, pH, env.seabed, env.surface)
    new{typeof.(ps)...}(ps..., ngrid, nmodes)
  end
end

function Base.show(io::IO, pm::PekerisModeSolver)
  print(io, "PekerisModeSolver(h=$(pm.h))")
end

function arrivals(pm::PekerisModeSolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver)
  p1 = location(tx)
  p2 = location(rx)
  R = sqrt(abs2(p1.x - p2.x) + abs2(p1.y - p2.y))
  ω = 2π * frequency(tx)
  log_a = log(absorption(frequency(tx), 1.0, pm.S, pm.T, pm.h / 2, pm.pH))  # nominal absorption
  k₁ = ω / pm.c
  if pm.seabed.c == 0                            # pressure release boundary
    M = floor(Int, k₁ * pm.h / π)
    0 < pm.nmodes < M && (M = pm.nmodes)
    return _mode.(1:M, ω, (1:M) .* (π / pm.h), k₁, pm.c, pm.ρ, pm.seabed.c, pm.seabed.ρ, pm.h, log_a)
  elseif isinf(pm.seabed.c)                      # rigid boundary
    M = floor(Int, k₁ * pm.h / π + 0.5)
    0 < pm.nmodes < M && (M = pm.nmodes)
    return _mode.(1:M, ω, ((1:M) .- 0.5) .* (π / pm.h), k₁, pm.c, pm.ρ, pm.seabed.c, pm.seabed.ρ, pm.h, log_a)
  else                                           # acousto-elastic boundary
    k₂ = ω / pm.seabed.c
    if k₂ > k₁
      # dummy computation only to get the correct type
      m = _mode(0, ω, 0.0, k₁, pm.c, pm.ρ, pm.seabed.c, pm.seabed.ρ, pm.h, log_a)
      return Vector{typeof(m)}(undef, 0)
    end
    dk² = k₁^2 - k₂^2
    ngrid = pm.ngrid > 0 ? pm.ngrid : 2 * ceil(Int, pm.h * k₁ / π) + 1
    γgrid = range(0, sqrt(dk²) - sqrt(eps()); length=ngrid)
    cost = map(γ -> _arrivals_cost(γ, (pm, dk²)), γgrid)
    ndx = findall(i -> sign(cost[i+1]) * sign(cost[i]) < 0, 1:length(γgrid)-1)
    0 < pm.nmodes < length(ndx) && (ndx = ndx[1:pm.nmodes])
    γ = [solve(IntervalNonlinearProblem{false}(_arrivals_cost, (γgrid[i], γgrid[i+1]), (pm, dk²))).u for i ∈ ndx]
    return _mode.(1:length(γ), ω, γ, k₁, pm.c, pm.ρ, pm.seabed.c, pm.seabed.ρ, pm.h, log_a)
  end
end

"""
    acoustic_field(pm::PekerisModeSolver, tx, rxs; mode=:coherent)

Compute the acoustic field at a receiver `rxs` due to a transmitter `tx` in the
Pekeris waveguide. The field can be computed incoherently if `mode=:incoherent`.
Otherwise, the field is computed coherently.
"""
function acoustic_field(pm::PekerisModeSolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; mode=:coherent)
  only(acoustic_field(pm, tx, [rx]; mode))
end

function acoustic_field(pm::PekerisModeSolver, tx::AbstractAcousticSource, rxs::AbstractArray{<:AbstractAcousticReceiver}; mode=:coherent)
  mode ∈ (:coherent, :incoherent) || error("Unknown mode :$mode")
  p1 = location(tx)
  # modes don't depend on the receiver, so we can compute based on any receiver
  modes = arrivals(pm, tx, first(rxs))
  kᵣ = [m.kᵣ for m ∈ modes]
  k = (2π * frequency(tx) / pm.c)
  k² = k^2
  γ = sqrt.(k² .- kᵣ.^2)
  map(rxs) do rx
    p2 = location(rx)
    R = sqrt(abs2(p1.x - p2.x))
    modal_terms = @. sin(γ * -p1.z) * sin(γ * -p2.z) * cis(-kᵣ * R) / sqrt(kᵣ)
    multiplier = cispi(-0.25) * 2 / pm.h * sqrt(2π / R)
    if mode === :coherent
      # conjugation required because we use the convention of cis(-kᵣR) to match with Kraken
      conj(sum(modal_terms) * db2amp(spl(tx)) * multiplier)
    else
      conj(complex(√sum(abs2, modal_terms)) * db2amp(spl(tx)) * multiplier)
    end
  end
end

"""
    impulse_response(pm::AbstractModePropagationModel, tx, rx, fs; kwargs...)

Compute the impulse response at the receiver `rx` due to the source `tx` using
propagation model `pm` at the given sampling frequency `fs`.

Several `kwargs` may be specified:

- If `abstime` is `true` (default: `false`), the result is in absolute time
  from the start of transmission. Otherwise, the result is relative to the
  earliest arrival time of the signal at the receiver (with some guard period
  to accommodate acausal response).
- `ntaps` (default: `nothing` for automatic) specifies the number of taps in
  the impulse response.
- `fmin` and `fmax` specifies the bandwidth of interest (default: `0` to `fs/2`).
  If the impulse response is used to convolve with bandlimited signals, it is
  recommended that the impulse response bandwidth be reduced to match the signal
  bandwidth to reduce computational load and improve numerical stability.
- `nmodes` (default: `nothing` for no limit) is the maximum of modes used in
  impulse response computation.
- `threshold` (default: `-60` dB) is used to drop attenuated modes to manage
  computational load.
- `acausal` (default: `0.02` s) controls the computation of acausal impulse
  response (before the estimated time of arrival of first mode).
- `taper` (default: `0.1`) applies a tukey window to the impulse response
  to limit it to the estimated delay spread.

The impulse response is computed only for positive frequencies. Such an impulse
response is suitable for convolution with passband complex analytic signals.
If convolved with real signals, the resulting signal is approximately equivalent
to converting the real signal to a complex analytic form and then convolving it
with the impulse response.
"""
function impulse_response(pm::AbstractModePropagationModel,
  tx::AbstractAcousticSource, rx::AbstractAcousticReceiver, fs;
  abstime=false, ntaps=nothing, fmin=0, fmax=fs/2, nmodes=nothing,
  acausal=0.02, threshold=-60.0, taper=0.1)
  # impulse response computation is designed for Pekeris mode solver, but should
  # work for any mode solver that returns ModeArrival, and hence is marked as
  # a fallback for any AbstractModePropagationModel
  0 ≤ fmin < fs/2 || error("fmin must be positive and less than fs/2")
  fmin < fmax ≤ fs/2 || error("fmax must be between fmin and fs/2")
  arr = arrivals(pm, tx, rx)
  isempty(arr) && return signal(ComplexF64[], fs)
  p1 = location(tx)
  p2 = location(rx)
  R = abs(p1.x - p2.x)
  Δt = R / maximum(a -> a.v, arr) - acausal
  apow = [abs2(cis(-R * a.kᵣ)) for a ∈ arr]
  θ = db2pow(-threshold)
  max_modes = findfirst(1:length(apow)) do i
    sum(apow[1:i]) > θ * sum(apow[i+1:end])
  end
  nmodes = something(nmodes, max_modes)
  nmodes = min(nmodes, max_modes, length(arr))
  mintaps = ceil(Int, (R / arr[nmodes].v - R / arr[1].v) * fs)
  N = nextfastfft(round(Int, (1 + taper) * mintaps))
  Δf = fs / N
  f = (0:N-1) .* Δf
  ndx = findall(fmin .< f .< fmax)
  X = zeros(ComplexF64, N)
  Threads.@threads for i ∈ ndx
    tx1 = AcousticSource(location(tx), f[i+1])
    X[i+1] = conj(acoustic_field(pm, tx1, rx))
  end
  X .*= cispi.(2f * Δt)
  x = √2 * ifft(X)  # factor of √2 to account for energy only in +ve freq
  taper > 0 && (x .*= tukey(length(x), 2 * taper))
  if abstime
    x = vcat(zeros(eltype(x), round(Int, Δt * fs)), x)
  end
  if ntaps !== nothing
    if length(x) ≥ ntaps
      x = x[1:ntaps]
    else
      x = vcat(x, zeros(eltype(x), ntaps - length(x)))
    end
  end
  signal(x, fs)
end

### helpers

function _group_velocity(ω, γ, kᵣ, c, ρ, cb, ρb, D)
  cb == 0 && return real(c * c * kᵣ / ω)
  isinf(cb) && return real(c * c * kᵣ / ω)
  k₁ = ω / c
  kb = ω / cb
  ζ = sqrt(k₁^2 - γ^2 - kb^2)
  sγD = sin(γ * D)
  cγD = cos(γ * D)
  ∂γ = -(c^2 - cb^2) * ρ * ω * sγD / (
    c^2 * cb^2 * D * ρ * cγD * γ^2 +
    c^2 * cb^2 * sγD * γ * (ρ + D * ρb * ζ) +
    cγD * (-cb^2 * D * ρ * ω^2 + c^2 * (D * ρ * ω^2 - cb^2 * ρb * ζ))
  )
  real(1 / (k₁ / (c * kᵣ) - γ / kᵣ * ∂γ))
end

# callable structure representing a mode
struct Mode{T}
  γ::T
  C::T
  zrange::Tuple{T,T}
  function Mode(γ, C, zrange)
    γ, C, z1, z2 = promote(γ, C, zrange...)
    new{typeof(γ)}(γ, C, (z1, z2))
  end
end

(m::Mode)(z) = m.C * sin(m.γ * -z)

function _mode(m, ω, γ, k, c, ρ, cb, ρb, D, log_a)
  kᵣ = complex(sqrt(k^2 - γ^2), log_a)
  v = m > 0 ? _group_velocity(ω, γ, kᵣ, c, ρ, cb, ρb, D) : 0.0
  ModeArrival(m, kᵣ, Mode(γ, sqrt(2/D), (zero(D), -D)), v, ω / real(kᵣ))
end

_arrivals_cost(γ, (pm, dk²)) = pm.seabed.ρ * γ * cos(γ * pm.h) + pm.ρ * sqrt(dk² - γ^2) * sin(γ * pm.h)
