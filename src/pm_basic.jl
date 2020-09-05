using Interpolations
using SignalAnalysis
using DSP: hilbert

export IsoSSP, MunkSSP, SampledSSP, ConstantDepth, SampledDepth, SampledAltitude
export ReflectionCoef, FlatSurface, Rayleigh, SurfaceLoss
export Rock, Pebbles, SandyGravel, CoarseSand, MediumSand, FineSand, VeryFineSand
export ClayeySand, CoarseSilt, SandySilt, Silt, FineSilt, SandyClay, SiltyClay, Clay
export Vacuum, SeaState0, SeaState1, SeaState2, SeaState3, SeaState4
export SeaState5, SeaState6, SeaState7, SeaState8, SeaState9
export AcousticReceiverGrid2D, AcousticReceiverGrid3D, NarrowbandAcousticSource
export RedGaussianNoise, Pinger

### sound speed profiles

"""
$(TYPEDEF)
Isovelocity sound speed profile.

---

    IsoSSP(c)

Create an isovelocity sound speed profile with sound speed `c`.
"""
struct IsoSSP{T} <: SoundSpeedProfile
  c::T
end

soundspeed(ssp::IsoSSP, x, y, z) = ssp.c

"""
$(TYPEDEF)
Munk sound speed profile.
"""
struct MunkSSP <: SoundSpeedProfile end

function soundspeed(::MunkSSP, x, y, z)
  ϵ = 0.00737
  z̃ = 2.0 * (-z - 1300.0) / 1300.0
  1500.0 * (1.0 + ϵ * (z̃ - 1 + exp(-z̃)))
end

"""
$(TYPEDEF)
Sound speed profile based on measurements at discrete depths.
"""
struct SampledSSP{T1,T2,T3} <: SoundSpeedProfile
  z::Vector{T1}
  c::Vector{T2}
  interp::Symbol
  f::T3
  function SampledSSP(depth, c, interp)
    if interp === :smooth
      depth isa AbstractRange || throw(ArgumentError("depth must be sampled uniformly and specified as an AbstractRange"))
      f = CubicSplineInterpolation(depth, c; extrapolation_bc=Line())
    elseif interp === :linear
      f = LinearInterpolation(depth, c; extrapolation_bc=Line())
    else
      throw(ArgumentError("Unknown interpolation"))
    end
    new{eltype(depth),eltype(c),typeof(f)}(-depth, c, interp, f)
  end
end

"""
    SampledSSP(depth, c)
    SampledSSP(depth, c, interp)

Create a sound speed profile based on measurements at discrete depths.
`interp` may be either `:linear` or `:smooth`, and defaults to `:linear`
if unspecified.
"""
SampledSSP(depth, c) = SampledSSP(depth, c, :linear)

soundspeed(ssp::SampledSSP, x, y, z) = ssp.f(-z)

function Base.show(io::IO, ssp::SampledSSP{T1,T2,T3}) where {T1,T2,T3}
  print(io, "SampledSSP{", T1, ",", T2, ",", ssp.interp, "}(", length(ssp.z), " points)")
end

### bathymetry models

"""
$(TYPEDEF)
Bathymetry with constant depth.

---

    ConstantDepth(depth)

Create a constant depth bathymetry.
"""
struct ConstantDepth{T} <: Bathymetry
  depth::T
end

depth(bathy::ConstantDepth, x, y) = bathy.depth
maxdepth(bathy::ConstantDepth) = bathy.depth

"""
$(TYPEDEF)
Bathymetry based on depth samples.
"""
struct SampledDepth{T1,T2,T3} <: Bathymetry
  x::Vector{T1}
  depth::Vector{T2}
  interp::Symbol
  f::T3
  function SampledDepth(x, depth, interp)
    if interp === :smooth
      x isa AbstractRange || throw(ArgumentError("x must be sampled uniformly and specified as an AbstractRange"))
      f = CubicSplineInterpolation(x, depth; extrapolation_bc=Line())
    elseif interp === :linear
      f = LinearInterpolation(x, depth; extrapolation_bc=Line())
    else
      throw(ArgumentError("Unknown interpolation"))
    end
    new{eltype(x),eltype(depth),typeof(f)}(x, depth, interp, f)
  end
end

"""
    SampledDepth(x, depth)
    SampledDepth(x, depth, interp)

Create a bathymetry given discrete depth measurements at locations given in `x`.
`interp` may be either `:linear` or `:smooth`, and defaults to `:linear`
if unspecified.
"""
SampledDepth(x, depth) = SampledDepth(x, depth, :linear)

depth(bathy::SampledDepth, x, y) = bathy.f(x)
maxdepth(bathy::SampledDepth) = maximum(bathy.depth)

function Base.show(io::IO, b::SampledDepth{T1,T2,T3}) where {T1,T2,T3}
  print(io, "SampledDepth{", T1, ",", T2, ",", b.interp, "}(", length(b.x), " points)")
end

### altimetry models

"""
$(TYPEDEF)
Altimetry for a flat surface with constant altitude of zero.
"""
struct FlatSurface <: Altimetry end

altitude(::FlatSurface, x, y) = 0.0

"""
$(TYPEDEF)
Altimetry based on altitude samples.
"""
struct SampledAltitude{T1,T2,T3} <: Altimetry
  x::Vector{T1}
  altitude::Vector{T2}
  interp::Symbol
  f::T3
  function SampledAltitude(x, altitude, interp)
    if interp === :smooth
      x isa AbstractRange || throw(ArgumentError("x must be sampled uniformly and specified as an AbstractRange"))
      f = CubicSplineInterpolation(x, altitude; extrapolation_bc=Line())
    elseif interp === :linear
      f = LinearInterpolation(x, altitude; extrapolation_bc=Line())
    else
      throw(ArgumentError("Unknown interpolation"))
    end
    new{eltype(x),eltype(altitude),typeof(f)}(x, altitude, interp, f)
  end
end

"""
    SampledAltitude(x, altitude)
    SampledAltitude(x, altitude, interp)

Create an altimetry given discrete altitude measurements at locations given in `x`.
`interp` may be either `:linear` or `:smooth`, and defaults to `:linear`
if unspecified.
"""
SampledAltitude(x, altitude) = SampledAltitude(x, altitude, :linear)

altitude(a::SampledAltitude, x, y) = a.f(x)

function Base.show(io::IO, a::SampledAltitude{T1,T2,T3}) where {T1,T2,T3}
  print(io, "SampledAltitude{", T1, ",", T2, ",", a.interp, "}(", length(a.x), " points)")
end

### reflection models

"""
$(TYPEDEF)
Reflection model for a surface with a constant reflection coefficient.

---

    ReflectionCoef(coef)

Create a reflection model for a surface with a constant reflection coefficient `coef`.
"""
struct ReflectionCoef{T<:Number} <: ReflectionModel
  coef::T
end

reflectioncoef(rm::ReflectionCoef, f, θ) = rm.coef

"""
$(TYPEDEF)
Reflection model for a surface with a Rayleigh reflection coefficient.

---

    Rayleigh(ρᵣ, cᵣ, δ)
    Rayleigh(ρᵣ, cᵣ)

Create a reflection model for a surface with a Rayleigh reflection coefficient
with relative density `ρᵣ`, relative sound speed `cᵣ`, and dimensionless
attentuation `δ`. If attentuation `δ` is unspecified, it is assumed to be zero.

See `reflectioncoef()` for more details.
"""
struct Rayleigh{T1,T2,T3} <: ReflectionModel
  ρᵣ::T1
  cᵣ::T2
  δ::T3
end

Rayleigh(ρᵣ, cᵣ) = Rayleigh(ρᵣ, cᵣ, 0.0)

# from APL-UW TR 9407 (1994), IV-6 Table 2
const Rock = Rayleigh(2.5, 2.5, 0.01374)
const Pebbles = Rayleigh(2.5, 1.8, 0.01374)
const SandyGravel = Rayleigh(2.492, 1.3370, 0.01705)
const CoarseSand = Rayleigh(2.231, 1.2503, 0.01638)
const MediumSand = Rayleigh(1.845, 1.1782, 0.01624)
const FineSand = Rayleigh(1.451, 1.1073, 0.01602)
const VeryFineSand = Rayleigh(1.268, 1.0568, 0.01875)
const ClayeySand = Rayleigh(1.224, 1.0364, 0.02019)
const CoarseSilt = Rayleigh(1.195, 1.0179, 0.02158)
const SandySilt = Rayleigh(1.169, 0.9999, 0.01261)
const Silt = Rayleigh(1.149, 0.9873, 0.00386)
const FineSilt = Rayleigh(1.148, 0.9861, 0.00306)
const SandyClay = Rayleigh(1.147, 0.9849, 0.00242)
const SiltyClay = Rayleigh(1.146, 0.9824, 0.00163)
const Clay = Rayleigh(1.145, 0.98, 0.00148)

reflectioncoef(rm::Rayleigh, f, θ) = reflectioncoef(θ, rm.ρᵣ, rm.cᵣ, rm.δ)

"""
$(TYPEDEF)
Reflection model for a water surface affected by wind.

---

    SurfaceLoss(windspeed)

Create a reflection model for a surface affected by wind. `windspeed` is
given in m/s.
"""
struct SurfaceLoss{T} <: ReflectionModel
  windspeed::T
end

reflectioncoef(rm::SurfaceLoss, f, θ) = complex(-surfaceloss(rm.windspeed, f, θ), 0.0)

const Vacuum = ReflectionCoef(complex(-1.0, 0.0))

# WMO sea states
# from APL-UW TR 9407 (1994), II-4 Table 2 (median windspeed)
const SeaState0 = SurfaceLoss(0.8)
const SeaState1 = SurfaceLoss(2.6)
const SeaState2 = SurfaceLoss(4.4)
const SeaState3 = SurfaceLoss(6.9)
const SeaState4 = SurfaceLoss(9.8)
const SeaState5 = SurfaceLoss(12.6)
const SeaState6 = SurfaceLoss(19.3)
const SeaState7 = SurfaceLoss(26.5)
const SeaState8 = SurfaceLoss(30.6)
const SeaState9 = SurfaceLoss(32.9)

### basic environmental model

Base.@kwdef struct BasicUnderwaterEnvironment{T1<:Altimetry, T2<:Bathymetry, T3<:SoundSpeedProfile, T4<:Number, T5<:ReflectionModel, T6<:ReflectionModel, T7} <: UnderwaterEnvironment
  altimetry::T1 = FlatSurface()
  bathymetry::T2 = ConstantDepth(20.0)
  ssp::T3 = IsoSSP(soundspeed())
  salinity::T4 = 35.0
  seasurface::T5 = SeaState1
  seabed::T6 = SandySilt
  noise::T7 = RedGaussianNoise(db2amp(120.0))
end

"""
    UnderwaterEnvironment(; altimetry, bathymetry, ssp, salinity, seasurface, seabed, noise)

Create an underwater environment.
"""
UnderwaterEnvironment(; kwargs...) = BasicUnderwaterEnvironment(; kwargs...)

altimetry(env::BasicUnderwaterEnvironment) = env.altimetry
bathymetry(env::BasicUnderwaterEnvironment) = env.bathymetry
ssp(env::BasicUnderwaterEnvironment) = env.ssp
salinity(env::BasicUnderwaterEnvironment) = env.salinity
seasurface(env::BasicUnderwaterEnvironment) = env.seasurface
seabed(env::BasicUnderwaterEnvironment) = env.seabed
noise(env::BasicUnderwaterEnvironment) = env.noise

### noise models

"""
$(TYPEDEF)
Ambient noise model with Gaussian noise with a `1/f²` power spectral density.

---

    RedGaussianNoise(σ)

Create an ambient noise model with variance `σ²` and `1/f²` power spectral density.
"""
struct RedGaussianNoise{T} <: NoiseModel
  σ::T
end

function record(noisemodel::RedGaussianNoise, duration, fs; start=0.0)
  analytic(signal(rand(RedGaussian(; n=round(Int, duration*fs), σ=noisemodel.σ)), fs))
end

### basic source & recevier models

"""
$(TYPEDEF)
Narrowband acoustic source.
"""
struct NarrowbandAcousticSource{T1,T2,T3,T4} <: AcousticSource
  pos::NTuple{3,T1}
  f::T2
  A::T3
  ϕ::T4
end

"""
$(SIGNATURES)
Create a narrowband acoustic source with frequency `f` Hz at location (`x`, `y`, `z`).
The `sourcelevel` is in µPa @ 1m. A phase `ϕ` may be optionlly specified.
"""
NarrowbandAcousticSource(x, y, z, f; sourcelevel=db2amp(180.0), ϕ=0.0) = NarrowbandAcousticSource(promote(x, y, z), f, sourcelevel, ϕ)

"""
$(SIGNATURES)
Create a narrowband acoustic source with frequency `f` Hz at location (`x`, z`).
The `sourcelevel` is in µPa @ 1m. A phase `ϕ` may be optionlly specified.
"""
NarrowbandAcousticSource(x, z, f; sourcelevel=db2amp(180.0), ϕ=0.0) = NarrowbandAcousticSource(promote(x, 0, z), f, sourcelevel, ϕ)

"""
$(SIGNATURES)
Create a narrowband acoustic source with frequency `f` Hz at location (`x`, `y`, `z`).
The `sourcelevel` is in µPa @ 1m. A phase `ϕ` may be optionlly specified.
"""
AcousticSource(x, y, z, f; sourcelevel=db2amp(180.0), ϕ=0.0) = NarrowbandAcousticSource(promote(x, y, z), f, sourcelevel, ϕ)

"""
$(SIGNATURES)
Create a narrowband acoustic source with frequency `f` Hz at location (`x`, `y`, `z`).
The `sourcelevel` is in µPa @ 1m. A phase `ϕ` may be optionlly specified.
"""
AcousticSource(x, z, f; sourcelevel=db2amp(180.0), ϕ=0.0) = NarrowbandAcousticSource(promote(x, 0, z), f, sourcelevel, ϕ)

location(tx::NarrowbandAcousticSource) = tx.pos
nominalfrequency(tx::NarrowbandAcousticSource) = tx.f
phasor(tx::NarrowbandAcousticSource) = tx.A * cis(tx.ϕ)
record(tx::NarrowbandAcousticSource, duration, fs; start=0.0) = signal(tx.A .* cis.(2π .* tx.f .* (start:1/fs:start+duration-1/fs) .+ tx.ϕ), fs)

"""
$(TYPEDEF)
Narrowband pulsed acoustic source.
"""
Base.@kwdef struct Pinger{T1,T2,T3,T4,T5,T6,T7,T8} <: AcousticSource
  pos::NTuple{3,T1}
  frequency::T2
  sourcelevel::T3 = db2amp(180.0)
  phase::T4 = 0.0
  duration::T5 = 0.02
  start::T6 = 0.0
  interval::T7 = 1.0
  window::T8 = nothing
end

"""
    Pinger(x, y, z, f; sourcelevel, phase, duration, start, interval, window)

Create a pulsed narrowband acoustic source with frequency `f` Hz at location (`x`, `y`, `z`).
Additional parameters that may be specified:
- `sourcelevel` in µPa @ 1m (default 180 dB)
- `phase` of the narrowband signal (default 0)
- `duration` of the pulse in seconds (default 20 ms)
- `start` time of one of the pulses in seconds (default 0)
- pulse repetition `interval` in seconds (default 1 second)
- `window` type (from `DSP.jl`) (default `nothing`)
"""
Pinger(x, y, z, f; kwargs...) = Pinger(; pos=promote(x, y, z), frequency=f, kwargs...)

"""
    Pinger(x, z, f; sourcelevel, phase, duration, start, interval, window)

Create a pulsed narrowband acoustic source with frequency `f` Hz at location (`x`, `z`).
Additional parameters that may be specified:
- `sourcelevel` in µPa @ 1m (default 180 dB)
- `phase` of the narrowband signal (default 0)
- `duration` of the pulse in seconds (default 20 ms)
- `start` time of one of the pulses in seconds (default 0)
- pulse repetition `interval` in seconds (default 1 second)
- `window` type (from `DSP.jl`) (default `nothing`)
"""
Pinger(x, z, f; kwargs...) = Pinger(; pos=promote(x, 0, z), frequency=f, kwargs...)

location(tx::Pinger) = tx.pos
nominalfrequency(tx::Pinger) = tx.frequency
phasor(tx::Pinger) = tx.sourcelevel * cis(tx.phase)

function record(pinger::Pinger, duration, fs; start=0.0)
  nsamples = round(Int, duration * fs)
  ping = cw(pinger.frequency, pinger.duration, fs; phase=pinger.phase, window=pinger.window)
  x = zeros(eltype(ping), nsamples)
  k1 = ceil(Int, (start - pinger.start - pinger.duration) / pinger.interval)
  k2 = floor(Int, (start - pinger.start + duration) / pinger.interval)
  n = length(ping)
  for k ∈ k1:k2
    ndx = round(Int, (pinger.start + k * pinger.interval - start) * fs) + 1
    if ndx > 0 && ndx + n ≤ nsamples
      x[ndx:ndx+n-1] .= ping
    elseif ndx > 0
      x[ndx:end] .= ping[1:nsamples-ndx+1]
    elseif ndx + n ≤ nsamples
      x[1:n+ndx-1] .= ping[2-ndx:end]
    else
      x .= ping[2-ndx:1-ndx+nsamples]
    end
  end
  signal(pinger.sourcelevel * x, fs)
end

"""
$(TYPEDEF)
Omnidirectional acoustic receiver.
"""
struct BasicAcousticReceiver{T} <: AcousticReceiver
  pos::NTuple{3,T}
end

"""
$(SIGNATURES)
Create an omnidirectional acoustic receiver at location (`x`, `y`, `z`).
"""
BasicAcousticReceiver(x, y, z) = BasicAcousticReceiver(promote(x, y, z))

"""
$(SIGNATURES)
Create an omnidirectional acoustic receiver at location (`x`, `z`).
"""
BasicAcousticReceiver(x, z) = BasicAcousticReceiver(promote(x, 0, z))

"""
$(SIGNATURES)
Create an omnidirectional acoustic receiver at location (`x`, `y`, `z`).
"""
AcousticReceiver(x, y, z) = BasicAcousticReceiver(promote(x, y, z))

"""
$(SIGNATURES)
Create an omnidirectional acoustic receiver at location (`x`, `z`).
"""
AcousticReceiver(x, z) = BasicAcousticReceiver(promote(x, 0, z))

location(rx::BasicAcousticReceiver) = rx.pos

### receiver grids for transmission loss computation

"""
$(TYPEDEF)
A 2D Cartesian grid of omnidirectional acoustic receivers.
"""
struct AcousticReceiverGrid2D{T} <: AbstractArray{BasicAcousticReceiver{T},2}
  xrange::StepRangeLen{T,T,T}
  zrange::StepRangeLen{T,T,T}
end

"""
$(SIGNATURES)
Create a 2D Cartesian grid of omnidirectional acoustic receivers with `nx` × `nz`
receviers starting (`xmin`, `zmin`) with step sizes `xstep` and `zstep`.
"""
function AcousticReceiverGrid2D(xmin, xstep, nx, zmin, zstep, nz)
  xmin, xstep, zmin, zstep = promote(xmin, xstep, zmin, zstep)
  AcousticReceiverGrid2D(StepRangeLen(xmin, xstep, nx), StepRangeLen(zmin, zstep, nz))
end

Base.size(g::AcousticReceiverGrid2D) = (g.xrange.len, g.zrange.len)
Base.getindex(g::AcousticReceiverGrid2D, I::Vararg{Int,2}) = AcousticReceiver(g.xrange[I[1]], g.zrange[I[2]])
Base.setindex!(g::AcousticReceiverGrid2D, v, I::Vararg{Int,2}) = throw(ArgumentError("AcousticReceiverGrid2D is readonly"))

"""
$(TYPEDEF)
A 3D Cartesian grid of omnidirectional acoustic receivers.
"""
struct AcousticReceiverGrid3D{T} <: AbstractArray{BasicAcousticReceiver{T},3}
  xrange::StepRangeLen{T,T,T}
  yrange::StepRangeLen{T,T,T}
  zrange::StepRangeLen{T,T,T}
end

"""
$(SIGNATURES)
Create a 3D Cartesian grid of omnidirectional acoustic receivers with `nx` × `ny` × `nz`
receviers starting (`xmin`, `ymin`, `zmin`) with step sizes `xstep`, `ystep`, and `zstep`.
"""
function AcousticReceiverGrid3D(xmin, xstep, nx, ymin, ystep, ny, zmin, zstep, nz)
  xmin, xstep, ymin, ystep, zmin, zstep = promote(xmin, xstep, ymin, ystep, zmin, zstep)
  AcousticReceiverGrid3D(StepRangeLen(xmin, xstep, nx), StepRangeLen(ymin, ystep, ny), StepRangeLen(zmin, zstep, nz))
end

Base.size(g::AcousticReceiverGrid3D) = (g.xrange.len, g.yrange.len, g.zrange.len)
Base.getindex(g::AcousticReceiverGrid3D, I::Vararg{Int,3}) = AcousticReceiver(g.xrange[I[1]], g.yrange[I[2]], g.zrange[I[3]])
Base.setindex!(g::AcousticReceiverGrid3D, v, I::Vararg{Int,3}) = throw(ArgumentError("AcousticReceiverGrid3D is readonly"))
