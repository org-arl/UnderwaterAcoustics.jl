import Random: default_rng
import SignalAnalysis: amp2db
import Printf: @printf

export UnderwaterEnvironment, transmission_loss, acoustic_field, arrivals
export impulse_response, channel, transmit, spl, frequency, location
export AcousticSource, AcousticReceiver, AcousticReceiverGrid2D, AcousticReceiverGrid3D

###############################################################################
### propagation model API

"""
Superclass for all propagation models.
"""
abstract type AbstractPropagationModel end

"""
Superclass for all acoustic arrivals.
"""
abstract type AbstractAcousticArrival end

"""
Type representing a single acoustic ray arrival.

Fields:
- `t`: arrival time (s)
- `ϕ`: complex amplitude
- `ns`: number of surface bounces
- `nb`: number of bottom bounces
- `θₛ`: launch angle at source (rad)
- `θᵣ`: arrival angle at receiver (rad)
- `path`: ray path (optional, vector of 3-tuples or missing)
"""
struct RayArrival{T1,T2,T3,T4,T5} <: AbstractAcousticArrival
  t::T1                 # arrival time
  ϕ::Complex{T2}        # complex amplitude
  ns::Int               # number of surface bounces
  nb::Int               # number of bottom bounces
  θₛ::T3                # launch angle at source
  θᵣ::T4                # arrival angle at receiver
  path::T5              # ray path
end

function Base.show(io::IO, a::RayArrival)
  @printf(io, "∠%5.1f° %2d⤒ %2d⤓ ∠%5.1f° | %6.2f ms | %5.1f dB ϕ%6.1f°%s",
    rad2deg(a.θₛ), a.ns, a.nb, rad2deg(a.θᵣ), 1000*a.t, amp2db(abs(a.ϕ)),
    rad2deg(angle(a.ϕ)), a.path === missing || length(a.path) == 0 ? "" : "↝")
end

"""
Type representing a single acoustic normal mode arrival.

Fields:
- `t`: arrival time (s)
- `m`: mode number
- `kz`: vertical wavenumber (rad/m)
- `kr`: horizontal wavenumber (rad/m)
- `θ`: propagation angle (rad)
- `ϕ`: complex multiplier
"""
struct ModeArrival{T1,T2,T3,T4,T5} <: AbstractAcousticArrival
  t::T1                 # arrival time
  m::Int                # mode number
  kz::T2                # vertical wavenumber
  kr::T3                # horizontal wavenumber
  θ::T4                 # propagation angle
  ϕ::Complex{T5}        # complex multiplier
end

function Base.show(io::IO, a::ModeArrival)
  @printf(io, "%3d | %6.3f→ %6.3f↑ ∠%5.1f° | %6.2f ms | %5.1f dB ϕ%6.1f°",
    a.m, a.kr, a.kz, rad2deg(a.θ), 1000*a.t, amp2db(abs(a.ϕ)), rad2deg(angle(a.ϕ)))
end

"""
    transmission_loss(pm, tx, rxs)

Compute the transmission loss from the source `tx` to the receivers `rxs`
using propagation model `pm`. If `rxs` denotes a single receiver, the result is a
scalar. If `rxs` is an `AbstractArray`, the result is an array of transmission
losses (in dB) with the same shape as `rxs`.
"""
function transmission_loss(pm, tx, rxs; kwargs...)
  spl(tx) .- 20 * log10.(abs.(acoustic_field(pm, tx, rxs; kwargs...)))
end

"""
    acoustic_field(pm, txs, rxs)

Compute the acoustic field at the receivers `rxs` due to the sources `txs`
using propagation model `pm`. If `rxs` denotes a single receiver, the result is a
complex scalar. If `rxs` is an `AbstractArray`, the result is an array of
complex numbers with the same shape as `rxs`. The amplitude of the field is
related to the transmission loss, and the angle is related to the acoustic phase
at the source frequency.
"""
function acoustic_field end

"""
    arrivals(pm, tx, rx)

Compute the arrivals at the receiver `rx` due to the source `tx` using
propagation model `pm`. Returns an array of arrivals.
"""
function arrivals end

"""
    impulse_response(pm, tx, rx, fs; abstime=false)

Compute the impulse response at the receiver `rx` due to the source `tx` using
propagation model `pm` at the given sampling frequency `fs`. If `abstime` is
`true`, the result is in absolute time from the start of transmission. Otherwise,
the result is relative to the earliest arrival time of the signal at the receiver.
"""
function impulse_response end

###############################################################################
### channel model API

"""
Superclass for all channel models.
"""
abstract type AbstractChannelModel end

"""
    channel(pm, txs, rxs; abstime=false)

Compute a channel model from the sources `txs` to the receivers `rxs` using
propagation model `pm`. The result is a channel model with the same number of
input channels as the number of sources and output channels as the number of
receivers.

If `abstime` is `true`, the resulting channel model returns signals with time
relative to the start of transmission. Otherwise, the result is relative to the
earliest arrival time of the signal at any receiver.
"""
function channel end

"""
    transmit(ch, x; fs=framerate(x), txs=:, rxs=:)

Simulate the transmission of the signal `x`, sampled at rate `fs`, through the
channel model `ch`. If `txs` is specified, it specifies the indices of the
sources active in the simulation. The number of sources must match the number
of channels in the input signal. If `rxs` is specified, it specifies the
indices of the receivers active in the simulation. Returns the received signal
at the specified (or all) receivers.
"""
function transmit end

###############################################################################
### boundary conditions

"""
Superclass for all boundary conditions.
"""
abstract type AbstractAcousticBoundary end

"""
    reflection_coef(bc::AbstractAcousticBoundary, frequency, θ)
    reflection_coef(bc::AbstractAcousticBoundary, frequency, θ, ρ, c)

Compute the complex reflection coefficient at a fluid-fluid boundary of type
`bc` at incidence angle `θ` and `frequency`. The density and sound speed in
the water are given by `ρ` and `c`, respectively.
"""
function reflection_coef(bc::AbstractAcousticBoundary, frequency, θ)
  reflection_coef(bc::AbstractAcousticBoundary, frequency, θ, water_density(), soundspeed())
end

###############################################################################
### environment model API

"""
    UnderwaterEnvironment(; kwargs...)

Create a generic underwater environment with the given parameters. The following
parameters are supported:

- `bathymetry` = bathymetry model
- `altimetry` = altimetry model
- `temperature` = temperature model
- `salinity` = salinity model
- `soundspeed` = sound speed profile model
- `density` = density model
- `seabed` = seabed sediment model
- `surface` = surface model
- `noise` = noise model

All parameters are optional and have default values. Parameters may be modified
after construction by setting the corresponding property in the environment.
"""
mutable struct UnderwaterEnvironment
  bathymetry::Any
  altimetry::Any
  temperature::Any
  salinity::Any
  soundspeed::Any
  density::Any
  seabed::Any
  surface::Any
  noise::Any
  function UnderwaterEnvironment(;
    bathymetry = 100.0,
    altimetry = 0.0,
    temperature = 25.0,
    salinity = 35.0,
    soundspeed = nothing,
    density = nothing,
    seabed = RigidBoundary,
    surface = PressureReleaseBoundary,
    noise = nothing
  )
    bathymetry isa Number && (bathymetry = in_units(u"m", bathymetry))
    altimetry isa Number && (altimetry = in_units(u"m", altimetry))
    temperature isa Number && (temperature = in_units(u"°C", temperature))
    salinity isa Number && (salinity = in_units(u"ppt", salinity))
    density isa Number && (density = in_units(u"kg/m^3", density))
    density = something(density, water_density(temperature, salinity))
    soundspeed isa Number && (soundspeed = in_units(u"m/s", soundspeed))
    soundspeed = something(soundspeed, UnderwaterAcoustics.soundspeed(temperature, salinity))
    new(bathymetry, altimetry, temperature, salinity, soundspeed, density, seabed, surface, noise)
  end
end

function Base.show(io::IO, env::UnderwaterEnvironment)
  println(io, "UnderwaterEnvironment(")
  println(io, "  bathymetry = $(env.bathymetry), ")
  println(io, "  altimetry = $(env.altimetry), ")
  println(io, "  temperature = $(env.temperature), ")
  println(io, "  salinity = $(env.salinity), ")
  println(io, "  soundspeed = $(env.soundspeed), ")
  println(io, "  density = $(env.density), ")
  println(io, "  seabed = $(env.seabed), ")
  println(io, "  surface = $(env.surface), ")
  println(io, "  noise = $(env.noise)")
  println(io, ")")
end

"""
    is_range_dependent(env)

Return `true` if any quantity (e.g. sound speed, bathymetry, etc) in the
environment `env` depends on the horizontal location, and `false` otherwise.
"""
function is_range_dependent(env::UnderwaterEnvironment)
  for p ∈ (:bathymetry, :altimetry, :temperature, :salinity, :soundspeed, :density, :seabed, :surface, :noise)
    is_range_dependent(getproperty(env, p)) && return true
  end
  false
end

"""
    isospeed(env)

Return `true` if the sound speed in the environment `env` is a constant.
"""
isospeed(env::UnderwaterEnvironment) = is_constant(env.soundspeed)

###############################################################################
### sources and receivers API

"""
Superclass for all acoustic sources.
"""
abstract type AbstractAcousticSource end

"""
Superclass for all acoustic receivers.
"""
abstract type AbstractAcousticReceiver end

function Base.show(io::IO, tx::AbstractAcousticSource)
  p = Tuple(location(tx))
  print(io, "TX$p")
end

function Base.show(io::IO, tx::AbstractAcousticReceiver)
  p = Tuple(location(tx))
  print(io, "RX$p")
end

"""
    AcousticSource(pos, frequency, spl=0)

An source at location `pos` with nominal `frequency` and source level `spl`
(dB re 1 µPa @ 1 m). The source is assumed to be omnidirectional and well
approximated by a point source. While the source may have some bandwidth,
the nominal frequency is used to estimate propagation effects such as
absorption, reflection coefficients, etc.

If the location of the source is unknown, it may be specified as `nothing`. This
is useful when the propagation model does not require the source location (e.g.,
data-driven models).
"""
struct AcousticSource{T1,T2,T3} <: AbstractAcousticSource
  pos::T1
  frequency::T2
  spl::T3
  function AcousticSource(pos, frequency, spl=0)
    p = XYZ(pos)
    f = in_units(u"Hz", frequency)
    s = in_units(u"dB", spl)
    new{typeof(p),typeof(f),typeof(s)}(p, f, s)
  end
end

"""
    AcousticReceiver(pos)
    AcousticReceiver(x, z)
    AcousticReceiver(x, y, z)

Receiver at location `pos`.
"""
struct AcousticReceiver{T1} <: AbstractAcousticReceiver
  pos::T1
  function AcousticReceiver(pos)
    p = XYZ(pos)
    new{typeof(p)}(p)
  end
end

AcousticReceiver(x, z) = AcousticReceiver(XYZ(x, z))
AcousticReceiver(x, y, z) = AcousticReceiver(XYZ(x, y, z))

"""
    AcousticReceiverGrid2D(xrange, zrange)

A 2D grid of receivers with the specified location ranges.
"""
struct AcousticReceiverGrid2D{T1,T2} <: AbstractArray{AcousticReceiver,2}
  xrange::T1
  zrange::T2
  function AcousticReceiverGrid2D(xrange, zrange)
    xrange isa Number && (xrange = StepRangeLen(xrange, 0, 1))
    zrange isa Number && (zrange = StepRangeLen(zrange, 0, 1))
    xrange = StepRangeLen(xrange)
    zrange = StepRangeLen(zrange)
    new{typeof(xrange),typeof(zrange)}(xrange, zrange)
  end
end

Base.size(g::AcousticReceiverGrid2D) = (g.xrange.len, g.zrange.len)
Base.getindex(g::AcousticReceiverGrid2D, I::Vararg{Int,2}) = AcousticReceiver(g.xrange[I[1]], g.zrange[I[2]])
Base.setindex!(g::AcousticReceiverGrid2D, v, I::Vararg{Int,2}) = error("AcousticReceiverGrid2D is readonly")

"""
    AcousticReceiverGrid3D(xrange, yrange, zrange)

A 3D grid of receivers with the specified location ranges.
"""
struct AcousticReceiverGrid3D{T1,T2,T3} <: AbstractArray{AcousticReceiver,3}
  xrange::T1
  yrange::T2
  zrange::T3
  function AcousticReceiverGrid3D(xrange, yrange, zrange)
    xrange isa Number && (xrange = StepRangeLen(xrange, 0, 1))
    yrange isa Number && (yrange = StepRangeLen(yrange, 0, 1))
    zrange isa Number && (zrange = StepRangeLen(zrange, 0, 1))
    xrange = StepRangeLen(xrange)
    yrange = StepRangeLen(yrange)
    zrange = StepRangeLen(zrange)
    new{typeof(xrange),typeof(yrange),typeof(zrange)}(xrange, yrange, zrange)
  end
end

Base.size(g::AcousticReceiverGrid3D) = (g.xrange.len, g.yrange.len, g.zrange.len)
Base.getindex(g::AcousticReceiverGrid3D, I::Vararg{Int,3}) = AcousticReceiver(g.xrange[I[1]], g.yrange[I[2]], g.zrange[I[3]])
Base.setindex!(g::AcousticReceiverGrid3D, v, I::Vararg{Int,3}) = error("AcousticReceiverGrid3D is readonly")

"""
    location(tx::AbstractAcousticSource)
    location(rx::AbstractAcousticReceiver)

Get the location of the source or receiver.
"""
location(tx::AcousticSource) = tx.pos
location(rx::AcousticReceiver) = rx.pos

"""
    frequency(tx::AbstractAcousticSource)

Get the nominal frequency of an acoustic source.
"""
frequency(tx::AcousticSource) = tx.frequency

"""
    spl(tx::AbstractAcousticSource)

Get the source level of an acoustic source.
"""
spl(tx::AcousticSource) = tx.spl

###############################################################################
### noise modeling API

"""
Superclass for all noise models.
"""
abstract type AbstractNoiseModel end

"""
    rand([rng::AbstractRNG, ] noise::AbstractNoiseModel, nsamples, fs)

Generate random noise samples from the noise model `noise` with the specified
size `nsamples` (2-tuple for multichannel noise). The noise is returned as a
signal sampled at `fs`. The optional `rng` argument specifies the random number
generator to use.
"""
function Base.rand(noise::AbstractNoiseModel, nsamples, fs)
  Base.rand(default_rng(), noise, nsamples, fs)
end
