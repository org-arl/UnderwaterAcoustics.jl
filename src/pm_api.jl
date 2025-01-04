export UnderwaterEnvironment, PekerisWaveguide
export FluidBoundary, RigidBoundary, PressureReleaseBoundary
export transmission_loss, acoustic_field, arrivals, impulse_response, channel, transmit
export NarrowbandAcousticSource, AcousticReceiver, AcousticReceiverGrid2D, AcousticReceiverGrid3D
export spl, frequency

###############################################################################
### propagation model API

"""
Superclass for all propagation models.
"""
abstract type AbstractPropagationModel end

"""
Superclass for all acoustic arrivals.
"""
abstract type AcousticArrival end

"""
    transmission_loss(pm, tx, rxs)

Compute the transmission loss from the source `tx` to the receivers `rxs`
using propagation model `pm`. If `rxs` denotes a single receiver, the result is a
scalar. If `rxs` is an `AbstractArray`, the result is an array of transmission
losses (in dB) with the same shape as `rxs`.
"""
function transmission_loss(pm, tx, rxs; kwargs...)
  20 * log10.(abs.(acoustic_field(pm, tx, rxs; kwargs...))) .- spl(tx)
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
    FluidBoundary(c, ρ, δ)

Create a fluid half-space boundary with sound speed `c`, density `ρ`, and
dimensionless absorption coefficient `δ`.
"""
struct FluidBoundary{T1,T2,T3} <: AbstractAcousticBoundary
  c::T1
  ρ::T2
  δ::T3
end

FluidBoundary(c, ρ) = FluidBoundary(c, ρ, 0.0)

function Base.show(io::IO, b::FluidBoundary)
  if b.c == Inf
    print(io, "RigidBoundary")
  elseif b.c == 0
    print(io, "PressureReleaseBoundary")
  elseif b.δ == 0
    print(io, "FluidBoundary(c=$(b.c), ρ=$(b.ρ))")
  else
    print(io, "FluidBoundary(c=$(b.c), ρ=$(b.ρ), δ=$(b.δ))")
  end
end

"""
Rigid boundary condition.
"""
const RigidBoundary = FluidBoundary(Inf, 0, 0)

"""
Pressure-release boundary condition.
"""
const PressureReleaseBoundary = FluidBoundary(0, 0, 0)

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
    density = 1025.0,
    seabed = RigidBoundary,
    surface = PressureReleaseBoundary,
    noise = nothing
  )
    bathymetry isa Number && (bathymetry = in_units(u"m", bathymetry))
    altimetry isa Number && (altimetry = in_units(u"m", altimetry))
    temperature isa Number && (temperature = in_units(u"°C", temperature))
    salinity isa Number && (salinity = in_units(u"ppt", salinity))
    density isa Number && (density = in_units(u"kg/m^3", density))
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
    PekerisWaveguide(; kwargs...)

Create a Pekeris waveguide environment with the given parameters. All parameters
are optional (have default values). Default values are specified below:

- `D`  = water depth (m)
- `c₁` = sound speed in water (m/s, computed from temperature and salinity)
- `ρ₁` = density of water (`1025` kg/m³)
- `c₂` = sound speed in seabed (`Inf` m/s)
- `ρ₂` = density of seabed (`2000` kg/m³)
- `δ₂` = dimensionless absorption coefficient for seabed (`0`)
- `σ₂` = seabed roughness (`0` m)
- `cₛ` = sound speed for surface (`0` m/s)
- `ρₛ` = density of surface (`0` kg/m³)
- `δₛ` = dimensionless absorption coefficient for surface (`0`)
- `σₛ` = surface roughness (`0` m)

Other parameters for `UnderwaterEnvironment` may be specified as well. For
example, `temperature`, `salinity`, `density` and `noise` may be specified.

A sound speed of `0` denotes a pressure release boundary. A sound speed of `Inf`
denotes a rigid boundary.

Returns an underwater environment with the specified parameters to ensure it
is a Pekeris waveguide.
"""
function PekerisWaveguide(;
  D=100.0, c₁=1500.0, ρ₁=1025.0, c₂=Inf, ρ₂=2000.0, δ₂=0.0, σ₂=0.0, cₛ=0.0,
  ρₛ=0.0, δₛ=0.0, σₛ=0.0, kwargs...
)
  UnderwaterEnvironment(
    bathymetry = D,
    altimetry = 0.0,
    soundspeed = c₁,
    density = ρ₁,
    seabed = FluidBoundary(c₂, ρ₂, δ₂),
    surface = FluidBoundary(cₛ, ρₛ, δₛ),
    noise = nothing,
    kwargs...
  )
end

"""
    is_range_dependent(env)

Return `true` if any quantity (e.g. sound speed, bathymetry, etc) in the
environment `env` depends on the horizontal position, and `false` otherwise.
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
  p = Tuple(position(tx))
  print(io, "TX$p")
end

function Base.show(io::IO, tx::AbstractAcousticReceiver)
  p = Tuple(position(tx))
  print(io, "RX$p")
end

"""
    NarrowbandAcousticSource(pos, frequency, spl=0)

Narrowband source at position `pos` with specified `frequency` and source level
`spl` (dB re 1 µPa @ 1 m).

If the position of the source is unknown, it may be specified as `nothing`. This
is useful when the propagation model does not require the source position (e.g.,
data-driven models).
"""
struct NarrowbandAcousticSource{T1,T2,T3} <: AbstractAcousticSource
  pos::T1
  frequency::T2
  spl::T3
  function NarrowbandAcousticSource(pos, frequency, spl=0)
    p = Position(pos)
    f = in_units(u"Hz", frequency)
    s = in_units(u"dB", spl)
    new{typeof(p),typeof(f),typeof(s)}(p, f, s)
  end
end

"""
    AcousticReceiver(pos)
    AcousticReceiver(x, z)
    AcousticReceiver(x, y, z)

Receiver at position `pos`.
"""
struct AcousticReceiver{T1} <: AbstractAcousticReceiver
  pos::T1
  function AcousticReceiver(pos)
    p = Position(pos)
    new{typeof(p)}(p)
  end
end

AcousticReceiver(x, z) = AcousticReceiver(Position(x, z))
AcousticReceiver(x, y, z) = AcousticReceiver(Position(x, y, z))

"""
    AcousticReceiverGrid2D(xrange, zrange)

A 2D grid of receivers with the specified position ranges.
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

A 3D grid of receivers with the specified position ranges.
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
    position(tx::AbstractAcousticSource)
    position(rx::AbstractAcousticReceiver)

Get the position of the source or receiver.
"""
Base.position(tx::NarrowbandAcousticSource) = tx.pos
Base.position(rx::AcousticReceiver) = rx.pos

"""
    frequency(tx::AbstractAcousticSource)

Get the nominal frequency of an acoustic source.
"""
frequency(tx::NarrowbandAcousticSource) = tx.frequency

"""
    spl(tx::AbstractAcousticSource)

Get the source level of an acoustic source.
"""
spl(tx::NarrowbandAcousticSource) = tx.spl
