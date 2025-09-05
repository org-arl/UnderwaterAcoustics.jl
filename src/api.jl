import Random: default_rng, AbstractRNG
import Logging: NullLogger, with_logger
import InteractiveUtils: subtypes
import SignalAnalysis: amp2db, SampledSignal, samples, signal, framerate, db2amp
import SignalAnalysis: nchannels, isanalytic, analytic, filt
import Printf: @printf

export UnderwaterEnvironment, transmission_loss, acoustic_field, arrivals, models
export impulse_response, channel, transmit, spl, frequency, location, samples, is_isovelocity
export AcousticSource, AcousticReceiver, AcousticReceiverGrid2D, AcousticReceiverGrid3D

###############################################################################
### propagation model API

"""
Superclass for all propagation models.
"""
abstract type AbstractPropagationModel end

"""
Superclass for all ray propagation models.
"""
abstract type AbstractRayPropagationModel <: AbstractPropagationModel end

"""
Superclass for all mode propagation models.
"""
abstract type AbstractModePropagationModel <: AbstractPropagationModel end

"""
    models()
    models(mtype::Type{<:AbstractPropagationModel})

Return a list of all available propagation models. If `mtype` is specified,
return only models of that type.
"""
function models(mtype::Type{<:AbstractPropagationModel}=AbstractPropagationModel)::Vector{Type{<:AbstractPropagationModel}}
  isabstracttype(mtype) || return Type{<:AbstractPropagationModel}[mtype]
  mapreduce(models, vcat, subtypes(mtype); init=Type{<:AbstractPropagationModel}[])
end

"""
Superclass for all acoustic arrivals.
"""
abstract type AbstractAcousticArrival end

"""
Type representing a single acoustic ray arrival.

Properties:
- `t` / `time`: arrival time (s)
- `ϕ` / `phasor`: complex amplitude
- `ns` / `surface_bounces`: number of surface bounces
- `nb` / `bottom_bounces`: number of bottom bounces
- `θₛ` / `launch_angle`: launch angle at source (rad)
- `θᵣ` / `arrival_angle`: arrival angle at receiver (rad)
- `path`: ray path (optional, vector of 3-tuples or missing)

The properties are accessible with the short names for brevity, and longer more
descriptive names where readability is desired or unicode symbols are undesired.
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

function Base.show(io::IO, ::MIME"text/plain", v::AbstractVector{<:RayArrival})
  println(io, "$(length(v)) element Vector{RayArrival}:")
  for a ∈ v
    println(io, a)
  end
end

# alternative property names for readability
Base.propertynames(a::RayArrival) = (:t, :ϕ, :ns, :nb, :θₛ, :θᵣ, :path,
  :time, :phasor, :surface_bounces, :bottom_bounces, :launch_angle, :arrival_angle)

function Base.getproperty(a::RayArrival, sym::Symbol)
  sym === :time && return getfield(a, :t)
  sym === :phasor && return getfield(a, :ϕ)
  sym === :surface_bounces && return getfield(a, :ns)
  sym === :bottom_bounces && return getfield(a, :nb)
  sym === :launch_angle && return getfield(a, :θₛ)
  sym === :arrival_angle && return getfield(a, :θᵣ)
  getfield(a, sym)
end

"""
Type representing a single acoustic mode arrival.

Properties:
- `m` / `mode`: mode number
- `kᵣ` / `hwavenumber`: horizontal wavenumber (rad/m)
- `ψ(z)` / `mode_function`: mode function
- `v` / `group_velocity`: group velocity (m/s)
- `vₚ` / `phase_velocity`: phase velocity (m/s)

The properties are accessible with the short names for brevity, and longer more
descriptive names where readability is desired or unicode symbols are undesired.
"""
struct ModeArrival{T1,T2,T3,T4} <: AbstractAcousticArrival
  m::Int                # mode number
  kᵣ::T1                # horizontal wavenumber
  ψ::T2                 # mode function
  v::T3                 # group velocity
  vₚ::T4                # phase velocity
end

function Base.show(io::IO, a::ModeArrival)
  @printf(io, "%8s: kᵣ = %s rad/m", "mode $(a.m)", string(round(ComplexF64(a.kᵣ); digits=6)))
  a.v isa Real && @printf(io, ", v = %0.2f m/s", a.v)
  a.vₚ isa Real && @printf(io, ", vₚ = %0.2f m/s", a.vₚ)
end

function Base.show(io::IO, ::MIME"text/plain", v::AbstractVector{<:ModeArrival})
  println(io, "$(length(v)) element Vector{ModeArrival}:")
  for a ∈ v
    println(io, a)
  end
end

# alternative property names for readability
Base.propertynames(a::ModeArrival) = (:m, :kᵣ, :ψ, :v, :vₚ,
  :mode, :hwavenumber, :mode_function, :group_velocity, :phase_velocity)

function Base.getproperty(a::ModeArrival, sym::Symbol)
  sym === :mode && return getfield(a, :m)
  sym === :hwavenumber && return getfield(a, :kᵣ)
  sym === :mode_function && return getfield(a, :ψ)
  sym === :group_velocity && return getfield(a, :v)
  sym === :phase_velocity && return getfield(a, :vₚ)
  getfield(a, sym)
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
    acoustic_field(pm, tx, rxs)

Compute the acoustic field at the receivers `rxs` due to the source `tx`
using propagation model `pm`. If `rxs` denotes a single receiver, the result is a
complex scalar. If `rxs` is an `AbstractArray`, the result is an array of
complex numbers with the same shape as `rxs`. The amplitude of the field is
related to the transmission loss, and the angle is related to the acoustic phase
at the source frequency.
"""
function acoustic_field(pm, tx, rxs::AbstractArray; kwargs...)
  tmap(rx -> acoustic_field(pm, tx, rx; kwargs...), rxs)
end

"""
    arrivals(pm, tx, rx; paths=true)

Compute the arrivals at the receiver `rx` due to the source `tx` using
propagation model `pm`. Returns an array of arrivals.

For ray models, eigenray paths are typically included in the arrivals.
However, if they are not needed, one may set `paths=false` to allow the
propagation model to avoid computing them.
"""
function arrivals end

"""
    impulse_response(pm, tx, rx, fs; abstime=false, ntaps=nothing)

Compute the impulse response at the receiver `rx` due to the source `tx` using
propagation model `pm` at the given sampling frequency `fs`. If `abstime` is
`true`, the result is in absolute time from the start of transmission. Otherwise,
the result is relative to the earliest arrival time of the signal at the receiver
(possibly with some guard period to accommodate acausal response). `ntaps`
specifies the number of taps in the impulse response. If not specified, the
number of taps is chosen automatically based on the arrival times.
"""
function impulse_response end

###############################################################################
### channel model API

"""
Superclass for all channel models.
"""
abstract type AbstractChannelModel end

struct SampledPassbandChannel{T1,T2} <: AbstractChannelModel
  irs::T1       # impulse responses
  noise::T2     # noise model
  fs::Float64   # sampling rate
  t0::Int       # start time
end

function Base.show(io::IO, ch::SampledPassbandChannel)
  print(io, "SampledPassbandChannel($(size(ch.irs,1))×$(size(ch.irs,2)), $(ch.fs) Hz)")
end

"""
    channel(pm, txs, rxs, fs; noise=nothing, kwargs...)

Compute a channel model from the sources `txs` to the receivers `rxs` using
propagation model `pm`. The result is a channel model with the same number of
input channels as the number of sources and output channels as the number of
receivers. The channel model accepts signals sampled at rate `fs` and returns
signals sampled at the same rate.

An additive noise model may be optionally specified as `noise`. If specified,
it is used to corrupt the received signals.

Propagation model specific keyword arguments `kwargs` supported by
`impulse_response()` can be passed through when generating a channel.
"""
function channel(pm, txs, rxs, fs; noise=nothing, kwargs...)
  txs isa AbstractArray || (txs = [txs])
  rxs isa AbstractArray || (rxs = [rxs])
  ch = [
    samples(impulse_response(pm, tx, rx, fs; abstime=true, kwargs...)) .* db2amp(spl(tx))
    for tx in vec(txs), rx in vec(rxs)
  ]
  t0 = minimum(findfirst(x -> abs(x) > 0, ch1) for ch1 ∈ ch)
  t0 === nothing && error("No arrivals found")
  ch = map(ch1 -> @view(ch1[t0:end]), ch)
  len = maximum(length.(ch))
  ch = map(ch1 -> vcat(ch1, zeros(eltype(ch1), len - length(ch1))), ch)
  SampledPassbandChannel(ch, noise, Float64(fs), t0)
end

"""
    transmit(ch, x; txs=:, rxs=:, abstime=false, noisy=true, fs=nothing)

Simulate the transmission of passband signal `x` through the channel model `ch`.
If `txs` is specified, it specifies the indices of the sources active in the
simulation. The number of sources must match the number of channels in the
input signal. If `rxs` is specified, it specifies the indices of the
receivers active in the simulation. Returns the received signal at the
specified (or all) receivers.

`fs` specifies the sampling rate of the input signal. The output signal is
sampled at the same rate. If `fs` is not specified but `x` is a `SampledSignal`,
the sampling rate of `x` is used. Otherwise, the signal is assumed to be
sampled at the channel's sampling rate.

If `abstime` is `true`, the returned signals begin at the start of transmission.
Otherwise, the result is relative to the earliest arrival time of the signal
at any receiver. If `noisy` is `true` and the channel has a noise model
associated with it, the received signal is corrupted by additive noise.
"""
function transmit(ch::SampledPassbandChannel, x; txs=:, rxs=:, abstime=false, noisy=true, fs=nothing)
  # defaults
  N, M = size(ch.irs)
  txs = txs === (:) ? (1:N) : txs
  rxs = rxs === (:) ? (1:M) : rxs
  # validate inputs
  fs = something(fs, x isa SampledSignal ? framerate(x) : ch.fs)
  fs != ch.fs && error("Mismatched sampling rate (expected $(ch.fs) Hz, actual $(fs) Hz)")
  all(tx ∈ 1:N for tx ∈ txs) || error("Invalid transmitter indices ($txs ⊄ 1:$N)")
  all(rx ∈ 1:M for rx ∈ rxs) || error("Invalid receiver indices ($rxs ⊄ 1:$M)")
  nchannels(x) == length(txs) || error("Mismatched number of sources (expected $(length(txs)), actual $(nchannels(x)))")
  # simulate transmission
  t0 = abstime ? ch.t0 : 1
  flen = size(ch.irs[1], 1)
  x̄ = analytic(x)
  x̄ = vcat(x̄, zeros(eltype(x̄), flen - 1, size(x̄,2)))
  ȳ = zeros(Base.promote_eltype(x̄, ch.irs[1]), size(x̄,1) + t0 - 1, length(rxs))
  let x̄ = x̄     # avoids boxing of x̄ and resulting type instability
    Threads.@threads for i ∈ eachindex(rxs)
      for j ∈ eachindex(txs)
        ȳ[t0:end, i] .+= filt(ch.irs[txs[j], rxs[i]], x̄[:,j])
      end
    end
  end
  # add noise
  if isanalytic(x)
    noisy && ch.noise !== nothing && (ȳ .+= analytic(rand(ch.noise, size(ȳ); fs)))
    signal(ȳ, fs)
  else
    y = real(ȳ)
    noisy && ch.noise !== nothing && (y .+= rand(ch.noise, size(y); fs))
    signal(y, fs)
  end
end

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
- `pH` = pH model
- `soundspeed` = sound speed profile model
- `density` = density model
- `seabed` = seabed sediment model
- `surface` = surface model

All parameters are optional and have default values.
"""
struct UnderwaterEnvironment{T1,T2,T3,T4,T5,T6,T7,T8,T9}
  bathymetry::T1
  altimetry::T2
  temperature::T3
  salinity::T4
  pH::T5
  soundspeed::T6
  density::T7
  seabed::T8
  surface::T9
end

function UnderwaterEnvironment(;
  bathymetry = 100.0,
  altimetry = 0.0,
  temperature = 27.0,
  salinity = 35.0,
  pH = 8.1,
  soundspeed = nothing,
  density = nothing,
  seabed = RigidBoundary,
  surface = PressureReleaseBoundary,
)
  bathymetry isa Number && (bathymetry = in_units(u"m", bathymetry))
  altimetry isa Number && (altimetry = in_units(u"m", altimetry))
  temperature isa Number && (temperature = in_units(u"°C", temperature))
  salinity isa Number && (salinity = in_units(u"ppt", salinity))
  density isa Number && (density = in_units(u"kg/m^3", density))
  density = something(density, water_density(temperature, salinity))
  soundspeed isa Number && (soundspeed = in_units(u"m/s", soundspeed))
  soundspeed = something(soundspeed, UnderwaterAcoustics.soundspeed(temperature, salinity))
  UnderwaterEnvironment(
    bathymetry, altimetry, temperature, salinity, pH, soundspeed, density, seabed, surface
  )
end

function Base.show(io::IO, env::UnderwaterEnvironment)
  println(io, "UnderwaterEnvironment(")
  println(io, "  bathymetry = $(env.bathymetry), ")
  println(io, "  altimetry = $(env.altimetry), ")
  println(io, "  temperature = $(env.temperature), ")
  println(io, "  salinity = $(env.salinity), ")
  println(io, "  pH = $(env.pH), ")
  println(io, "  soundspeed = $(env.soundspeed), ")
  println(io, "  density = $(env.density), ")
  println(io, "  seabed = $(env.seabed), ")
  println(io, "  surface = $(env.surface), ")
  println(io, ")")
end

"""
    is_range_dependent(env)

Return `true` if any quantity (e.g. sound speed, bathymetry, etc) in the
environment `env` depends on the horizontal location, and `false` otherwise.
"""
function is_range_dependent(env::UnderwaterEnvironment)
  for p ∈ (:bathymetry, :altimetry, :temperature, :salinity, :pH, :soundspeed, :density, :seabed, :surface)
    is_range_dependent(getproperty(env, p)) && return true
  end
  false
end

"""
    is_isovelocity(env)

Return `true` if the sound speed in the environment `env` is a constant.
"""
is_isovelocity(env::UnderwaterEnvironment) = is_constant(env.soundspeed)

"""
    env_type(env)

Return the base number type for the environment. Typically, this is a `Float64`,
but could differ if the environment was constructed using other number types.
"""
function env_type(env::UnderwaterEnvironment)
  a = value(env.altimetry, (0,0,0))
  b = value(env.bathymetry, (0,0,0))
  c = value(env.soundspeed, (0,0,0))
  s = value(env.salinity, (0,0,0))
  pH = value(env.pH, (0,0,0))
  d = value(env.density, (0,0,0))
  r1 = real(reflection_coef(env.surface, 1000.0, 0.0, d, c))
  r2 = real(reflection_coef(env.seabed, 1000.0, 0.0, d, c))
  promote_type(typeof(a), typeof(b), typeof(c), typeof(s), typeof(pH), typeof(d), typeof(r1), typeof(r2))
end

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
  f = frequency(tx)
  s = spl(tx)
  p = Tuple(location(tx))
  print(io, "TX[$f Hz, $s dB]$p")
end

function Base.show(io::IO, tx::AbstractAcousticReceiver)
  p = Tuple(location(tx))
  print(io, "RX$p")
end

"""
    AcousticSource(pos, frequency; spl=0)
    AcousticSource(x, z, frequency; spl=0)
    AcousticSource(x, y, z, frequency; spl=0)

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
  function AcousticSource(pos, frequency; spl=0)
    p = xyz(pos)
    f = in_units(u"Hz", frequency)
    s = in_units(u"dB", spl)
    new{typeof(p),typeof(f),typeof(s)}(p, f, s)
  end
end

AcousticSource(x, z, frequency; spl=0) = AcousticSource(xyz(x, z), frequency; spl)
AcousticSource(x, y, z, frequency; spl=0) = AcousticSource(xyz(x, y, z), frequency; spl)

"""
    AcousticReceiver(pos)
    AcousticReceiver(x, z)
    AcousticReceiver(x, y, z)

Receiver at location `pos`.
"""
struct AcousticReceiver{T1} <: AbstractAcousticReceiver
  pos::T1
  function AcousticReceiver(pos)
    p = xyz(pos)
    new{typeof(p)}(p)
  end
end

AcousticReceiver(x, z) = AcousticReceiver(xyz(x, z))
AcousticReceiver(x, y, z) = AcousticReceiver(xyz(x, y, z))

"""
    AcousticReceiverGrid2D(xrange, zrange)

A 2D grid of receivers with the specified location ranges.
"""
struct AcousticReceiverGrid2D{T1,T2,T3} <: AbstractArray{AcousticReceiver{T3},2}
  xrange::T1
  zrange::T2
  function AcousticReceiverGrid2D(xrange, zrange)
    xrange isa Number && (xrange = StepRangeLen(xrange, 0, 1))
    zrange isa Number && (zrange = StepRangeLen(zrange, 0, 1))
    xrange = StepRangeLen(xrange)
    zrange = StepRangeLen(zrange)
    T3 = typeof(xyz(xrange[1], zrange[1]))
    new{typeof(xrange),typeof(zrange),T3}(xrange, zrange)
  end
end

Base.size(g::AcousticReceiverGrid2D) = (g.xrange.len, g.zrange.len)
Base.getindex(g::AcousticReceiverGrid2D, I::Vararg{Int,2}) = AcousticReceiver(g.xrange[I[1]], g.zrange[I[2]])
Base.setindex!(g::AcousticReceiverGrid2D, v, I::Vararg{Int,2}) = error("AcousticReceiverGrid2D is readonly")

"""
    AcousticReceiverGrid3D(xrange, yrange, zrange)

A 3D grid of receivers with the specified location ranges.
"""
struct AcousticReceiverGrid3D{T1,T2,T3,T4} <: AbstractArray{AcousticReceiver{T4},3}
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
    T4 = typeof(xyz(xrange[1], yrange[1], zrange[1]))
    new{typeof(xrange),typeof(yrange),typeof(zrange),T4}(xrange, yrange, zrange)
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
    rand([rng::AbstractRNG, ] noise::AbstractNoiseModel, nsamples; fs)
    rand([rng::AbstractRNG, ] noise::AbstractNoiseModel, nsamples, nchannels; fs)

Generate random noise samples from the noise model `noise` with the specified
size `nsamples` (can be a 2-tuple for multichannel noise). The noise is returned
as a signal sampled at `fs`. The optional `rng` argument specifies the random
number generator to use.
"""
function Base.rand(noise::AbstractNoiseModel, nsamples::Integer; fs)
  Base.rand(default_rng(), noise, nsamples; fs)
end

function Base.rand(noise::AbstractNoiseModel, nsamples::Integer, nch::Integer; fs)
  Base.rand(default_rng(), noise, nsamples, nch; fs)
end

function Base.rand(noise::AbstractNoiseModel, nsamples::NTuple{2,Int64}; fs)
  rand(default_rng(), noise, nsamples[1], nsamples[2]; fs)
end

function Base.rand(rng::AbstractRNG, noise::AbstractNoiseModel, nsamples::NTuple{2,Int64}; fs)
  rand(rng, noise, nsamples[1], nsamples[2]; fs)
end
