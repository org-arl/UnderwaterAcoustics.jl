using SignalAnalysis: signal
using FFTW: ifft!
using DSP: nextfastfft
using ToeplitzMatrices: TriangularToeplitz

export SoundSpeedProfile, soundspeed
export Bathymetry, depth, maxdepth
export Altimetry, altitude
export ReflectionModel, reflectioncoef
export UnderwaterEnvironment, altimetry, bathymetry, ssp, salinity, seasurface, seabed, noise
export AcousticSource, AcousticReceiver, location, nominalfrequency, phasor, record, recorder
export PropagationModel, arrivals, transfercoef, transmissionloss, eigenrays, rays
export NoiseModel
export impulseresponse, channelmatrix

### interface: SoundSpeedProfile

abstract type SoundSpeedProfile end

"""
    soundspeed(ssp::SoundSpeedProfile, x, y, z)

Get sound speed at location (`x`, `y`, `z`). If a sound speed profile is range
independent, `x` and `y` may be ignored. `z` is generally negative, since the
sea surface is the datum and z-axis points upwards.
"""
function soundspeed end

### interface: Bathymetry

abstract type Bathymetry end

"""
    depth(bathy::Bathymetry, x, y)

Get water depth at location (`x`, `y`).
"""
function depth end

"""
    maxdepth(bathy::Bathymetry)

Get the maximum water depth.
"""
function maxdepth end

### interface: Altimetry

abstract type Altimetry end

"""
    altitude(alt::Altimetry, x, y)

Get water surface altitude at location (`x`, `y`). The nominal water surface
is considered to have an altitude of zero. However, the water surface may
not be flat, and the Altimetry provisions for variations of altitude around the
nominal altitutde of zero.
"""
function altitude end

### interface: ReflectionModel

abstract type ReflectionModel end

"""
    reflectioncoef(rm::ReflectionModel, f, θ)

Get complex reflection coefficient at frequency `f` Hz and incidence angle `θ`
(w.r.t. the surface normal).
"""
function reflectioncoef end

### interface: UnderwaterEnvironment

abstract type UnderwaterEnvironment end

"""
    altimetry(env::UnderwaterEnvironment)::Altimetry

Get the altimetry for the underwater environment.
"""
function altimetry end

"""
    bathymetry(env::UnderwaterEnvironment)::Bathymetry

Get the bathymetry for the underwater environment.
"""
function bathymetry end

"""
    ssp(env::UnderwaterEnvironment)::SoundSpeedProfile

Get the sound speed profile for the underwater environment.
"""
function ssp end

"""
    salinity(env::UnderwaterEnvironment)

Get the salinity of the underwater environment.
"""
function salinity end

"""
    waterdensity(env::UnderwaterEnvironment)

Get the nominal water density for the underwater environment.
"""
function waterdensity end

"""
    seasurface(env::UnderwaterEnvironment)::ReflectionModel

Get the sea surface reflection model for the underwater environment.
"""
function seasurface end

"""
    seabed(env::UnderwaterEnvironment)::ReflectionModel

Get the seabed reflection model for the underwater environment.
"""
function seabed end

"""
    noise(env::UnderwaterEnvironment)::NoiseModel

Get the noise model for the underwater environment.
"""
function noise end

### interface: AcousticSource

abstract type AcousticSource end

"""
    location(src::AcousticSource)
    location(src::AcousticReceiver)

Get the location of an acoustic source or receiver as a 3-tuple (`x`, `y`, `z`).
"""
function location end

"""
    nominalfrequency(src::AcousticSource)

Get the nominal frequency of an acoustic source in Hz.
"""
function nominalfrequency end

"""
    phasor(src::AcousticSource)

Get the complex phasor representation (amplitude & phase) of a narrowband
acoustic source at the nominal frequency.
"""
function phasor end

"""
    record(src::AcousticSource, duration, fs; start=0.0)
    record(noise::NoiseModel, duration, fs; start=0.0)
    record(model::PropagationModel, tx, rx, duration, fs; start=0.0)
    record(model::PropagationModel, tx, rx, sig; reltime=true)

Make a recording of an acoustic source or ambient noise. The `start` time and
`duration` are specified in seconds, and the recording is made at a sampling
rate of `fs` Hz.

For an recording of an acoustic source, free space propagation is assumed,
and the recording is made at a nominal range of 1 meter from the acoustic
center of the source.

For a recording through a propagation model, `tx` and `rx` may be single `AcousticSource`
and `AcousticReceiver`, or an array each. The returned signal is always complex,
irrespective of whether the source is real or complex.

When a `sig` is specified, the sources are assumed to transmit the sampled signal in `sig`.
The number of channels in `sig` must match the number of sources. The returned signal
is the same type as the input signal (real or complex). If `reltime` is `true`, the
recorded signal starts at the first arrival, otherwise it starts at the beginning of the
transmission.
"""
function record end

### interface: AcousticReceiver

abstract type AcousticReceiver end
function location end

### interface: NoiseModel

abstract type NoiseModel end
function record end

### interface: PropagationModel

abstract type PropagationModel{T<:UnderwaterEnvironment} end

"""
    environment(pm::PropagationModel)

Get the environment associated with the propagation model.
"""
function environment end

"""
    check(pm::Type{<:PropagationModel}, env::UnderwaterEnvironment)
    check(pm::Type{<:PropagationModel}, env=missing)

Check if an propagation model is available, and can simulate the specified
environment. Returns the environment if it can be simulated, or throws an
error with a descriptive error message if it cannot be simulated.

This function is internally used by the propagation modeling toolbox to choose
a model or offer a selection of models to the user.
"""
function check end

"""
    arrivals(pm::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver)

Compute the arrivals from `tx1` to `rx1`. Returns an array of `Arrival` structs.
"""
function arrivals end

"""
    transfercoef(pm::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
    transfercoef(pm::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver}; mode=:coherent)

Compute the complex transfer coefficients from `tx1` to `rx1` or all receivers in `rx`.
The mode may be `:coherent` or `:incoherent`.
"""
function transfercoef end

"""
    transmissionloss(pm::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
    transmissionloss(pm::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver}; mode=:coherent)

Compute the transmission loss in dB from `tx1` to `rx1` or all receivers in `rx`.
The mode may be `:coherent` or `:incoherent`.
"""
function transmissionloss end

"""
    eigenrays(pm::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver)

Compute the eigenrays from `tx1` to `rx1`. Returns an array of `RayArrival` structs.
"""
function eigenrays end

"""
    rays(pm::PropagationModel, tx1::AcousticSource, θ::Real, rmax)

Compute the rays from `tx1` launched at angle `θ` (or all angles in `θ`, if it is a vector).
Returns an array of `RayArrival` datatypes. `rmax` is the maximum horizontal range in meters
to track the rays over.
"""
function rays end

function record end

"""
    recorder(model::PropagationModel, tx, rx)

Create a recorder function that may be called later to make an acoustic recording
of sources in `tx` at receviers `rx`. `tx` and `rx` may be single `AcousticSource` and
`AcousticReceiver`, or an array each.

The recorder function may be called later with `duration`, `fs`, and optionally a `start`
time. Alternatively, the recorder function may also be called with a sampled signal.
It functions in a similar way as the `record()` function.

# Examples:
```julia-repl
julia> rec = recorder(pm, tx, rx);
julia> s = rec(1.0, 44100.0; start=0.0);  # make a recording of 1 second at 44.1 kHz
julia> s = rec(signal(randn(44100), 44100));  # transmit a random 1 second signal
```
"""
function recorder end

### arrival types

abstract type Arrival end

struct RayArrival{T1,T2} <: Arrival
  time::T1
  phasor::Complex{T1}
  surface::Int
  bottom::Int
  launchangle::T1
  arrivalangle::T1
  raypath::Union{Vector{NTuple{3,T2}},Missing}
end

function RayArrival(time::T1, phasor::Complex{T1}, surface::Int, bottom::Int, launchangle::T1, arrivalangle::T1, raypath::Vector{NTuple{3,T2}}) where {T1,T2}
  RayArrival{T1,T2}(time, phasor, surface, bottom, launchangle, arrivalangle, raypath)
end

function RayArrival(time, phasor, surface::Int, bottom::Int, launchangle, arrivalangle, raypath::Vector{NTuple{3,T}}) where T
  t, r, i, θ1, θ2 = promote(time, real(phasor), imag(phasor), launchangle, arrivalangle)
  RayArrival{typeof(t),T}(t, Complex(r, i), surface, bottom, θ1, θ2, raypath)
end

function RayArrival(time, phasor, surface::Int, bottom::Int, launchangle, arrivalangle)
  t, r, i, θ1, θ2 = promote(time, real(phasor), imag(phasor), launchangle, arrivalangle)
  RayArrival{typeof(t),Missing}(t, Complex(r, i), surface, bottom, θ1, θ2, missing)
end

phasortype(::Type{RayArrival{T1,T2}}) where {T1,T2} = T1

### fallbacks & helpers

location(x::NTuple{3,T}) where T = x
location(x::NTuple{2,T}) where T = (x[1], zero(T), x[2])

check(model::PropagationModel, env) = env
environment(model::PropagationModel) = model.env

function transfercoef(model::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
  arr = arrivals(model, tx1, rx1)
  length(arr) == 0 && return zero(phasortype(eltype(arr)))
  if mode === :coherent
    f = nominalfrequency(tx1)
    tc = sum(a.phasor * cis(2π * a.time * f) for a ∈ arr)
  elseif mode === :incoherent
    tc = Complex(√sum(abs2(a.phasor) for a ∈ arr), 0)
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  tc
end

function transfercoef(model::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver}; mode=:coherent)
  tmap(rx1 -> transfercoef(model, tx1, rx1; mode), rx)
end

function rays(model::PropagationModel, tx1::AcousticSource, θ::AbstractArray, rmax)
  tmap(θ1 -> rays(model, tx1, θ1, rmax), θ)
end

transmissionloss(model, tx, rx; mode=:coherent) = -amp2db.(abs.(transfercoef(model, tx, rx; mode)))

struct Recorder{T1,T2,T3,T4}
  noisemodel::T1
  tx::T2  # always an array of sources
  rx::T3  # could be an array or a single recevier
  arr::Matrix{T4}
end

function (rec::Recorder)(duration, fs; start=0.0)
  mindelay, maxdelay = extrema(Iterators.flatten([[a1.time for a1 ∈ a] for a ∈ rec.arr]))
  src = [analytic(record(tx1, duration + (maxdelay-mindelay) + 1/fs, fs; start=start-maxdelay)) for tx1 ∈ rec.tx]
  nsamples = round(Int, duration * fs)
  x = zeros(Base.promote_eltype(src...), nsamples, size(rec.arr, 2))
  for j = 1:size(rec.arr, 1)
    for k = 1:size(rec.arr, 2)
      for a ∈ rec.arr[j,k]
        t = round(Int, (maxdelay - a.time) * fs)
        x[:,k] .+= a.phasor .* src[j][t+1:t+nsamples]
      end
    end
  end
  if rec.noisemodel !== missing && rec.noisemodel !== nothing
    for k = 1:size(rec.arr, 2)
      x[:,k] .+= record(rec.noisemodel, duration, fs; start)
    end
  end
  x̄ = rec.rx isa AbstractArray ? x : dropdims(x; dims=2)
  signal(x̄, fs)
end

function (rec::Recorder)(sig; fs=framerate(sig), reltime=true)
  nchannels(sig) == length(rec.tx) || throw(ArgumentError("Input signal must have $(length(rec.tx)) channel(s)"))
  mindelay, maxdelay = extrema(Iterators.flatten([[a1.time for a1 ∈ a] for a ∈ rec.arr]))
  n1 = round(Int, maxdelay * fs)
  n2 = round(Int, (maxdelay - mindelay) * fs) + 1
  src = [analytic(vcat(zeros(eltype(sig), n1), sig[:,i], zeros(eltype(sig), n2))) for i ∈ eachindex(rec.tx)]
  nsamples = nframes(sig) + n1
  x = zeros(Base.promote_eltype(src...), nsamples, size(rec.arr, 2))
  for j = 1:size(rec.arr, 1)
    for k = 1:size(rec.arr, 2)
      for a ∈ rec.arr[j,k]
        t = round(Int, (maxdelay - a.time) * fs)
        x[:,k] .+= a.phasor .* src[j][t+1:t+nsamples]
      end
    end
  end
  if rec.noisemodel !== missing && rec.noisemodel !== nothing
    for k = 1:size(rec.arr, 2)
      x[:,k] .+= record(rec.noisemodel, nsamples/fs, fs)
    end
  end
  n3 = reltime ? round(Int, mindelay * fs) : 1
  x̄ = rec.rx isa AbstractArray ? @view(x[n3:end,:]) : dropdims(@view x[n3:end,:]; dims=2)
  isanalytic(sig) ? signal(x̄, fs) : signal(convert.(eltype(sig), real.(x̄) .* √2), fs)
end

"""
    channelmatrix(rec::Recorder, fs, ntaps=0; tx=1, rx=1, approx=false)
    channelmatrix(rec::Vector{<:Arrival}, fs, ntaps=0; approx=false)

Generate a sampled channel matrix at a sampling rate of `fs` Hz. If `ntaps` is zero,
the number of taps of the channel matrix are chosen automatically.

If `approx` is `true`, a fast algorithm is used to generate a sparse channel matrix
that assigns an arrival to the nearest sampling time.
"""
function channelmatrix(rec::Recorder, fs, ntaps=0; tx=1, rx=1, approx=false)
  ir = impulseresponse(rec.arr[tx,rx], fs, ntaps; reltime=true, approx)
  TriangularToeplitz(ir, :L)
end

function channelmatrix(arrivals::Vector{<:Arrival}, fs, ntaps=0; approx=false)
  ir = impulseresponse(arrivals, fs, ntaps; reltime=true, approx)
  TriangularToeplitz(ir, :L)
end

function recorder(model::PropagationModel, tx::AcousticSource, rx::AcousticReceiver)
  arr = [arrivals(model, tx1, rx1) for tx1 ∈ [tx], rx1 ∈ [rx]]
  Recorder(noise(environment(model)), [tx], rx, arr)
end

function recorder(model::PropagationModel, tx::AcousticSource, rx::AbstractArray{<:AcousticReceiver})
  arr = [arrivals(model, tx1, rx1) for tx1 ∈ [tx], rx1 ∈ rx]
  Recorder(noise(environment(model)), [tx], rx, arr)
end

function recorder(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx::AcousticReceiver)
  arr = [arrivals(model, tx1, rx1) for tx1 ∈ tx, rx1 ∈ [rx]]
  Recorder(noise(environment(model)), tx, rx, arr)
end

function recorder(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver})
  arr = [arrivals(model, tx1, rx1) for tx1 ∈ tx, rx1 ∈ rx]
  Recorder(noise(environment(model)), tx, rx, arr)
end

# delegate recording task to an ephemeral recorder
record(model::PropagationModel, tx, rx, duration, fs; start=0.0) = recorder(model, tx, rx)(duration, fs; start)
record(model::PropagationModel, tx, rx, sig; reltime=true) = recorder(model, tx, rx)(sig; reltime)

"""
$(SIGNATURES)
Convert a vector of arrivals to a sampled impulse response time series at a
sampling rate of `fs` Hz. If `ntaps` is zero, the number of taps of the impulse
response are chosen automatically.

If `reltime` is `true`, the impulse response start time is relative to the
first arrival, otherwise it is relative to the absolute time. If `approx`
is `true`, a fast algorithm is used to generate a sparse impulse response
that assigns an arrival to the nearest sampling time.
"""
function impulseresponse(arrivals::Vector{<:Arrival}, fs, ntaps=0; reltime=false, approx=false)
  length(arrivals) == 0 && throw(ArgumentError("No arrivals"))
  mintime, maxtime = extrema(a.time for a ∈ arrivals)
  reltime || (mintime = zero(typeof(arrivals[1].time)))
  mintaps = ceil(Int, (maxtime-mintime) * fs) + 1
  if approx
    ntaps == 0 && (ntaps = mintaps)
    ir = zeros(typeof(arrivals[1].phasor), ntaps)
    for a ∈ arrivals
      ndx = round(Int, (a.time - mintime) * fs) + 1
      ndx ≤ ntaps && (ir[ndx] = a.phasor)
    end
  else
    N = nextfastfft(4*max(256, mintaps))
    x = zeros(typeof(arrivals[1].phasor), N)
    for a ∈ arrivals
      δ = (a.time - mintime) * fs
      for i ∈ 1:N
        x[i] += a.phasor * cis(-2π * (i-1)/N * δ)
      end
    end
    ifft!(x)
    if ntaps == 0
      ndx = findfirst(abs.(x[mintaps+1:end]) .≤ abs(arrivals[1].phasor)/100)
      ndx === nothing && (ndx = argmin(abs.(x[mintaps+1:end])))
      ir = x[1:mintaps+ndx]
    elseif ntaps ≤ N/2
      ir = x[1:ntaps]
    else
      Nby2 = floor(Int, N/2)
      ir = vcat(x[1:Nby2], zeros(typeof(arrivals[1].phasor), ntaps-Nby2))
    end
  end
  ir
end

# fast threaded map, assuming all entries have the same result type
function tmap(f, x)
  x1 = first(x)
  y1 = f(x1)
  y = Array{typeof(y1)}(undef, size(x))
  Threads.@threads for i ∈ eachindex(x)
    y[i] = x1 === x[i] ? y1 : f(x[i])
  end
  y
end

function envrealtype(env::UnderwaterEnvironment)
  a = altitude(altimetry(env), 0.0, 0.0)
  b = depth(bathymetry(env), 0.0, 0.0)
  c = soundspeed(ssp(env), 0.0, 0.0, 0.0)
  s = salinity(env)
  d = waterdensity(env)
  r1 = real(reflectioncoef(seasurface(env), 1000.0, 0.0))
  r2 = real(reflectioncoef(seabed(env), 1000.0, 0.0))
  promote_type(typeof(a), typeof(b), typeof(c), typeof(s), typeof(d), typeof(r1), typeof(r2))
end

### pretty printing

function Base.show(io::IO, env::UnderwaterEnvironment)
  s = replace(replace(string(typeof(env)), r"\{.*$" => ""), r"^[^\.]*\." => "")
  println(io, s, ":")
  println(io, "  altimetry = ", altimetry(env))
  println(io, "  bathymetry = ", bathymetry(env))
  println(io, "  ssp = ", ssp(env))
  println(io, "  salinity = ", salinity(env))
  println(io, "  waterdensity = ", waterdensity(env))
  println(io, "  seasurface = ", seasurface(env))
  println(io, "  seabed = ", seabed(env))
  println(io, "  noise = ", noise(env))
end

function Base.show(io::IO, model::PropagationModel)
  s = replace(string(typeof(model)), r"\{.*$" => "")
  print(io, s, " with ")
  show(io, environment(model))
end

function Base.show(io::IO, a::Arrival)
  if isnan(a.time) || isnan(a.phasor)
    @printf(io, "∠%5.1f° %2d↑ %2d↓%s",
      rad2deg(a.launchangle), a.surface, a.bottom, a.raypath === missing || length(a.raypath) == 0 ? "" : " ⤷")
  else
    @printf(io, "∠%5.1f° %2d↑ %2d↓ ∠%5.1f° | %6.2f ms | %5.1f dB ϕ%6.1f°%s",
      rad2deg(a.launchangle), a.surface, a.bottom, rad2deg(a.arrivalangle),
      1000*a.time, amp2db(abs(a.phasor)), rad2deg(angle(a.phasor)),
      a.raypath === missing || length(a.raypath) == 0 ? "" : " ⤷")
  end
end

function Base.show(io::IO, rec::Recorder)
  s = replace(string(typeof(rec)), r"\{.*$" => "")
  print(io, "$s($(size(rec.arr, 1)) => $(size(rec.arr, 2)))")
end

