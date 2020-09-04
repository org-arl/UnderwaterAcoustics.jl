using SignalAnalysis: signal
using FFTW: ifft!
using DSP: nextfastfft

export SoundSpeedProfile, soundspeed
export Bathymetry, depth, maxdepth
export Altimetry, altitude
export ReflectionModel, reflectioncoef
export UnderwaterEnvironment, altimetry, bathymetry, ssp, salinity, seasurface, seabed, noise
export AcousticSource, AcousticReceiver, location, nominalfrequency, phasor, record, recorder
export PropagationModel, arrivals, transfercoef, transmissionloss, eigenrays, rays
export NoiseModel
export impulseresponse

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
    depth(bathy:Bathymetry, x, y)

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
    altitude(alt:Altimetry, x, y)

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
    nominalfrequency(src::AcousticSource)

Get the complex phasor representation (amplitude & phase) of a narrowband
acoustic source at the nominal frequency.
"""
function phasor end

"""
    record(src::AcousticSource, duration, fs; start=0.0)
    record(noise::NoiseModel, duration, fs; start=0.0)
    record(model::PropagationModel, tx, rx, duration, fs; start=0.0)

Make a recording of an acoustic source or ambient noise. The `start` time and
`duration` are specified in seconds, and the recording is made at a sampling
rate of `fs` Hz.

For an recording of an acoustic source, free space propagation is assumed,
and the recording is made at a nominal range of 1 meter from the acoustic
center of the source.

For a recording through a propagation model, `tx` and `rx` may be single `AcousticSource`
and `AcousticReceiver`, or an array each.
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
time. It functions in a similar way as the `record()` function.

# Examples:
```julia-repl
julia> rec = recorder(pm, tx, rx);
julia> s = rec(1.0, 44100.0; start=0.0);  # make a recording of 1 second at 44.1 kHz
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
  # threaded version of [transfercoef(model, tx1, rx1; mode=mode) for rx1 ∈ rx], seems to be faster than tmap()
  rx1 = first(rx)
  tc1 = transfercoef(model, tx1, rx1; mode=mode)
  tc = Array{typeof(tc1)}(undef, size(rx))
  Threads.@threads for i ∈ eachindex(rx)
    tc[i] = rx1 === rx[i] ? tc1 : transfercoef(model, tx1, rx[i]; mode=mode)
  end
  tc
end

function rays(model::PropagationModel, tx1::AcousticSource, θ::AbstractArray, rmax)
  tmap(θ1 -> rays(model, tx1, θ1, rmax), θ)
end

transmissionloss(model, tx, rx; mode=:coherent) = -amp2db.(abs.(transfercoef(model, tx, rx; mode=mode)))

function recorder(model::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver)
  f = recorder(model, [tx1], [rx1])
  function rec(duration, fs; start=0.0)
    x = f(duration, fs; start=start)
    signal(dropdims(samples(x); dims=2), framerate(x))
  end
end

function recorder(model::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver})
  recorder(model, [tx1], rx)
end

function recorder(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx1::AcousticReceiver)
  f = recorder(model, tx, [rx1])
  function rec(duration, fs; start=0.0)
    x = f(duration, fs; start=start)
    signal(dropdims(samples(x); dims=2), framerate(x))
  end
end

function recorder(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver})
  arr = [arrivals(model, tx1, rx1) for tx1 ∈ tx, rx1 ∈ rx]
  mindelay, maxdelay = extrema(Iterators.flatten([[a1.time for a1 ∈ a] for a ∈ arr]))
  function rec(duration, fs; start=0.0)
    src = [record(tx1, duration + (maxdelay-mindelay) + 1/fs, fs; start=start-maxdelay) for tx1 ∈ tx]
    nsamples = round(Int, duration * fs)
    x = zeros(Base.promote_eltype(src...), nsamples, length(rx))
    for j = 1:length(tx)
      for k = 1:length(rx)
        for a ∈ arr[j,k]
          t = round(Int, (maxdelay - a.time) * fs)
          x[:,k] .+= a.phasor .* src[j][t+1:t+nsamples]
        end
      end
    end
    noisemodel = noise(environment(model))
    if noisemodel !== missing
      for k = 1:length(rx)
        x[:,k] .+= record(noisemodel, duration, fs; start=start)
      end
    end
    signal(x, fs)
  end
end

function record(model::PropagationModel, tx, rx, duration, fs; start=0.0)
  recorder(model, tx, rx)(duration, fs; start=start)
end

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

function tmap(f, itr)
  refs = [Threads.@spawn(f(i)) for i ∈ itr]
  fetch.(refs)
end

### pretty printing

function Base.show(io::IO, env::UnderwaterEnvironment)
  println(io, split(string(typeof(env).name), ".")[end], ":")
  println(io, "  altimetry = ", altimetry(env))
  println(io, "  bathymetry = ", bathymetry(env))
  println(io, "  ssp = ", ssp(env))
  println(io, "  salinity = ", salinity(env))
  println(io, "  seasurface = ", seasurface(env))
  println(io, "  seabed = ", seabed(env))
  println(io, "  noise = ", noise(env))
end

function Base.show(io::IO, model::PropagationModel)
  print(io, split(string(typeof(model).name), ".")[end], " with ")
  show(io, environment(model))
end

function Base.show(io::IO, a::Arrival)
  if isnan(a.time) || isnan(a.phasor)
    @printf(io, "∠%5.1f° %2d↑ %2d↓%s",
      rad2deg(a.launchangle), a.surface, a.bottom, a.raypath === missing ? "" : " ⤷")
  else
    @printf(io, "∠%5.1f° %2d↑ %2d↓ ∠%5.1f° | %6.2f ms | %5.1f dB ϕ%6.1f°%s",
      rad2deg(a.launchangle), a.surface, a.bottom, rad2deg(a.arrivalangle),
      1000*a.time, amp2db(abs(a.phasor)), rad2deg(angle(a.phasor)),
      a.raypath === missing ? "" : " ⤷")
  end
end
