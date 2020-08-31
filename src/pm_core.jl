using SignalAnalysis: signal

export SoundSpeedProfile, soundspeed
export Bathymetry, depth, maxdepth
export Altimetry, altitude
export ReflectionModel, reflectioncoef
export UnderwaterEnvironment, altimetry, bathymetry, ssp, salinity, seasurface, seabed, noise
export AcousticSource, AcousticReceiver, location, nominalfrequency, phasor, record
export PropagationModel, arrivals, transfercoef, transmissionloss, eigenrays, rays
export NoiseModel
export impulseresponse

### interfaces

abstract type SoundSpeedProfile end
function soundspeed end

abstract type Bathymetry end
function depth end
function maxdepth end

abstract type Altimetry end
function altitude end

abstract type ReflectionModel end
function reflectioncoef end

abstract type UnderwaterEnvironment end
function altimetry end
function bathymetry end
function ssp end
function salinity end
function seasurface end
function seabed end
function noise end

abstract type AcousticSource end
function location end
function nominalfrequency end
function phasor end
function record end

abstract type AcousticReceiver end
function location end

abstract type NoiseModel end
function record end

abstract type PropagationModel{T<:UnderwaterEnvironment} end
function environment end
function check end
function arrivals end
function transfercoef end
function transmissionloss end
function eigenrays end
function rays end
function record end

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

# TODO: define scatterers

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

function record(model::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver, duration, fs; start=0.0)
  record(model, [tx1], [rx1], duration, fs; start=start)
end

function record(model::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver)
  record(model, [tx1], [rx1])
end

function record(model::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver}, duration, fs; start=0.0)
  record(model, [tx1], rx, duration, fs; start=start)
end

function record(model::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver})
  record(model, [tx1], rx)
end

function record(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx1::AcousticReceiver, duration, fs; start=0.0)
  signal(dropdims(record(model, tx, [rx1], duration, fs; start=start), 2), fs)
end

function record(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx1::AcousticReceiver)
  dropdims(record(model, tx, [rx1]), 2)
end

function record(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver})
  arr = [arrivals(model, tx1, rx1) for tx1 ∈ tx, rx1 ∈ rx]
  mindelay, maxdelay = extrema(Iterators.flatten([[a1.time for a1 ∈ a] for a ∈ arr]))
  function rec(duration, fs; start=0.0)
    src = [record(tx1, duration + (maxdelay-mindelay), fs; start=start-maxdelay) for tx1 ∈ tx]
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

function record(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver}, duration, fs; start=0.0)
  record(model, tx, rx)(duration, fs; start=start)
end

function impulseresponse(arrivals::Vector{<:Arrival}, fs; reltime=true)
  length(arrivals) == 0 && throw(ArgumentError("No arrivals"))
  mintime, maxtime = extrema(a.time for a ∈ arrivals)
  reltime || (mintime = zero(typeof(arrivals[1].time)))
  ntaps = ceil(Int, (maxtime-mintime) * fs) + 1
  ir = zeros(typeof(arrivals[1].phasor), ntaps)
  for a ∈ arrivals
    ndx = round(Int, (a.time - mintime) * fs) + 1
    ir[ndx] = a.phasor
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
