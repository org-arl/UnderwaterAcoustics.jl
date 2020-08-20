export SoundSpeedProfile, soundspeed
export Bathymetry, depth, maxdepth
export Altimetry, altitude
export ReflectionModel, reflectioncoef
export UnderwaterEnvironment, altimetry, bathymetry, ssp, salinity, seasurface, seabed
export AcousticSource, AcousticReceiver, location, nominalfrequency, phasor, record
export PropagationModel, arrivals, transfercoef, transmissionloss, eigenrays, rays
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

abstract type AcousticSource end
function location end
function nominalfrequency end
function phasor end
function record end

abstract type AcousticReceiver end
function location end

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

struct RayArrival{T1,T2,T3,T4} <: Arrival
  time::T1
  phasor::T2
  surface::Int
  bottom::Int
  launchangle::T3
  arrivalangle::T3
  raypath::Union{Vector{NTuple{3,T4}},Missing}
end

function RayArrival(time::T1, phasor::T2, surface::Int, bottom::Int, launchangle::T3, arrivalangle::T3, raypath::Vector{NTuple{3,T4}}) where {T1,T2,T3,T4}
  RayArrival{T1,T2,T3,T4}(time, phasor, surface, bottom, launchangle, arrivalangle, raypath)
end

function RayArrival(time::T1, phasor::T2, surface::Int, bottom::Int, launchangle::T3, arrivalangle::T3, ) where {T1,T2,T3}
  RayArrival{T1,T2,T3,Missing}(time, phasor, surface, bottom, launchangle, arrivalangle, missing)
end

### fallbacks

location(x::NTuple{3,T}) where T = x
location(x::NTuple{2,T}) where T = (x[1], zero(T), x[2])

check(model::PropagationModel, env) = env
environment(model::PropagationModel) = model.env

function transfercoef(model::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver}; mode=:coherent)
  # threaded version of [transfercoef(model, tx1, rx1; mode=mode) for rx1 ∈ rx]
  rx1 = first(rx)
  tc1 = transfercoef(model, tx1, rx1; mode=mode)
  tc = Array{typeof(tc1)}(undef, size(rx))
  Threads.@threads for i ∈ eachindex(rx)
    tc[i] = rx1 === rx[i] ? tc1 : transfercoef(model, tx1, rx[i]; mode=mode)
  end
  tc
end

function rays(model::PropagationModel, tx1::AcousticSource, θ::AbstractArray, rmax)
  θ1 = first(θ)
  r1 = rays(model, tx1, θ1, rmax)
  r = Array{typeof(r1)}(undef, size(θ))
  Threads.@threads for i ∈ eachindex(θ)
    r[i] = θ1 === θ[i] ? r1 : rays(model, tx1, θ[i], rmax)
  end
  r
end

transmissionloss(model, tx, rx; mode=:coherent) = -amp2db.(abs.(transfercoef(model, tx, rx; mode=mode)))

function record(model::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver, duration, fs; start=0.0)
  record(model, [tx1], [rx1], duration, fs; start=start)
end

function record(model::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver}, duration, fs; start=0.0)
  record(model, [tx1], rx, duration, fs; start=start)
end

function record(model::PropagationModel, tx::AbstractArray{<:AcousticSource}, rx1::AcousticReceiver, duration, fs; start=0.0)
  record(model, tx, [rx1], duration, fs; start=start)
end

### core implementation

function record(model::PropagationModel, tx::AbstractArray{AcousticSource}, rx::AbstractArray{AcousticReceiver}, duration, fs; start=0.0)
  # TODO implement record
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

### pretty printing

function Base.show(io::IO, env::UnderwaterEnvironment)
  println(io, split(string(typeof(env).name), ".")[end], ":")
  println(io, "  altimetry = ", altimetry(env))
  println(io, "  bathymetry = ", bathymetry(env))
  println(io, "  ssp = ", ssp(env))
  println(io, "  salinity = ", salinity(env))
  println(io, "  seasurface = ", seasurface(env))
  println(io, "  seabed = ", seabed(env))
end

function Base.show(io::IO, model::PropagationModel)
  print(io, split(string(typeof(model).name), ".")[end], " with ")
  show(io, environment(model))
end

function Base.show(io::IO, a::Arrival)
  if a.time === missing || a.phasor === missing
    @printf(io, "∠%5.1f° %2d↑ %2d↓%s",
      rad2deg(a.launchangle), a.surface, a.bottom, a.raypath === missing ? "" : " ⤷")
  else
    @printf(io, "∠%5.1f° %2d↑ %2d↓ ∠%5.1f° | %6.2f ms | %5.1f dB ϕ%6.1f°%s",
      rad2deg(a.launchangle), a.surface, a.bottom, rad2deg(a.arrivalangle),
      1000*a.time, amp2db(abs(a.phasor)), rad2deg(angle(a.phasor)),
      a.raypath === missing ? "" : " ⤷")
  end
end
