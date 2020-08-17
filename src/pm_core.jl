export SoundSpeedProfile, soundspeed
export Bathymetry, depth
export Altimetry, altitude
export ReflectionModel, reflectioncoef
export UnderwaterEnvironment, altimetry, bathymetry, ssp, salinity, seasurface, seabed
export AcousticSource, AcousticReceiver, location, nominalfrequency, phasor, record
export PropagationModel, arrivals, transfercoef, transmissionloss, eigenrays
export impulseresponse

### interfaces

abstract type SoundSpeedProfile end
function soundspeed end

abstract type Bathymetry end
function depth end

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
function checkenvironment end
function arrivals end
function transfercoef end
function transmissionloss end
function eigenrays end
function record end

struct Arrival{T1,T2,T3}
  time::T1
  phasor::T2
  surface::Int
  bottom::Int
  raypath::Union{Vector{NTuple{3,T3}},Missing}
end

function Arrival(time::T1, phasor::T2, surface::Int, bottom::Int, raypath::Vector{NTuple{3,T3}}) where {T1, T2, T3}
  Arrival{T1,T2,T3}(time, phasor, surface, bottom, raypath)
end

function Arrival(time::T1, phasor::T2, surface::Int, bottom::Int) where {T1, T2}
  Arrival{T1,T2,Missing}(time, phasor, surface, bottom, missing)
end

### fallbacks

location(x::NTuple{3,T}) where T = x
location(x::NTuple{2,T}) where T = (x[1], zero(T), x[2])

checkenv(model, env) = env

function transfercoef(model::PropagationModel, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver}; mode=:coherent)
  [transfercoef(model, tx1, rx1; mode=mode) for rx1 ∈ rx]
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
  # TODO
end

function impulseresponse(arrivals::Vector{<:Arrival}, fs; reltime=true)
  length(arrivals) == 0 && throw(ArgumentError("No arrivals"))
  mintime, maxtime = extrema(a.time for a ∈ arrivals)
  reltime || (mintime = zero(typeof(arrivals[1].time)))
  ntaps = ceil(Int, (maxtime-mintime) * fs)
  ir = zeros(typeof(arrivals[1].phasor), ntaps)
  for a ∈ arrivals
    # TODO: think about whether nearest point is the best approach
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
  print(io, "(", a.surface, "/", a.bottom, ") ", round(amp2db(abs(a.phasor)); digits=1),
    " dB ∠", round(rad2deg(angle(a.phasor))), "° @ ", round(1000*a.time; digits=2), " ms")
end
