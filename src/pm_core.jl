export SoundSpeedProfile, soundspeed
export Bathymetry, depth
export Altimetry, altitude
export ReflectionModel, reflectioncoef
export UnderwaterEnvironment, altimetry, bathymetry, ssp, salinity, seasurface, seabed
export AcousticSource, AcousticReceiver, location, nominalfrequency, phasor, record
export PropagationModel, arrivals, transfercoef, eigenrays

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
function checkenv end
function arrivals end
function transfercoef end
function eigenrays end
function record end

### fallbacks

checkenv(model, env) = env

function transfercoef(model::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
  transfercoef(model, tx1, [rx1]; mode=mode)
end

function record(model::PropagationModel, tx1::AcousticSource, rx1::AcousticReceiver, duration, fs; start=0.0)
  record(model, [tx1], [rx1], duration, fs; start=start)
end

function record(model::PropagationModel, tx1::AcousticSource, rx::AbstractArray{AcousticReceiver}, duration, fs; start=0.0)
  record(model, [tx1], rx, duration, fs; start=start)
end

function record(model::PropagationModel, tx::AbstractArray{AcousticSource}, rx1::AcousticReceiver, duration, fs; start=0.0)
  record(model, tx, [rx1], duration, fs; start=start)
end
