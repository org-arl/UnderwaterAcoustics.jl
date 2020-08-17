export IsoSSP, ConstantDepth, ReflectionCoef, FlatSurface, Rayleigh, SurfaceLoss
export Rock, Pebbles, SandyGravel, CoarseSand, MediumSand, FineSand, VeryFineSand
export ClayeySand, CoarseSilt, SandySilt, Silt, FineSilt, SandyClay, SiltyClay, Clay
export SeaState0, SeaState1, SeaState2, SeaState3, SeaState4
export SeaState5, SeaState6, SeaState7, SeaState8, SeaState9

### sound speed profiles

struct IsoSSP{T} <: SoundSpeedProfile
  c::T
end

soundspeed(ssp::IsoSSP, x, y, z) = ssp.c

### bathymetry models

struct ConstantDepth{T} <: Bathymetry
  depth::T
end

depth(bathy::ConstantDepth, x, y) = bathy.depth

### altimetry models

struct FlatSurface <: Altimetry end

altitude(::FlatSurface, x, y) = 0.0

### reflection models

struct ReflectionCoef{T<:Number} <: ReflectionModel
  coef::T
end

reflectioncoef(rm::ReflectionCoef, f, θ) = rm.coef

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

struct SurfaceLoss{T} <: ReflectionModel
  windspeed::T
end

reflectioncoef(rm::SurfaceLoss, f, θ) = -surfaceloss(rm.windspeed, f, θ)

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

Base.@kwdef struct BasicUnderwaterEnvironment{T1<:Altimetry, T2<:Bathymetry, T3<:SoundSpeedProfile, T4<:Number, T5<:ReflectionModel, T6<:ReflectionModel} <: UnderwaterEnvironment
  altimetry::T1 = FlatSurface()
  bathymetry::T2 = ConstantDepth(20.0)
  ssp::T3 = IsoSSP(soundspeed())
  salinity::T4 = 35.0
  seasurface::T5 = SeaState1
  seabed::T6 = SandySilt
end

UnderwaterEnvironment(; kwargs...) = BasicUnderwaterEnvironment(; kwargs...)

altimetry(env::BasicUnderwaterEnvironment) = env.altimetry
bathymetry(env::BasicUnderwaterEnvironment) = env.bathymetry
ssp(env::BasicUnderwaterEnvironment) = env.ssp
salinity(env::BasicUnderwaterEnvironment) = env.salinity
seasurface(env::BasicUnderwaterEnvironment) = env.seasurface
seabed(env::BasicUnderwaterEnvironment) = env.seabed

### basic source & recevier models

struct NarrowbandAcousticSource{T1,T2,T3,T4} <: AcousticSource
  pos::NTuple{3,T1}
  f::T2
  A::T3
  ϕ::T4
end

NarrowbandAcousticSource(x::T, y::T, z::T, f; sourcelevel=1.0, ϕ=0.0) where T = NarrowbandAcousticSource((x, y, z), f, sourcelevel, ϕ)
NarrowbandAcousticSource(x::T, z::T, f; sourcelevel=1.0, ϕ=0.0) where T = NarrowbandAcousticSource((x, zero(T), z), f, sourcelevel, ϕ)
AcousticSource(x::T, y::T, z::T, f; sourcelevel=1.0, ϕ=0.0) where T = NarrowbandAcousticSource((x, y, z), f, sourcelevel, ϕ)
AcousticSource(x::T, z::T, f; sourcelevel=1.0, ϕ=0.0) where T = NarrowbandAcousticSource((x, zero(T), z), f, sourcelevel, ϕ)

location(tx::NarrowbandAcousticSource) = tx.pos
nominalfrequency(tx::NarrowbandAcousticSource) = tx.f
phasor(tx::NarrowbandAcousticSource) = tx.A * cis(tx.ϕ)
record(tx::NarrowbandAcousticSource, start, duration, fs) = tx.A .* sin(2π .* tx.f .* (start:1/fs:start+duration) .+ tx.ϕ)

struct BasicAcousticReceiver{T} <: AcousticReceiver
  pos::NTuple{3,T}
end

BasicAcousticReceiver(x::T, y::T, z::T) where T = BasicAcousticReceiver((x, y, z))
BasicAcousticReceiver(x::T, z::T) where T = BasicAcousticReceiver((x, zero(T), z))
AcousticReceiver(x::T, y::T, z::T) where T = BasicAcousticReceiver((x, y, z))
AcousticReceiver(x::T, z::T) where T = BasicAcousticReceiver((x, zero(T), z))

location(rx::BasicAcousticReceiver) = rx.pos
