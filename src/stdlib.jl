import Random: AbstractRNG
import SignalAnalysis: RedGaussian, signal, db2amp
import Interpolations: interpolate, extrapolate, Gridded, Flat, Linear, Throw
import Interpolations: linear_interpolation, cubic_spline_interpolation

export FluidBoundary, RigidBoundary, PressureReleaseBoundary, Rock, Pebbles
export SandyGravel, CoarseSand, MediumSand, FineSand, VeryFineSand, ClayeySand
export CoarseSilt, SandySilt, Silt, FineSilt, SandyClay, SiltyClay, Clay
export SeaState0, SeaState1, SeaState2, SeaState3, SeaState4, SeaState5
export SeaState6, SeaState7, SeaState8, SeaState9, WhiteGaussianNoise
export RedGaussianNoise, PekerisWaveguide, WindySurface, SampledField

################################################################################
# standard environments

"""
    PekerisWaveguide(; kwargs...)

Create a Pekeris waveguide environment with the given parameters. All parameters
are optional (have default values). Default values are specified below:

- `h`  = water depth (m)
- `c` = sound speed in water (m/s, computed from temperature and salinity)
- `ρ` = density of water (kg/m³, computed from temperature and salinity)

Other parameters for `UnderwaterEnvironment` may be specified as well. For
example, `seabed`, `surface`, `temperature`, `salinity` and `density` may be
specified.

Returns an underwater environment with the specified parameters to ensure it
is a Pekeris waveguide.
"""
function PekerisWaveguide(; h=100.0, c=nothing, ρ=nothing, kwargs...)
  UnderwaterEnvironment(;
    bathymetry = h,
    soundspeed = c,
    density = ρ,
    kwargs...
  )
end

################################################################################
# boundary conditions

"""
    FluidBoundary(ρ, c, δ)

Create a fluid half-space boundary with density `ρ`, sound speed `c`, and
dimensionless absorption coefficient `δ`.
"""
struct FluidBoundary{T} <: AbstractAcousticBoundary
  ρ::T
  c::T
  δ::T
  function FluidBoundary(ρ, c, δ)
    ρ = in_units(u"kg/m^3", ρ)
    c = in_units(u"m/s", c)
    isinf(c) && return new{typeof(c)}(zero(c), c, zero(c))
    new{typeof(ρ)}(promote(ρ, c, δ)...)
  end
end

FluidBoundary(ρ, c) = FluidBoundary(ρ, c, 0)

function Base.show(io::IO, b::FluidBoundary)
  if isinf(b.c)
    print(io, "RigidBoundary")
  elseif b.c == 0
    print(io, "PressureReleaseBoundary")
  elseif b.δ == 0
    print(io, "FluidBoundary(ρ=$(b.ρ), c=$(b.c))")
  else
    print(io, "FluidBoundary(ρ=$(b.ρ), c=$(b.c), δ=$(b.δ))")
  end
end

function reflection_coef(bc::FluidBoundary, frequency, θ, ρ, c)
  isinf(bc.c) && return 1.0 + 0im
  bc.c == 0 && return -1.0 + 0im
  θ = in_units(u"rad", θ)
  ρ = in_units(u"kg/m^3", ρ)
  c = in_units(u"m/s", c)
  reflection_coef(θ, bc.ρ / ρ, bc.c / c, bc.δ)
end

"""
Rigid boundary condition.
"""
const RigidBoundary = FluidBoundary(0.0, Inf, 0.0)

"""
Pressure-release boundary condition.
"""
const PressureReleaseBoundary = FluidBoundary(0.0, 0.0, 0.0)

# seabeds from APL-UW TR 9407 (1994), IV-6 Table 2
const Rock = FluidBoundary(2.5 * 1023, 2.5 * 1528, 0.01374)
const Pebbles = FluidBoundary(2.5 * 1023, 1.8 * 1528, 0.01374)
const SandyGravel = FluidBoundary(2.492 * 1023, 1.3370 * 1528, 0.01705)
const CoarseSand = FluidBoundary(2.231 * 1023, 1.2503 * 1528, 0.01638)
const MediumSand = FluidBoundary(1.845 * 1023, 1.1782 * 1528, 0.01624)
const FineSand = FluidBoundary(1.451 * 1023, 1.1073 * 1528, 0.01602)
const VeryFineSand = FluidBoundary(1.268 * 1023, 1.0568 * 1528, 0.01875)
const ClayeySand = FluidBoundary(1.224 * 1023, 1.0364 * 1528, 0.02019)
const CoarseSilt = FluidBoundary(1.195 * 1023, 1.0179 * 1528, 0.02158)
const SandySilt = FluidBoundary(1.169 * 1023, 0.9999 * 1528, 0.01261)
const Silt = FluidBoundary(1.149 * 1023, 0.9873 * 1528, 0.00386)
const FineSilt = FluidBoundary(1.148 * 1023, 0.9861 * 1528, 0.00306)
const SandyClay = FluidBoundary(1.147 * 1023, 0.9849 * 1528, 0.00242)
const SiltyClay = FluidBoundary(1.146 * 1023, 0.9824 * 1528, 0.00163)
const Clay = FluidBoundary(1.145 * 1023, 0.98 * 1528, 0.00148)

"""
    WindySurface(windspeed)

Reflection model for a water surface affected by wind. `windspeed` is given
in m/s.
"""
struct WindySurface{T} <: AbstractAcousticBoundary
  windspeed::T
  function WindySurface(windspeed)
    windspeed = in_units(u"m/s", windspeed)
    new{typeof(windspeed)}(windspeed)
  end
end

function reflection_coef(bc::WindySurface, frequency, θ, ρ, c)
  frequency = in_units(u"Hz", frequency)
  θ = in_units(u"rad", θ)
  surface_reflection_coef(bc.windspeed, frequency, θ)
end

# WMO sea states
# from APL-UW TR 9407 (1994), II-4 Table 2 (median windspeed)
const SeaState0 = WindySurface(0.8)
const SeaState1 = WindySurface(2.6)
const SeaState2 = WindySurface(4.4)
const SeaState3 = WindySurface(6.9)
const SeaState4 = WindySurface(9.8)
const SeaState5 = WindySurface(12.6)
const SeaState6 = WindySurface(19.3)
const SeaState7 = WindySurface(26.5)
const SeaState8 = WindySurface(30.6)
const SeaState9 = WindySurface(32.9)

################################################################################
# noise models

"""
    WhiteGaussianNoise(σ)
    WhiteGaussianNoise(psd, fs)

Create a white Gaussian ambient noise model with variance `σ²` µPa² or with
power spectral density `psd` µPa²/Hz and bandwidth `fs/2` Hz.
"""
struct WhiteGaussianNoise{T<:AbstractFloat} <: AbstractNoiseModel
  σ::T
  function WhiteGaussianNoise(σ)
    σ = in_units(u"µPa", σ)
    new{typeof(σ)}(σ)
  end
end

function WhiteGaussianNoise(psd, fs)
  psd = in_units(u"dB", psd)
  fs = in_units(u"Hz", fs)
  WhiteGaussianNoise(db2amp(psd) * sqrt(fs / 2))
end

function Base.rand(rng::AbstractRNG, noise::WhiteGaussianNoise, nsamples, fs)
  fs = in_units(u"Hz", fs)
  signal(randn(rng, nsamples) .* noise.σ, fs)
end

"""
    RedGaussianNoise(σ)

Create an ambient noise model with variance `σ²` µPa² and `1/f²` variation
in power spectral density.
"""
struct RedGaussianNoise{T<:AbstractFloat} <: AbstractNoiseModel
  σ::T
  function RedGaussianNoise(σ)
    σ = in_units(u"µPa", σ)
    new{typeof(σ)}(σ)
  end
end

function Base.rand(rng::AbstractRNG, noise::RedGaussianNoise, nsamples, fs)
  fs = in_units(u"Hz", fs)
  if length(nsamples) == 1
    signal(rand(rng, RedGaussian(σ=noise.σ, n=nsamples[1])), fs)
  else
    x = mapreduce(hcat, 1:nsamples[2]) do _
      rand(rng, RedGaussian(σ=noise.σ, n=nsamples[1]))
    end
    signal(x, fs)
  end
end

################################################################################
# fields

struct SampledFieldZ{T1,T2} <: DepthDependent
  f::T1
  zrange::T2
  interp::Symbol
end

struct SampledFieldX{T1,T2} <: PositionDependent
  f::T1
  xrange::T2
  interp::Symbol
end

struct SampledFieldXZ{T1,T2} <: PositionDependent
  f::T1
  xrange::T2
  zrange::T2
  interp::Symbol
end

struct SampledFieldXY{T1,T2} <: PositionDependent
  f::T1
  xrange::T2
  yrange::T2
  interp::Symbol
end

struct SampledFieldXYZ{T1,T2} <: PositionDependent
  f::T1
  xrange::T2
  yrange::T2
  zrange::T2
  interp::Symbol
end

(v::SampledFieldZ)(z::Number) = v.f(z)
(v::SampledFieldX)(x::Number) = v.f(x)
(v::SampledFieldXZ)(x, z) = v.f(x, z)
(v::SampledFieldXY)(x, y) = v.f(x, y)
(v::SampledFieldXYZ)(x, y, z) = v.f(x, y, z)

(v::SampledFieldZ)(pos::NamedTuple{(:x,:y,:z)}) = v.f(pos.z)
(v::SampledFieldX)(pos::NamedTuple{(:x,:y,:z)}) = v.f(pos.x)
(v::SampledFieldXZ)(pos::NamedTuple{(:x,:y,:z)}) = v.f(pos.x, pos.z)
(v::SampledFieldXY)(pos::NamedTuple{(:x,:y,:z)}) = v.f(pos.x, pos.y)
(v::SampledFieldXYZ)(pos::NamedTuple{(:x,:y,:z)}) = v.f(pos.x, pos.y, pos.z)

Base.show(io::IO, v::SampledFieldZ) = print(io, "SampledField(z-varying, $(length(v.zrange)) samples)")
Base.show(io::IO, v::SampledFieldX) = print(io, "SampledField(x-varying, $(length(v.xrange)) samples)")
Base.show(io::IO, v::SampledFieldXZ) = print(io, "SampledField(xz-varying, $(length(v.xrange))×$(length(v.zrange)) samples)")
Base.show(io::IO, v::SampledFieldXY) = print(io, "SampledField(xy-varying, $(length(v.xrange))×$(length(v.yrange)) samples)")
Base.show(io::IO, v::SampledFieldXYZ) = print(io, "SampledField(xyz-varying, $(length(v.xrange))×$(length(v.yrange))×$(length(v.zrange)) samples)")

"""
    SampledField(v; x)
    SampledField(v; z)
    SampledField(v; x, y)
    SampledField(v; x, z)
    SampledField(v; x, y, z)

Create a sampled field from a data `v` that may depend on position. For 1D
fields, the `x` or `z` coordinate is required, and `v` is a vector. For 2D
fields, the `x` and `y` coordinates or the `x` and `z` coordinates are
required, and `v` is a matrix. For 3D fields, the `x`, `y`, and `z` coordinates
are required, and `v` is a 3D array.

Keyword argument `interp` is used to specify the interpolation method. `:linear`
interpolation is supported for 1D, 2D and 3D fields. For 2D and 3D fields, the
data must be sampled on a regular grid. For uniformly sampled `1D` fields,
`:cubic` interpolation is also supported.
"""
function SampledField(v; x=nothing, y=nothing, z=nothing, interp=:linear)
  # TODO: support PCHIP interpolation (see HLS-2021-01.pdf in OALIB distribution)
  if x === nothing && y === nothing && z !== nothing
    v = float(v)
    z = float(z)
    ndx = sortperm(z)
    if interp === :cubic
      z isa AbstractRange || error("Cubic interpolation requires `z` to be an `AbstractRange`")
      f = cubic_spline_interpolation(z[ndx], v[ndx]; extrapolation_bc=Flat())
    elseif interp === :linear
      f = linear_interpolation(z[ndx], v[ndx]; extrapolation_bc=Flat())
    else
      error("Unsupported interpolation")
    end
    SampledFieldZ(f, z, interp)
  elseif x !== nothing && y === nothing && z === nothing
    interp ∈ (:linear, :cubic) || error("Unsupported interpolation")
    v = float(v)
    x = float(x)
    ndx = sortperm(x)
    if interp === :cubic
      x isa AbstractRange || error("Cubic interpolation requires `x` to be an `AbstractRange`")
      f = cubic_spline_interpolation(x[ndx], v[ndx]; extrapolation_bc=Flat())
    elseif interp === :linear
      f = linear_interpolation(x[ndx], v[ndx]; extrapolation_bc=Flat())
    else
      error("Unsupported interpolation")
    end
    SampledFieldX(f, x, interp)
  elseif x !== nothing && y === nothing && z !== nothing
    interp === :linear || error("Unsupported interpolation")
    v = float(v)
    x = float(x)
    z = float(z)
    f = extrapolate(interpolate((x, z), v, Gridded(Linear())), Flat())
    SampledFieldXZ(f, x, z, interp)
  elseif x !== nothing && y !== nothing && z === nothing
    interp === :linear || error("Unsupported interpolation")
    v = float(v)
    x = float(x)
    y = float(y)
    f = extrapolate(interpolate((x, y), v, Gridded(Linear())), Flat())
    SampledFieldXY(f, x, y, interp)
  elseif x !== nothing && y !== nothing && z !== nothing
    interp === :linear || error("Unsupported interpolation")
    v = float(v)
    x = float(x)
    y = float(y)
    z = float(z)
    f = extrapolate(interpolate((x, y, z), v, Gridded(Linear())), Flat())
    SampledFieldXYZ(f, x, y, z, interp)
  else
    error("Only x, z, x-z, x-y, or x-y-z interpolation supported")
  end
end

