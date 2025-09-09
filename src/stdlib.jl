import Random: AbstractRNG
import SignalAnalysis: RedGaussian, signal, db2amp
import Interpolations: interpolate, extrapolate, Gridded, Flat, Throw, Interpolations
import Interpolations: linear_interpolation, cubic_spline_interpolation

export FluidBoundary, RigidBoundary, PressureReleaseBoundary
export Rock, Pebbles, SandyGravel, VeryCoarseSand, MuddySandyGravel, CoarseSand
export GravellyMuddySand, MediumSand, MuddyGravel, FineSand, MuddySand, VeryFineSand
export ClayeySand, CoarseSilt, SandySilt, MediumSilt, SandyMud, FineSilt, SandyClay
export VeryFineSilt, SiltyClay, Clay, MultilayerElasticBoundary, Linear, CubicSpline
export SeaState0, SeaState1, SeaState2, SeaState3, SeaState4, SeaState5
export SeaState6, SeaState7, SeaState8, SeaState9, WhiteGaussianNoise
export RedGaussianNoise, WindySurface, SampledField, ElasticBoundary

################################################################################
# boundary conditions

"""
    FluidBoundary(ρ, c, δ=0, σ=0)

Create a fluid half-space boundary with density `ρ`, sound speed `c`,
dimensionless absorption coefficient `δ`, and interfacial roughness `σ`.

# Examples
```julia-repl
julia> FluidBoundary(1200, 1500)
FluidBoundary(ρ=1200.0, c=1500.0)

julia> FluidBoundary(1200, 1500, 0.1)
FluidBoundary(ρ=1200.0, c=1500.0, δ=0.1)

julia> FluidBoundary(1200, 1500, 0.1, 0.1)
FluidBoundary(ρ=1200.0, c=1500.0, δ=0.1, σ=0.1)

julia> FluidBoundary(ρ=1200, c=1500)
FluidBoundary(ρ=1200.0, c=1500.0)

julia> FluidBoundary(ρ=1200, c=1500, δ=0.1, σ=0.1)
FluidBoundary(ρ=1200.0, c=1500.0, δ=0.1, σ=0.1)

julia> FluidBoundary(ρ=1200.0, c=0.0)
PressureReleaseBoundary

julia> FluidBoundary(ρ=1200.0, c=Inf)
RigidBoundary

julia> FluidBoundary(1.2u"g/cm^3", 1500u"m/s")
FluidBoundary(ρ=1200.0, c=1500.0)
```
"""
struct FluidBoundary{T} <: AbstractAcousticBoundary
  ρ::T
  c::T
  δ::T
  σ::T
  function FluidBoundary(ρ, c, δ, σ)
    ρ = in_units(u"kg/m^3", ρ)
    c = in_units(u"m/s", c)
    σ = in_units(u"m", σ)
    isinf(c) && return new{typeof(c)}(zero(c), c, zero(c), zero(c))
    ρ, c, δ, σ = float.(promote(ρ, c, δ, σ))
    new{typeof(ρ)}(ρ, c, δ, σ)
  end
end

FluidBoundary(ρ, c) = FluidBoundary(ρ, c, 0, 0)
FluidBoundary(ρ, c, δ) = FluidBoundary(ρ, c, δ, 0)
FluidBoundary(; ρ, c, δ=0, σ=0) = FluidBoundary(ρ, c, δ, σ)

function Base.show(io::IO, b::FluidBoundary)
  if isinf(b.c)
    print(io, "RigidBoundary")
  elseif b.c == 0
    print(io, "PressureReleaseBoundary")
  else
    print(io, "FluidBoundary(ρ=$(b.ρ), c=$(b.c)")
    b.δ == 0 || print(io, ", δ=$(b.δ)")
    b.σ == 0 || print(io, ", σ=$(b.σ)")
    print(io, ")")
  end
end

function reflection_coef(bc::FluidBoundary, frequency, θ, ρ, c)
  isinf(bc.c) && return 1.0 + 0im
  bc.c == 0 && return -1.0 + 0im
  θ = in_units(u"rad", θ)
  ρ = in_units(u"kg/m^3", ρ)
  c = in_units(u"m/s", c)
  bc.σ == 0 || @warn "Roughness is currently ignored in computing reflection coefficient" maxlog=1
  reflection_coef(θ, bc.ρ / ρ, bc.c / c, bc.δ)
end

"""
Rigid boundary condition.
"""
const RigidBoundary = FluidBoundary(Inf, Inf, 0.0)

"""
Pressure-release boundary condition.
"""
const PressureReleaseBoundary = FluidBoundary(0.0, 0.0, 0.0)

# seabeds from APL-UW TR 9407 (1994), IV-6 Table 2
# ratios evaluated wrt reference ρ=1023, c=1528
const Rock              = FluidBoundary(2557.50, 3820.00, 0.01374)
const Pebbles           = FluidBoundary(2557.50, 2750.40, 0.01374)
const SandyGravel       = FluidBoundary(2549.32, 2073.50, 0.01705)
const VeryCoarseSand    = FluidBoundary(2455.20, 1996.64, 0.01667)
const MuddySandyGravel  = FluidBoundary(2367.22, 1952.48, 0.0163)
const CoarseSand        = FluidBoundary(2282.31, 1910.46, 0.01638)
const GravellyMuddySand = FluidBoundary(2200.47, 1870.42, 0.01645)
const MediumSand        = FluidBoundary(1887.44, 1800.29, 0.01624)
const MuddyGravel       = FluidBoundary(1652.15, 1741.31, 0.0161)
const FineSand          = FluidBoundary(1484.37, 1691.95, 0.01602)
const MuddySand         = FluidBoundary(1369.80, 1650.24, 0.01728)
const VeryFineSand      = FluidBoundary(1297.16, 1614.79, 0.01875)
const ClayeySand        = FluidBoundary(1252.15, 1583.62, 0.02019)
const CoarseSilt        = FluidBoundary(1222.49, 1555.35, 0.02158)
const SandySilt         = FluidBoundary(1195.89, 1527.85, 0.01261)
const MediumSilt        = FluidBoundary(1175.43, 1510.43, 0.00676)
const SandyMud          = FluidBoundary(1175.43, 1508.59, 0.00386)
const FineSilt          = FluidBoundary(1174.40, 1506.76, 0.00306)
const SandyClay         = FluidBoundary(1173.38, 1504.93, 0.00242)
const VeryFineSilt      = FluidBoundary(1173.38, 1503.09, 0.00194)
const SiltyClay         = FluidBoundary(1172.36, 1501.11, 0.00163)
const Clay              = FluidBoundary(1171.34, 1497.44, 0.00148)

"""
    ElasticBoundary(ρ, cₚ, cₛ)
    ElasticBoundary(ρ, cₚ, cₛ, δₚ, δₛ)
    ElasticBoundary(ρ, cₚ, cₛ, δₚ, δₛ, σ)
    ElasticBoundary(b::FluidBoundary)
    ElasticBoundary(b::FluidBoundary, δₛ)
    ElasticBoundary(b::FluidBoundary, cₛ, δₛ)

Create a solid half-space boundary with density `ρ`, compressional sound
speed `cₚ`, shear wave speed `cₛ`, dimensionless compressional absorption
coefficient `δₚ`, dimensionless shear absorption coefficient `δₛ`, and
interfacial roughness `σ`. If the absorption coefficients or interfacial
roughness are unspecified, they are assumed to be 0.

An `ElasticBoundary` may also be constructed from a `b::FluidBoundary` by adding
shear wave speed `cₛ` and shear absorption coefficient `δₛ`. If `cₛ` is not
specified, it is computed using `shearspeed(b.cₚ)`.

# Examples
```julia-repl
julia> ElasticBoundary(1200, 1500, 500)
ElasticBoundary(ρ=1200.0, cₚ=1500.0, cₛ=500.0)

julia> ElasticBoundary(1200, 1500, 500, 0.1, 0.2)
ElasticBoundary(ρ=1200.0, cₚ=1500.0, cₛ=500.0, δₚ=0.1, δₛ=0.2)

julia> ElasticBoundary(1200, 1500, 500, 0.1, 0.2, 0.1)
ElasticBoundary(ρ=1200.0, cₚ=1500.0, cₛ=500.0, δₚ=0.1, δₛ=0.2, σ=0.1)

julia> ElasticBoundary(ρ=1200, cₚ=1500, cₛ=500)
ElasticBoundary(ρ=1200.0, cₚ=1500.0, cₛ=500.0)

julia> ElasticBoundary(ρ=1200, cₚ=1500, cₛ=500, δₚ=0.1, δₛ=0.2)
ElasticBoundary(ρ=1200.0, cₚ=1500.0, cₛ=500.0, δₚ=0.1, δₛ=0.2)

julia> ElasticBoundary(ρ=1200, cₚ=1500, cₛ=500, δₚ=0.1, δₛ=0.2, σ=0.1)
ElasticBoundary(ρ=1200.0, cₚ=1500.0, cₛ=500.0, δₚ=0.1, δₛ=0.2, σ=0.1)

julia> ElasticBoundary(1.2u"g/cm^3", 1500u"m/s", 500u"m/s")
ElasticBoundary(ρ=1200.0, cₚ=1500.0, cₛ=500.0)

julia> ElasticBoundary(FineSand)
ElasticBoundary(ρ=1484.373, cₚ=1691.954, cₛ=414.413, δₚ=0.01602)

julia> ElasticBoundary(FineSand, 0.1)
ElasticBoundary(ρ=1484.373, cₚ=1691.954, cₛ=414.413, δₚ=0.01602, δₛ=0.1)

julia> ElasticBoundary(FineSand, 500, 0.1)
ElasticBoundary(ρ=1484.373, cₚ=1691.954, cₛ=500.0, δₚ=0.01602, δₛ=0.1)
```
"""
struct ElasticBoundary{T} <: AbstractAcousticBoundary
  ρ::T
  cₚ::T
  cₛ::T
  δₚ::T
  δₛ::T
  σ::T
  function ElasticBoundary(ρ, cₚ, cₛ, δₚ, δₛ, σ)
    ρ = in_units(u"kg/m^3", ρ)
    cₚ = in_units(u"m/s", cₚ)
    cₛ = in_units(u"m/s", cₛ)
    σ = in_units(u"m", σ)
    ρ, cₚ, cₛ, δₚ, δₛ, σ = float.(promote(ρ, cₚ, cₛ, δₚ, δₛ, σ))
    new{typeof(ρ)}(ρ, cₚ, cₛ, δₚ, δₛ, σ)
  end
end

ElasticBoundary(ρ, cₚ, cₛ) = ElasticBoundary(ρ, cₚ, cₛ, 0, 0, 0)
ElasticBoundary(ρ, cₚ, cₛ, δₚ, δₛ) = ElasticBoundary(ρ, cₚ, cₛ, δₚ, δₛ, 0)
ElasticBoundary(; ρ, cₚ, cₛ, δₚ=0, δₛ=0, σ=0) = ElasticBoundary(ρ, cₚ, cₛ, δₚ, δₛ, σ)
ElasticBoundary(b::FluidBoundary) = ElasticBoundary(b.ρ, b.c, shearspeed(b.c), b.δ, 0, b.σ)
ElasticBoundary(b::FluidBoundary, δₛ) = ElasticBoundary(b.ρ, b.c, shearspeed(b.c), b.δ, δₛ, b.σ)
ElasticBoundary(b::FluidBoundary, cₛ, δₛ) = ElasticBoundary(b.ρ, b.c, cₛ, b.δ, δₛ, b.σ)

function Base.show(io::IO, b::ElasticBoundary)
  print(io, "ElasticBoundary(ρ=$(b.ρ), cₚ=$(b.cₚ), cₛ=$(b.cₛ)")
  b.δₚ == 0 || print(io, ", δₚ=$(b.δₚ)")
  b.δₛ == 0 || print(io, ", δₛ=$(b.δₛ)")
  b.σ == 0 || print(io, ", σ=$(b.σ)")
  print(io, ")")
end

function reflection_coef(bc::ElasticBoundary, frequency, θ, ρ, c)
  isinf(bc.cₚ) && return 1.0 + 0im
  isinf(bc.cₛ) && return 1.0 + 0im
  bc.σ == 0 || @warn "Roughness is currently ignored in computing reflection coefficient" maxlog=1
  θ = in_units(u"rad", θ)
  ρ = in_units(u"kg/m^3", ρ)
  c = in_units(u"m/s", c)
  reflection_coef(θ, bc.ρ / ρ, bc.cₚ / c, bc.cₛ / c, bc.δₚ, bc.δₛ)
end

"""
    MultilayerElasticBoundary([(h, ρ, cₚ, cₛ), ...])
    MultilayerElasticBoundary([(h, ρ, cₚ, cₛ, δₚ, δₛ), ...])
    MultilayerElasticBoundary([(h, ρ, cₚ, cₛ, δₚ, δₛ, σ), ...])
    MultilayerElasticBoundary([(h, b::FluidBoundary), ...])
    MultilayerElasticBoundary([(h, b::ElasticBoundary), ...])

Create a multilayer solid boundary with layers defined by a vector of tuples
specifying the layer thickness `h`, density `ρ`, compressional sound speed `cₚ`,
shear wave speed `cₛ`, dimensionless compressional absorption coefficient `δₚ`,
dimensionless shear absorption coefficient `δₛ`, and interfacial roughness `σ`
for each layer. If the absorption coefficients or interfacial roughness is
unspecified, they are assumed to be 0. The last entry in the vector is
considered the bottom half-space and has an infinite thickness. By convention,
we specify a value of `Inf` for `h` for that layer.

`ρ`, `cₚ`, and `cₛ` may also be specified with linear variation within a layer
by specifying the value as a 2-tuple, with the first entry being the value at
the top of the layer and the second entry being the value at the bottom.

The properties of a layer may be specified by giving a `FluidBoundary` or
`ElasticBoundary` layer instead of `ρ`, `cₚ`, `cₛ`, `δₚ`, `δₛ` and `σ`.

# Examples
```julia-repl
julia> MultilayerElasticBoundary([
         (h = 5.2, FineSand),
         (h = Inf, ρ = 2000, cₚ = 2500, cₛ = 500)
       ])
MultilayerElasticBoundary(2 layers):
  (h = 5.2, ρ = 1484.373, cₚ = 1691.954, cₛ = 0.0, δₚ = 0.01602, δₛ = 0.0, σ = 0.0)
  (h = Inf, ρ = 2000.0, cₚ = 2500.0, cₛ = 500.0, δₚ = 0.0, δₛ = 0.0, σ = 0.0)

julia> MultilayerElasticBoundary([
         (5.2, 1300, 1700, 100),
         (Inf, 2000, 2500, 500, 0.1, 0.2)
       ])
MultilayerElasticBoundary(2 layers):
  (h = 5.2, ρ = 1300.0, cₚ = 1700.0, cₛ = 100.0, δₚ = 0.0, δₛ = 0.0, σ = 0.0)
  (h = Inf, ρ = 2000.0, cₚ = 2500.0, cₛ = 500.0, δₚ = 0.1, δₛ = 0.2, σ = 0.0)

julia> MultilayerElasticBoundary([
         (5.2, 1300, (1700,2000), 100),
         (Inf, 2000, 2500, 500, 0.1, 0.2)
       ])
MultilayerElasticBoundary(2 layers):
  (h = 5.2, ρ = 1300.0, cₚ = (1700.0, 2000.0), cₛ = 100.0, δₚ = 0.0, δₛ = 0.0, σ = 0.0)
  (h = Inf, ρ = 2000.0, cₚ = 2500.0, cₛ = 500.0, δₚ = 0.1, δₛ = 0.2, σ = 0.0)

julia> MultilayerElasticBoundary([
         (5.2, 1300, 1700, 100, 0, 0, 0.1),
         (Inf, 2000, 2500, 500, 0.1, 0.2, 0.1)
       ])
MultilayerElasticBoundary(2 layers):
  (h = 5.2, ρ = 1300.0, cₚ = 1700.0, cₛ = 100.0, δₚ = 0.0, δₛ = 0.0, σ = 0.1)
  (h = Inf, ρ = 2000.0, cₚ = 2500.0, cₛ = 500.0, δₚ = 0.1, δₛ = 0.2, σ = 0.1)

julia> MultilayerElasticBoundary([
         (h = 5.2, ρ = 1300, cₚ = 1700, cₛ = 100, δₚ = 0.1, δₛ = 0.2, σ = 0.1),
         (h = Inf, ρ = 2000, cₚ = 2500, cₛ = 500)
       ])
MultilayerElasticBoundary(2 layers):
  (h = 5.2, ρ = 1300.0, cₚ = 1700.0, cₛ = 100.0, δₚ = 0.1, δₛ = 0.2, σ = 0.1)
  (h = Inf, ρ = 2000.0, cₚ = 2500.0, cₛ = 500.0, δₚ = 0.0, δₛ = 0.0, σ = 0.0)

julia> MultilayerElasticBoundary([
         (5.2u"m", 1.3u"g/cm^3", 1700u"m/s", 100u"m/s", 0.1, 0.2),
         (Inf, 2u"g/cm^3", 2500u"m/s", 500u"m/s")
       ])
MultilayerElasticBoundary(2 layers):
  (h = 5.2, ρ = 1300.0, cₚ = 1700.0, cₛ = 100.0, δₚ = 0.1, δₛ = 0.2, σ = 0.0)
  (h = Inf, ρ = 2000.0, cₚ = 2500.0, cₛ = 500.0, δₚ = 0.0, δₛ = 0.0, σ = 0.0)
```
"""
struct MultilayerElasticBoundary{T} <: AbstractAcousticBoundary
  layers::Vector{@NamedTuple{h::T, ρ::Union{T,Tuple{T,T}}, cₚ::Union{T,Tuple{T,T}}, cₛ::Union{T,Tuple{T,T}}, δₚ::T, δₛ::T, σ::T}}
  function MultilayerElasticBoundary(layers::AbstractVector)
    isempty(layers) && error("No layers specified")
    layers = map(layers) do l
      if length(l) == 2
        if l[2] isa FluidBoundary
          l = (l[1], l[2].ρ, l[2].c, 0, l[2].δ, 0, l[2].σ)
        elseif l[2] isa ElasticBoundary
          l = (l[1], l[2].ρ, l[2].cₚ, l[2].cₛ, l[2].δₚ, l[2].δₛ, l[2].σ)
        else
          error("Unknown layer type: $(typeof(l[2]))")
        end
      end
      4 ≤ length(l) ≤ 7 || error("Each layer should be defined by a tuple with 4-7 numbers")
      if l isa NamedTuple
        length(l) == 4 && keys(l) != (:h, :ρ, :cₚ, :cₛ) && error("Invalid named tuple definition")
        length(l) == 5 && keys(l) != (:h, :ρ, :cₚ, :cₛ, :δₚ) && error("Invalid named tuple definition")
        length(l) == 6 && keys(l) != (:h, :ρ, :cₚ, :cₛ, :δₚ, :δₛ) && error("Invalid named tuple definition")
        length(l) == 7 && keys(l) != (:h, :ρ, :cₚ, :cₛ, :δₚ, :δₛ, :σ) && error("Invalid named tuple definition")
      end
      l = (float(in_units(u"m", l[1])), in_units.(u"kg/m^3", l[2]),
           in_units.(u"m/s", l[3]), in_units.(u"m/s", l[4]),
           length(l) > 4 ? l[5] : 0, length(l) > 5 ? l[6] : 0,
           length(l) > 6 ? l[7] : 0)
      T = promote_type(map(typeof, Iterators.flatten(l))...)
      l = map(l1 -> T.(l1), l)
      (h=l[1], ρ=l[2], cₚ=l[3], cₛ=l[4], δₚ=l[5], δₛ=l[6], σ=l[7])
    end
    isinf(layers[end][1]) || error("Last layer must have an infinite thickness")
    new{typeof(layers[1][1])}(layers)
  end
end

function Base.show(io::IO, b::MultilayerElasticBoundary)
  print(io, "MultilayerElasticBoundary($(length(b.layers)) layer")
  length(b.layers) > 1 && print(io, "s")
  print(io, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", b::MultilayerElasticBoundary)
  print(io, "MultilayerElasticBoundary($(length(b.layers)) layer")
  length(b.layers) > 1 && print(io, "s")
  println(io, "):")
  for l ∈ b.layers
    println(io, "  (", join(map(k -> "$k = $(repr(l[k]))", keys(l)), ", "), ")")
  end
end

function reflection_coef(bc::MultilayerElasticBoundary, frequency, θ, ρ, c)
  if length(bc.layers) > 1
    # TODO: implement multilayer reflection coefficient
    @warn "Multilayer boundary reflection not implemented, using top layer only..." maxlog=1
  end
  l1 = bc.layers[1]
  bc1 = ElasticBoundary(l1.ρ, l1.cₚ, l1.cₛ, l1.δₚ, l1.δₛ, l1.σ)
  reflection_coef(bc1, frequency, θ, ρ, c)
end

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

function Base.rand(rng::AbstractRNG, noise::WhiteGaussianNoise, nsamples::Integer; fs)
  fs = in_units(u"Hz", fs)
  signal(randn(rng, typeof(noise.σ), nsamples) .* noise.σ, fs)
end

function Base.rand(rng::AbstractRNG, noise::WhiteGaussianNoise, nsamples::Integer, nch::Integer; fs)
  fs = in_units(u"Hz", fs)
  signal(randn(rng, typeof(noise.σ), nsamples, nch) .* noise.σ, fs)
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

function Base.rand(rng::AbstractRNG, noise::RedGaussianNoise, nsamples::Integer; fs)
  fs = in_units(u"Hz", fs)
  signal(rand(rng, RedGaussian(σ=noise.σ, n=nsamples)), fs)
end

function Base.rand(rng::AbstractRNG, noise::RedGaussianNoise, nsamples::Integer, nch::Integer; fs)
  fs = in_units(u"Hz", fs)
  x = Array{typeof(noise.σ), 2}(undef, nsamples, nch)
  dist = RedGaussian(σ=noise.σ, n=nsamples)
  for i ∈ 1:nch
    x[:,i] .= rand(rng, dist)
  end
  signal(x, fs)
end

################################################################################
# fields

abstract type Interpolation end

"""
Linear interpolation.
"""
struct Linear <: Interpolation end

"""
Cubic spline interpolation.
"""
struct CubicSpline <: Interpolation end

struct SampledFieldZ{T1,T2,T3} <: DepthDependent
  f::T1
  zrange::T2
  interp::T3
end

struct SampledFieldX{T1,T2,T3} <: PositionDependent
  f::T1
  xrange::T2
  interp::T3
end

struct SampledFieldXZ{T1,T2,T3} <: PositionDependent
  f::T1
  xrange::T2
  zrange::T2
  interp::T3
end

struct SampledFieldXY{T1,T2,T3} <: PositionDependent
  f::T1
  xrange::T2
  yrange::T2
  interp::T3
end

struct SampledFieldXYZ{T1,T2,T3} <: PositionDependent
  f::T1
  xrange::T2
  yrange::T2
  zrange::T2
  interp::T3
end

(v::SampledFieldZ)(z::Number) = v.f(z)
(v::SampledFieldX)(x::Number) = v.f(x)
(v::SampledFieldXZ)(x, z) = v.f(x, z)
(v::SampledFieldXY)(x, y) = v.f(x, y)
(v::SampledFieldXYZ)(x, y, z) = v.f(x, y, z)

(v::SampledFieldZ)(pos::XYZ) = v.f(pos.z)
(v::SampledFieldX)(pos::XYZ) = v.f(pos.x)
(v::SampledFieldXZ)(pos::XYZ) = v.f(pos.x, pos.z)
(v::SampledFieldXY)(pos::XYZ) = v.f(pos.x, pos.y)
(v::SampledFieldXYZ)(pos::XYZ) = v.f(pos.x, pos.y, pos.z)

Base.show(io::IO, v::SampledFieldZ) = print(io, "SampledField(z-varying, $(length(v.zrange)) samples)")
Base.show(io::IO, v::SampledFieldX) = print(io, "SampledField(x-varying, $(length(v.xrange)) samples)")
Base.show(io::IO, v::SampledFieldXZ) = print(io, "SampledField(xz-varying, $(length(v.xrange))×$(length(v.zrange)) samples)")
Base.show(io::IO, v::SampledFieldXY) = print(io, "SampledField(xy-varying, $(length(v.xrange))×$(length(v.yrange)) samples)")
Base.show(io::IO, v::SampledFieldXYZ) = print(io, "SampledField(xyz-varying, $(length(v.xrange))×$(length(v.yrange))×$(length(v.zrange)) samples)")

Base.minimum(v::SampledFieldZ) = minimum(v.f.(v.zrange))
Base.maximum(v::SampledFieldZ) = maximum(v.f.(v.zrange))
Base.minimum(v::SampledFieldX) = minimum(v.f.(v.xrange))
Base.maximum(v::SampledFieldX) = maximum(v.f.(v.xrange))
Base.minimum(v::SampledFieldXZ) = minimum(v.f(x, z) for x ∈ v.xrange for z ∈ v.zrange)
Base.maximum(v::SampledFieldXZ) = maximum(v.f(x, z) for x ∈ v.xrange for z ∈ v.zrange)
Base.minimum(v::SampledFieldXY) = minimum(v.f(x, y) for x ∈ v.xrange for y ∈ v.yrange)
Base.maximum(v::SampledFieldXY) = maximum(v.f(x, y) for x ∈ v.xrange for y ∈ v.yrange)
Base.minimum(v::SampledFieldXYZ) = minimum(v.f(x, y, z) for x ∈ v.xrange for y ∈ v.yrange for z ∈ v.zrange)
Base.maximum(v::SampledFieldXYZ) = maximum(v.f(x, y, z) for x ∈ v.xrange for y ∈ v.yrange for z ∈ v.zrange)

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

Keyword argument `interp` is used to specify the interpolation method. `Linear()`
interpolation is supported for 1D, 2D and 3D fields. For 2D and 3D fields, the
data must be sampled on a regular grid. For uniformly sampled `1D` fields,
`CubicSpline()` interpolation is also supported.
"""
function SampledField(v; x=nothing, y=nothing, z=nothing, interp=Linear())
  # TODO: support PCHIP interpolation (see HLS-2021-01.pdf in OALIB distribution)
  if x === nothing && y === nothing && z !== nothing
    v = float(v)
    z = float(z)
    ndx = sortperm(z)
    if interp === CubicSpline()
      z isa AbstractRange || error("Cubic interpolation requires `z` to be an `AbstractRange`")
      f = cubic_spline_interpolation(z[ndx], v[ndx]; extrapolation_bc=Flat())
    elseif interp === Linear()
      f = linear_interpolation(z[ndx], v[ndx]; extrapolation_bc=Flat())
    else
      error("Unsupported interpolation")
    end
    SampledFieldZ(f, z, interp)
  elseif x !== nothing && y === nothing && z === nothing
    interp ∈ (Linear(), CubicSpline()) || error("Unsupported interpolation")
    v = float(v)
    x = float(x)
    ndx = sortperm(x)
    if interp === CubicSpline()
      x isa AbstractRange || error("Cubic interpolation requires `x` to be an `AbstractRange`")
      f = cubic_spline_interpolation(x[ndx], v[ndx]; extrapolation_bc=Flat())
    elseif interp === Linear()
      f = linear_interpolation(x[ndx], v[ndx]; extrapolation_bc=Flat())
    else
      error("Unsupported interpolation")
    end
    SampledFieldX(f, x, interp)
  elseif x !== nothing && y === nothing && z !== nothing
    interp === Linear() || error("Unsupported interpolation")
    v = float(v)
    x = float(x)
    z = float(z)
    f = extrapolate(interpolate((x, z), v, Gridded(Interpolations.Linear())), Flat())
    SampledFieldXZ(f, x, z, interp)
  elseif x !== nothing && y !== nothing && z === nothing
    interp === Linear() || error("Unsupported interpolation")
    v = float(v)
    x = float(x)
    y = float(y)
    f = extrapolate(interpolate((x, y), v, Gridded(Interpolations.Linear())), Flat())
    SampledFieldXY(f, x, y, interp)
  elseif x !== nothing && y !== nothing && z !== nothing
    interp === Linear() || error("Unsupported interpolation")
    v = float(v)
    x = float(x)
    y = float(y)
    z = float(z)
    f = extrapolate(interpolate((x, y, z), v, Gridded(Interpolations.Linear())), Flat())
    SampledFieldXYZ(f, x, y, z, interp)
  else
    error("Only x, z, x-z, x-y, or x-y-z interpolation supported")
  end
end

