using Interpolations
import Unitful: ustrip, Quantity, Units, @u_str

export @u_str, SampledField

################################################################################
# unit conversion utilities, with default units if none are specified

"""
    in_units(u, x)

Get the numerical value of `x` in units of `u`. If `x` is a `Unitful.Quantity`,
it is converted to `u`. If `x` is a number, it is assumed to be in `u` and is
returned as is.

If `u` is `u"dB"`, the `x` may be specified as a number of a `u"dB"` quantity.
In both cases, the numerical value of `x` in dB is returned.

### Examples:
```julia
julia> in_units(u"Hz", 10u"kHz")
10000

julia> in_units(u"m", 10)
10

julia> in_units(u"m", 1u"cm")
1//100

julia> in_units(u"dB", 3)
3

julia> in_units(u"dB", 3u"dB")
3
```
"""
in_units(u::Units, x::Real) = x
in_units(u::Units, x::Quantity) = ustrip(u, x)
in_units(u::typeof(u"dB"), x::Real) = x
in_units(u::typeof(u"dB"), x::typeof(1u"dB")) = ustrip(x)

################################################################################
# locations are represented as named tuples with coordinates in meters

"""
    XYZ(pos)
    XYZ(x, y, z)
    XYZ(x, z)
    XYZ(z)

Convert a location to a named tuple with fields `x`, `y`, and `z`. If any of the
coordinates is not provided, they are assumed to be zero. If the coordinates have
units, they are converted to meters.
"""
XYZ(xyz::NTuple{3,Real}) = NamedTuple{(:x,:y,:z)}(float.(in_units.(u"m", promote(xyz...))))
XYZ(xz::NTuple{2,Real}) = XYZ((xz[1], 0, xz[2]))
XYZ(xyz::NamedTuple{(:x,:y,:z)}) = XYZ(xyz.x, xyz.y, xyz.z)
XYZ(xz::NamedTuple{(:x,:z)}) = XYZ(xz.x, xz.z)
XYZ(z::NamedTuple{(:z,)}) = XYZ(z.z)
XYZ(x::Number, y::Number, z::Number) = XYZ((x, y, z))
XYZ(x::Number, z::Number) = XYZ((x, 0, z))
XYZ(z::Number) = XYZ((0, 0, z))
XYZ(::Nothing) = nothing
XYZ(::Missing) = missing

################################################################################
# fields - types and utilities for quantities that may depend on position
#
# Quantities that are independent of position are represented as scalars.
# Quantities that depend on position are represented as callable structs that
# take a position as argument. Such structs should be tagged with the abstract
# types `DepthDependent` or `PositionDependent` to indicate whether they depend
# only on depth, or on depth and range.

"""
Quantity that may vary with depth, but not with range (`x` or `y` coordinate).
"""
abstract type DepthDependent end

"""
Quantity that may vary with depth and range (`x`, `y` and/or `z` coordinate).
"""
abstract type PositionDependent <: DepthDependent end

"""
    is_range_dependent(q)

Return `true` if the quantity `q` may be range-dependent, `false` if it is
guaranteed to not depend on `x` or `y` coordinate.
"""
is_range_dependent(q::PositionDependent) = true
is_range_dependent(q) = false

"""
    is_constant(q)

Return `true` if the quantity `q` is a constant, `false` if it could depend
on position.
"""
is_constant(q::DepthDependent) = false
is_constant(q) = true

"""
    value(q)
    value(q, pos)

Get the value of the varying quantity `q` at the given position `pos`. `pos`
may be specified as a `(x, y, z)` tuple, a `(x, z)` tuple, a `z` value.
"""
value(q, pos) = q
value(q::DepthDependent, pos) = q(XYZ(pos))
value(q) = value(q, nothing)
value(q::DepthDependent) = error("position not specified")

################################################################################
# interpolants

struct SampledFieldZ{T1,T2} <: DepthDependent
  f::T1
  zrange::T2
end

struct SampledFieldX{T1,T2} <: DepthDependent
  f::T1
  xrange::T2
end

struct SampledFieldXZ{T1,T2} <: PositionDependent
  f::T1
  xrange::T2
  zrange::T2
end

struct SampledFieldXY{T1,T2} <: PositionDependent
  f::T1
  xrange::T2
  yrange::T2
end

struct SampledFieldXYZ{T1,T2} <: PositionDependent
  f::T1
  xrange::T2
  yrange::T2
  zrange::T2
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

Keyword `interp` is used to specify the interpolation method. Only `:linear`
interpolation is supported at present. For 2D and 3D fields, the data must be
sampled on a regular grid.
"""
function SampledField(v; x=nothing, y=nothing, z=nothing, interp=:linear)
  interp === :linear || error("Only linear interpolation supported")
  if x === nothing && y === nothing && z !== nothing
    f = linear_interpolation(z, v; extrapolation_bc=Flat())
    SampledFieldZ(f, z)
  elseif x !== nothing && y === nothing && z === nothing
    f = linear_interpolation(x, v; extrapolation_bc=Flat())
    SampledFieldX(f, x)
  elseif x !== nothing && y === nothing && z !== nothing
    f = extrapolate(interpolate((x, z), v, Gridded(Linear())), Flat())
    SampledFieldXZ(f, x, z)
  elseif x !== nothing && y !== nothing && z === nothing
    f = extrapolate(interpolate((x, y), v, Gridded(Linear())), Flat())
    SampledFieldXY(f, x, y)
  elseif x !== nothing && y !== nothing && z !== nothing
    f = extrapolate(interpolate((x, y, z), v, Gridded(Linear())), Flat())
    SampledFieldXYZ(f, x, y, z)
  else
    error("Only x, z, x-z, x-y, or x-y-z interpolation supported")
  end
end

################################################################################
# general utilities

# fast threaded map, assuming all entries have the same result type
function tmap(f, x)
  x1 = first(x)
  y1 = f(x1)
  y = Array{typeof(y1)}(undef, size(x))
  y[1] = y1
  Threads.@threads for i ∈ eachindex(x)
    i == firstindex(x) && continue
    y[i] = f(x[i])
  end
  y
end
