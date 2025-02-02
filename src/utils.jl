import Unitful: ustrip, Quantity, Units, @u_str
export @u_str

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
XYZ(xyz::NTuple{3,Number}) = NamedTuple{(:x,:y,:z)}(float.(in_units.(u"m", xyz)))
XYZ(xz::NTuple{2,Number}) = XYZ(promote(in_units(u"m",xz[1]), 0, in_units(u"m",xz[2])))
XYZ(xyz::NamedTuple{(:x,:y,:z)}) = XYZ(xyz.x, xyz.y, xyz.z)
XYZ(xz::NamedTuple{(:x,:z)}) = XYZ(xz.x, xz.z)
XYZ(z::NamedTuple{(:z,)}) = XYZ(z.z)
XYZ(x::Number, y::Number, z::Number) = XYZ(promote(in_units.(u"m", (x, y, z))...))
XYZ(x::Number, z::Number) = XYZ(promote(in_units(u"m", x), 0, in_units(u"m", z)))
XYZ(z::Number) = XYZ(promote(0, 0, in_units(u"m", z)))
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
value(q::DepthDependent) = error("Position not specified")

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
