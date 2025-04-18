import Unitful: ustrip, Quantity, Units, @u_str
export @u_str, °, XYZ, is_constant, is_range_dependent, value

################################################################################
# unit conversion utilities, with default units if none are specified

# define degree symbol
const ° = u"°"

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
# positions are represented as named tuples with coordinates in meters

# position type
const XYZ = NamedTuple{(:x,:y,:z)}

"""
    xyz(pos)
    xyz(x, y, z)
    xyz(x, z)
    xyz(z)

Convert a position to a named tuple with fields `x`, `y`, and `z`. If any of the
coordinates is not provided, they are assumed to be zero. If the coordinates have
units, they are converted to meters.
"""
xyz(pos::NTuple{3,Number}) = NamedTuple{(:x,:y,:z)}(float.(in_units.(u"m", pos)))
xyz(xz::NTuple{2,Number}) = xyz(promote(in_units(u"m",xz[1]), 0, in_units(u"m",xz[2])))
xyz(pos::NamedTuple{(:x,:y,:z)}) = xyz(pos.x, pos.y, pos.z)
xyz(xz::NamedTuple{(:x,:z)}) = xyz(xz.x, xz.z)
xyz(z::NamedTuple{(:z,)}) = xyz(z.z)
xyz(x::Number, y::Number, z::Number) = xyz(promote(in_units.(u"m", (x, y, z))...))
xyz(x::Number, z::Number) = xyz(promote(in_units(u"m", x), 0, in_units(u"m", z)))
xyz(z::Number) = xyz(promote(0, 0, in_units(u"m", z)))
xyz(::Nothing) = nothing
xyz(::Missing) = missing

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
may be specified as a `(x, y, z)` tuple, a `(x, z)` tuple, or a `z` value.

# Examples
```julia
value(q)             # get value of a constant quantity
value(q, -10)        # get value of a depth-dependent quantity at z=-10
value(q, (1000,-10)) # get value of a position-dependent quantity at x=1000, z=-10
value(q, (0,0,-10))  # get value of a position-dependent quantity at (0,0,-10)
```
"""
value(q, pos) = q
value(q::DepthDependent, pos) = q(xyz(pos))
value(q) = value(q, nothing)
value(q::DepthDependent) = error("Position not specified")

"""
    minimum(q)

Get the minimum value of a field quantity `q`.
"""
Base.minimum

"""
    maximum(q)

Get the maximum value of a field quantity `q`.
"""
Base.maximum

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
