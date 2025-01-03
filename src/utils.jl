using Unitful: ustrip, Quantity, Units, @u_str

"""
Quantity that does not vary with position.
"""
const Constant = Number

"""
Quantity that may vary with depth, but not with range (`x` or `y` coordinate).
"""
abstract type RangeIndependent end

"""
Quantity that may vary with depth and range (`x`, `y` and/or `z` coordinate).
"""
abstract type RangeDependent end

"""
    is_range_dependent(q)

Return `true` if the quantity `q` may be range-dependent, `false` if it is
guaranteed to not depend on `x` or `y` coordinate.
"""
is_range_dependent(q::Constant) = false
is_range_dependent(q::RangeIndependent) = false
is_range_dependent(q::RangeDependent) = true

"""
    is_constant(q)

Return `true` if the quantity `q` is a constant, `false` if it could depend
on position.
"""
is_constant(q::Constant) = true
is_constant(q::RangeIndependent) = false
is_constant(q::RangeDependent) = false

"""
    value(q)
    value(q, pos)

Get the value of the varying quantity `q` at the given position `pos`. `pos`
may be specified as a `(x, y, z)` tuple, a `(x, z)` tuple, a `z` value.
"""
value(q::Number, pos) = q
value(q, pos) = q(Position(pos))
value(q) = value(q, Position())

"""
    Position(pos)
    Position(x, y, z)
    Position(x, z)
    Position(z)
    Position()

Convert a position to a named tuple with fields `x`, `y`, and `z`. If any of the
coordinates is not provided, they are assumed to be zero.
"""
Position(pos::NTuple{3,Real}) = NamedTuple{(:x,:y,:z)}(in_units.(u"m", promote(pos...)))
Position(pos::NTuple{2,Real}) = Position((pos[1], 0, pos[2]))
Position(pos::NamedTuple{(:x,:y,:z)}) = Position(pos.x, pos.y, pos.z)
Position(pos::NamedTuple{(:x,:z)}) = Position(pos.x, pos.z)
Position(pos::NamedTuple{(:z,)}) = Position(pos.z)
Position(x::Number, y::Number, z::Number) = Position((x, y, z))
Position(x::Number, z::Number) = Position((x, 0, z))
Position(z::Number) = Position((0, 0, z))
Position() = Position((0, 0, 0))

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
