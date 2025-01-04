# General utilities

```@meta
CurrentModule = UnderwaterAcoustics
```

These utilities are not exported, but are primarily designed for use within the
package, and by related packages that extend `UnderwaterAcoustics.jl`.

## Unit conversion

Support for `Unitful.jl` units:

```@docs
in_units
```

## Positions

Creation of named tuples representing positions:

```@docs
Position
```

## Fields

Fields are quantities that may vary with position. The quantities that do not
vary with position are considered `Constant`. If they vay only with depth, but
not with horizontal position, they are considered `RangeIndependent`. The ones
that vary with any positional coordinate are considered `RangeDependent`.

Scalars are, by definition, `Constant`. Other data types should be tagged with
the appropriate abstract type to help the propagation modeling API check/choose
a model's validity.

```@docs
Constant
RangeIndependent
RangeDependent
is_range_dependent(::Constant)
is_constant
value
```