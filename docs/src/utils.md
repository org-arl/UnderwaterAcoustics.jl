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

## Locations

Creation of named tuples representing locations:
```@docs
XYZ
```

## Fields

Fields are quantities that may vary with position (e.g. sound speed profile,
bathymetry, etc). The quantities that do not vary with position are considered
constant. If they vay only with depth, but not with horizontal position,
they are considered `DepthDependent`. The ones that vary with any positional
coordinate are considered `PositionDependent`.

Data types that vary with position should be tagged with the appropriate
abstract type to help the propagation modeling API check/choose a model's
validity.

```@docs
DepthDependent
PositionDependent
is_range_dependent(::PositionDependent)
is_constant
value
```

Implementations:
```@docs
SampledField
```

