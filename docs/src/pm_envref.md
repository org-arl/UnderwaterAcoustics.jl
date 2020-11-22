# Environmental model reference

```@meta
CurrentModule = UnderwaterAcoustics
```

## Environmental model

*Interface*:

- `abstract type UnderwaterEnvironment`
- [`altimetry`](@ref)`()`
- [`bathymetry`](@ref)`()`
- [`ssp`](@ref)`()`
- [`salinity`](@ref)`()`
- [`seasurface`](@ref)`()`
- [`seabed`](@ref)`()`
- [`noise`](@ref)`()`

*Standard models:*

- [`UnderwaterEnvironment`](@ref)

## Sound speed profiles

*Interface*:

- `abstract type SoundSpeedProfile`
- [`soundspeed`](@ref)`()`

*Standard models:*

- [`IsoSSP`](@ref)
- [`MunkSSP`](@ref)
- [`SampledSSP`](@ref)

## Bathymetry

*Interface*:

- `abstract type Bathymetry`
- [`depth`](@ref)`()`
- [`maxdepth`](@ref)`()`

*Standard models:*

- [`ConstantDepth`](@ref)
- [`SampledDepth`](@ref)

## Altimetry

*Interface*:

- `abstract type Altimetry`
- [`altitude`](@ref)`()`

*Standard models:*

- [`FlatSurface`](@ref)
- [`SampledAltitude`](@ref)

## Sea surface

*Interface*:

- `abstract type ReflectionModel`
- [`reflectioncoef`](@ref)`()`

*Standard models:*

- [`ReflectionCoef`](@ref)
- [`RayleighReflectionCoef`](@ref)
- [`SurfaceLoss`](@ref)
- `const Vacuum`
- `const SeaState0`
- `const SeaState1`
- `const SeaState2`
- `const SeaState3`
- `const SeaState4`
- `const SeaState5`
- `const SeaState6`
- `const SeaState7`
- `const SeaState8`
- `const SeaState9`

## Seabed

*Interface*:

- `abstract type ReflectionModel`
- [`reflectioncoef`](@ref)`()`

*Standard models:*

- [`ReflectionCoef`](@ref)
- [`RayleighReflectionCoef`](@ref)
- `const Rock`
- `const Pebbles`
- `const SandyGravel`
- `const CoarseSand`
- `const MediumSand`
- `const FineSand`
- `const VeryFineSand`
- `const ClayeySand`
- `const CoarseSilt`
- `const SandySilt`
- `const Silt`
- `const FineSilt`
- `const SandyClay`
- `const SiltyClay`
- `const Clay`

## Acoustic sources

*Interface*:

- [`AcousticSource`](@ref)
- [`location`](@ref)`()`
- [`nominalfrequency`](@ref)`()`
- [`phasor`](@ref)`()`
- [`record`](@ref)`()`
- [`recorder`](@ref)`()`

*Standard models:*

- [`NarrowbandAcousticSource`](@ref)
- [`Pinger`](@ref)
- [`SampledAcousticSource`](@ref)

# Acoustic receivers

*Interface*:

- [`AcousticReceiver`](@ref)
- [`location`](@ref)`()`

*Standard models:*

- [`AcousticReceiver`](@ref)
- [`AcousticReceiverGrid2D`](@ref)
- [`AcousticReceiverGrid3D`](@ref)

# Ambient noise

*Interface*:

- `abstact type NoiseModel`
- [`record`](@ref)`()`

*Standard models:*

- [`RedGaussianNoise`](@ref)
- any other distribution that works with `rand()`
