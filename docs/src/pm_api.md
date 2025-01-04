# Propagation modeling API

```@meta
CurrentModule = UnderwaterAcoustics
```

## Propagation modeling

```@docs
AbstractPropagationModel
PekerisRayTracer
AcousticArrival
RayArrival
ModeArrival
transmission_loss
acoustic_field
arrivals
impulse_response
```

## Channel modeling and simulation

```@docs
AbstractChannelModel
channel
transmit
```

## Environment modeling

```@docs
UnderwaterEnvironment
PekerisWaveguide
is_range_dependent(::UnderwaterEnvironment)
isospeed
```

## Boundary conditions

```@docs
AbstractAcousticBoundary
reflection_coef(::AbstractAcousticBoundary, ::Any, ::Any, ::Any, ::Any)
RigidBoundary
PressureReleaseBoundary
FluidBoundary
```

Other fluid seabed boundary conditions based on APL-UW Technical Report 9407:
```julia
UnderwaterAcoustics.Rock — Constant
UnderwaterAcoustics.Pebbles — Constant
UnderwaterAcoustics.SandyGravel — Constant
UnderwaterAcoustics.CoarseSand — Constant
UnderwaterAcoustics.MediumSand — Constant
UnderwaterAcoustics.FineSand — Constant
UnderwaterAcoustics.VeryFineSand — Constant
UnderwaterAcoustics.ClayeySand — Constant
UnderwaterAcoustics.CoarseSilt — Constant
UnderwaterAcoustics.SandySilt — Constant
UnderwaterAcoustics.Silt — Constant
UnderwaterAcoustics.FineSilt — Constant
UnderwaterAcoustics.SandyClay — Constant
UnderwaterAcoustics.SiltyClay — Constant
UnderwaterAcoustics.Clay — Constant
```

```@docs
WindySurface
```

Various surface boundary conditions based on APL-UW Technical Report 9407:
```julia
SeaState0 — Constant
SeaState1 — Constant
SeaState2 — Constant
SeaState3 — Constant
SeaState4 — Constant
SeaState5 — Constant
SeaState6 — Constant
SeaState7 — Constant
SeaState8 — Constant
SeaState9 — Constant
```

## Sources and receivers

```@docs
AbstractAcousticSource
AbstractAcousticReceiver
AcousticSource
AcousticReceiver
AcousticReceiverGrid2D
AcousticReceiverGrid3D
location
frequency
spl
```

## Noise models

```@docs
AbstractNoiseModel
WhiteGaussianNoise
RedGaussianNoise
rand
```
