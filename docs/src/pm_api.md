# Propagation modeling API

```@meta
CurrentModule = UnderwaterAcoustics
```

## Propagation modeling

```@docs
AbstractPropagationModel
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
RigidBoundary
PressureReleaseBoundary
FluidBoundary
```

## Sources and receivers

```@docs
AbstractAcousticSource
AbstractAcousticReceiver
NarrowbandAcousticSource
AcousticReceiver
AcousticReceiverGrid2D
AcousticReceiverGrid3D
position
frequency
spl
```
