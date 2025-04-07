# Developing your own propagation or channel models

```@meta
CurrentModule = UnderwaterAcoustics
```

`UnderwaterAcoustics.jl` is designed to allow the community to extend it by adding propagation models and channel models.

## Propagation models

A new propagation model (let's call it `MyPropagationModel` for illustration) should define a type that extends one of:
```@docs; canonical=false
AbstractPropagationModel
AbstractRayPropagationModel
AbstractModePropagationModel
```

The constructor for `MyPropagationModel` usually will take in an environmental description and optionally, keyword options that control the model:
```julia
MyPropagationModel(env::UnderwaterEnvironment; kwargs...)
```
However, for data-driven models, the constructor might take in data or partial environmental information and data:
```julia
MyPropagationModel(data; kwargs...)
MyPropagationModel(env, data; kwargs...)
```

The following methods should be defined for `MyPropagationModel`:
```@docs; canonical=false
acoustic_field
```
The acoustic field is represented by complex numbers with amplitude that is related to the source level (`spl`) and transmission loss, and angle that is related to the acoustic phase at the source frequency.

!!! tip "Additional options"
    The `acoustic_field()`, `transmission_loss()`, `arrivals()` and `impulse_response()` methods may support propagation model specific keyword arguments (options) to control finer details of the propagation model.

A `transmission_loss()` method may be optionally defined for `MyPropagationModel`:
```@docs; canonical=false
transmission_loss
```
If it is not defined, the transmission loss is automatically computed as `20 * log10(acoustic_field(...))`.

A propagation model should also typically define:
```@docs; canonical=false
arrivals
```
The returned arrivals should be an array of arrivals that extend:
```@docs; canonical=false
AbstractAcousticArrival
```
The information held in an arrival is propagation model dependent. For example, ray models may return arrivals that contain ray information such as time of arrival, angle of arrival, amplitude and phase of arrival, eigenpath, etc. On the other hand, models based on normal modes may return arrivals containing mode information such as mode number, horizontal and vertical wavenumber, etc. For ray and mode arrivals, the following concrete subtypes should be used when possible:
```@docs; canonical=false
RayArrival
ModeArrival
```

Some propagation models may be able to generate an impulse response. If so, they should define:
```@docs; canonical=false
impulse_response
```
If defined, the impulse response may be used to generate a channel model automatically (by calling `channel()`).

## Channel model

If a propagation model can estimate a received signal from a transmit signal without having to compute an impulse response and convolve it, it may wish to implement the channel modeling API directly:
```@docs; canonical=false
channel
```
The returned channel model must extend:
```@docs; canonical=false
AbstractChannelModel
```
and support:
```@docs; canonical=false
transmit
```

In some cases (e.g. channel replay techniques), a channel model may be defined without the need to derive it from a propagation model. In such a case, one may extend `UnderwaterAcoustics.jl` by directly defining the channel model.
