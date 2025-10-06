import MAT: matread
import SignalAnalysis: duration, nchannels, SampledSignal, samples, signal
import SignalAnalysis: framerate, nframes, resample, isanalytic, analytic, padded

export BasebandReplayChannel

struct BasebandReplayChannel{T1,T2} <: AbstractChannelModel
  h::Array{Complex{T1},3}   # channel impulse responses (delay × rx × time)
  θ::Matrix{T1}             # phase estimates (time × rx)
  fs::T1                    # sampling frequency (Sa/s)
  fc::T1                    # carrier frequency (Hz)
  step::Int                 # step size for h time axis (fs ÷ step IRs/s)
  noise::T2                 # noise model
  function BasebandReplayChannel(h, θ, fs, fc, step::Int=1; noise=nothing)
    h = ComplexF32.(h)
    θ = Float32.(θ)
    new{Float32,typeof(noise)}(h, θ, Float32(fs), Float32(fc), step, noise)
  end
end

function Base.show(io::IO, ch::BasebandReplayChannel)
  print(io, "BasebandReplayChannel($(size(ch.h,2)) × $(round(size(ch.h,3)/ch.fs*ch.step; digits=1)) s, $(ch.fc) Hz, $(ch.fs) Sa/s)")
end

"""
    BasebandReplayChannel(h, θ, fs, fc, step=1; noise=nothing)
    BasebandReplayChannel(h, fs, fc, step=1; noise=nothing)

Construct a baseband replay channel with impulse responses `h` and phase
estimates `θ`. The phase estimates are optional. `fs` is the sampling
frequency in Sa/s, `fc` is the carrier frequency in Hz, and `step` is the
decimation rate for the time axis of `h`. The effective sampling frequency
of the impulse responses is `fs ÷ step` impulse responses per second.

An additive noise model may be optionally specified as `noise`. If specified,
it is used to corrupt the received signals.
"""
function BasebandReplayChannel(h, fs, fc, step::Int=1; noise=nothing)
  fs = in_units(u"Hz", fs)
  fc = in_units(u"Hz", fc)
  θ = Matrix{Float64}(undef, size(h,2), 0)
  BasebandReplayChannel(h, θ, fs, fc, step; noise)
end

"""
    BasebandReplayChannel(filename; upsample=false, rxs=:, noise=nothing)

Load a baseband replay channel from a file.

If `upsample` is `true`, the impulse responses are upsampled to the delay axis
sampling rate. This makes applying the channel faster but requires more memory.
`rxs` controls which receivers to load from the file. By default, all receivers
are loaded.

An additive noise model may be optionally specified as `noise`. If specified,
it is used to corrupt the received signals.

Supported formats:
- `.mat` (MATLAB) file in underwater acoustic channel repository (UACR) format.
  See https://github.com/uwa-channels/ for details.
"""
function BasebandReplayChannel(filename::AbstractString; upsample=false, rxs=:, noise=nothing)
  # TODO: support UACR noise models
  if endswith(filename, ".mat")
    data = matread(filename)
    all(["version", "h_hat", "theta_hat", "params"] .∈ Ref(keys(data))) || error("Bad channel file format")
    data["version"] == 1.0 || @warn "Unsupported channel file version"
    h = reverse(data["h_hat"]; dims=1)
    M = size(h, 2)
    rxs === (:) && (rxs = 1:M)
    ndims(rxs) == 0 && (rxs = [rxs])
    h = h[:,rxs,:]
    θ = data["theta_hat"]
    size(θ, 1) == M || error("Invalid phase estimates size")
    θ = transpose(θ[rxs,:])
    fs = data["params"]["fs_delay"]
    fc = data["params"]["fc"]
    if upsample && fs != data["params"]["fs_time"]
      step = 1
      h = resample(h, fs / data["params"]["fs_time"]; dims=3)
    else
      step = round(Int, fs / data["params"]["fs_time"])
    end
    return BasebandReplayChannel(h, θ, fs, fc, step; noise)
  else
    error("Unsupported file format")
  end
end

"""
    transmit(ch::BasebandReplayChannel, x; rxs=:, abstime=false, noisy=true, fs=nothing, start=nothing)

Simulate the transmission of passband signal `x` through the channel model `ch`.
If `txs` is specified, it specifies the indices of the sources active in the
simulation. The number of sources must match the number of channels in the
input signal. If `rxs` is specified, it specifies the indices of the
receivers active in the simulation. Returns the received signal at the
specified (or all) receivers.

`fs` specifies the sampling rate of the input signal. The output signal is
sampled at the same rate. If `fs` is not specified but `x` is a `SampledSignal`,
the sampling rate of `x` is used. Otherwise, the signal is assumed to be
sampled at the channel's sampling rate.

If `abstime` is `true`, the returned signals begin at the start of transmission.
Otherwise, the result is relative to the earliest arrival time of the signal
at any receiver. If `noisy` is `true` and the channel has a noise model
associated with it, the received signal is corrupted by additive noise.

If `start` is specified, it specifies the starting time index in the replay channel.
If not specified, a random start time is chosen.
"""
function transmit(ch::BasebandReplayChannel, x; txs=:, rxs=:, abstime=false, noisy=true, fs=nothing, start=nothing)
  # TODO: add test cases for BasebandReplayChannel
  fs === nothing && x isa SampledSignal && (fs = framerate(x))
  L, M, T = size(ch.h)
  maxtime = (T - 1) / ch.fs * ch.step
  txs === (:) && (txs = 1)
  rxs === (:) && (rxs = 1:M)
  ndims(rxs) == 0 && (rxs = [rxs])
  nchannels(x) == 1 || error("Replay channel has only one transmitter")
  length(txs) == 1 || error("Replay channel has only one transmitter")
  only(txs) == 1 || error("Replay channel has only one transmitter")
  abstime && error("Replay channels do not support absolute time")
  all(rx ∈ 1:M for rx ∈ rxs) || error("Invalid receiver indices ($rxs ⊄ 1:$M)")
  fs === nothing && error("Sampling rate must be specified")
  fs < 2 * ch.fc && error("Signal sampling rate ($fs Hz) is too low for carrier frequency ($(ch.fc) Hz)")
  input_was_analytic = isanalytic(x)
  x = analytic(signal(samples(x), fs))
  duration(x) ≤ maxtime || error("Signal duration ($(round(duration(x); digits=1)) s) exceeds replay channel duration ($(round(maxtime; digits=1)) s)")
  # convert to baseband and downsample
  x̄ = samples(resample(x .* cispi.(-2 * ch.fc * (0:nframes(x)-1) ./ fs), ch.fs/fs))
  # choose a random start time if not specified
  Treq = ceil(Int, (nframes(x̄) + L - 1) / ch.step) + 1
  start = something(start, rand(1:T-Treq))
  # apply the channel
  ȳ = similar(x̄, nframes(x̄) + L - 1, length(rxs))
  h = @view ch.h[:,rxs,start:start+Treq]
  _apply_tvir!(ȳ, x̄, ch.step == 1 ? h : resample(h, ch.step; dims=3))
  if size(ch.θ, 2) > 0
    # adjust phase
    i = start * ch.step
    ȳ .*= cis.(@view(ch.θ[i:i+nframes(ȳ)-1,rxs]))
    # TODO: apply drift
  end
  # resample to original sampling rate and upconvert to passband
  y = resample(ȳ, fs/ch.fs; dims=1)
  y .*= cispi.(2 * ch.fc * (0:nframes(y)-1) ./ fs)
  y = signal(input_was_analytic ? y : real(y), fs)
  # add noise
  if noisy && ch.noise !== nothing
    if input_was_analytic
      y .+= analytic(rand(ch.noise, size(y); fs))
    else
      y .+= rand(ch.noise, size(y); fs)
    end
  end
  y
end

# helpers

function _apply_tvir!(y, x, h)
  L = size(h, 1)
  x = padded(x, L - 1)
  for i ∈ 1:size(y,1)
    y[i,:] .= @views transpose(h[:,:,i]) * x[i-L+1:i]
  end
  y
end
