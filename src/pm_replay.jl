import MAT: matread
import SignalAnalysis: duration, nchannels, SampledSignal, samples, signal, framerate, nframes, resample, isanalytic, analytic

export BasebandReplayChannel

struct BasebandReplayChannel{T1} <: AbstractChannelModel
  h::Array{Complex{T1},3}   # channel impulse responses (delay x rx × time)
  θ::Matrix{T1}             # phase estimates (rx × time)
  fs::T1                    # sampling frequency (Sa/s)
  fc::T1                    # carrier frequency (Hz)
  step::Int                 # step size for h time axis (fs ÷ step IRs/s)
  function BasebandReplayChannel(h, θ, fs, fc, step=1)
    h = ComplexF32.(h)
    θ = Float32.(θ)
    new{Float32}(h, θ, Float32(fs), Float32(fc), step)
  end
end

function Base.show(io::IO, ch::BasebandReplayChannel)
  print(io, "BasebandReplayChannel($(size(ch.h,2)) × $(round(size(ch.h,3)/ch.fs*ch.step; digits=1)) s, $(ch.fc) Hz, $(ch.fs) Sa/s)")
end

"""
    BasebandReplayChannel(h, θ, fs, fc, step=1)
    BasebandReplayChannel(h, fs, fc, step=1)

Construct a baseband replay channel with impulse responses `h` and phase
estimates `θ`. The phase estimates are optional. `fs` is the sampling
frequency in Sa/s, `fc` is the carrier frequency in Hz, and `step` is the
decimation rate for the time axis of `h`. The effective sampling frequency
of the impulse responses is `fs ÷ step` impulse responses per second.
"""
function BasebandReplayChannel(h, fs, fc, step=1)
  fs = in_units(u"Hz", fs)
  fc = in_units(u"Hz", fc)
  θ = Matrix{Float64}(undef, size(h,2), 0)
  BasebandReplayChannel(h, θ, fs, fc, step)
end

"""
    BasebandReplayChannel(filename; upsample=false, rxs=:)

Load a baseband replay channel from a file.

If `upsample` is `true`, the impulse responses are upsampled to the delay axis
sampling rate. This makes applying the channel faster but requires more memory.
`rxs` controls which receivers to load from the file. By default, all receivers
are loaded.

Supported formats:
- `.mat` (MATLAB) file in underwater acoustic channel repository format.
  See https://github.com/uwa-channels/replay for details.
"""
function BasebandReplayChannel(filename::AbstractString; upsample=false, rxs=:)
  if endswith(filename, ".mat")
    data = matread(filename)
    all(["version", "h_hat", "theta_hat", "params"] .∈ Ref(keys(data))) || error("Bad channel file format")
    data["version"] == 1.0 || @warn "Unsupported channel file version"
    h = data["h_hat"]
    M = size(h, 2)
    rxs === (:) && (rxs = 1:M)
    ndims(rxs) == 0 && (rxs = [rxs])
    h = h[:,rxs,:]
    θ = data["theta_hat"]
    size(θ, 1) == M || error("Invalid phase estimates size")
    θ = θ[rxs,:]
    fs = data["params"]["fs_delay"]
    fc = data["params"]["fc"]
    if upsample && fs != data["params"]["fs_time"]
      step = 1
      h = resample(h, fs / data["params"]["fs_time"]; dims=3)
    else
      step = round(Int, fs / data["params"]["fs_time"])
    end
    return BasebandReplayChannel(h, θ, fs, fc, step)
  else
    error("Unsupported file format")
  end
end

# TODO: support noise models
function transmit(ch::BasebandReplayChannel, x; txs=:, rxs=:, abstime=false, noisy=true, fs=framerate(x), start=nothing)
  D, M, T = size(ch.h)
  maxtime = (T - 1) / ch.fs * ch.step
  txs === (:) && (txs = 1)
  rxs === (:) && (rxs = 1:M)
  ndims(rxs) == 0 && (rxs = [rxs])
  nchannels(x) == 1 || error("Replay channel has only one transmitter")
  length(txs) == 1 || error("Replay channel has only one transmitter")
  only(txs) == 1 || error("Replay channel has only one transmitter")
  abstime && error("Replay channels do not support absolute time")
  all(rx ∈ 1:M for rx ∈ rxs) || error("Invalid receiver indices ($rxs ⊄ 1:$M)")
  fs < 2 * ch.fc && error("Signal sampling rate ($fs Hz) is too low for carrier frequency ($(ch.fc) Hz)")
  input_was_analytic = isanalytic(x)
  x = analytic(signal(samples(x), fs))
  duration(x) ≤ maxtime || error("Signal duration ($(round(duration(x); digits=1)) s) exceeds replay channel duration ($(round(maxtime; digits=1)) s)")
  # convert to baseband and downsample
  x̄ = samples(resample(x .* cispi.(-2 * ch.fc * (0:nframes(x)-1) ./ fs), ch.fs/fs))
  # choose a random start time if not specified
  Treq = ceil(Int, length(x̄) / ch.step) + 1
  start = something(start, rand(1:T-Treq))
  # apply the channel
  ȳ = Array{ComplexF64}(undef, nframes(x̄) + D - 1, length(rxs))
  h = @view ch.h[:,rxs,start:start+Treq]
  _apply_tvir!(ȳ, x̄, ch.step == 1 ? h : resample(h, ch.step; dims=3))
  if size(ch.θ, 2) > 0
    # TODO: apply drift
  end
  # resample to original sampling rate and upconvert to passband
  y = resample(ȳ, fs/ch.fs; dims=1)
  y .*= cispi.(2 * ch.fc * (0:nframes(y)-1) ./ fs)
  signal(input_was_analytic ? y : real(y), fs)
end

# helpers

function _apply_tvir!(y, x, h)
  # TODO: apply time-varying impulse response
end
