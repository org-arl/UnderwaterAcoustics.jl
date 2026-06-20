import MAT: matread
import SignalAnalysis: duration, nchannels, SampledSignal, samples, signal
import SignalAnalysis: framerate, nframes, resample, isanalytic, analytic, padded
import Interpolations: interpolate, BSpline, Cubic, Line, OnGrid, scale, extrapolate

export BasebandReplayChannel

struct BasebandReplayChannel{T1,T2} <: AbstractChannelModel
  h::Array{Complex{T1},3}
  θ::Matrix{Float64}        # was Matrix{T1}
  φ::Matrix{Float64}        # was Matrix{T1}
  fs::T1
  fc::Float64               # was T1
  step::Int
  f_resamp::T1
  noise::T2
  function BasebandReplayChannel(h, θ, φ, fs, fc, step::Int=1, f_resamp=1.0; noise=nothing)
    h = ComplexF32.(h)
    θ = Float64.(θ)
    φ = Float64.(φ)
    new{Float32,typeof(noise)}(h, θ, φ, Float32(fs), Float64(fc), step, Float32(f_resamp), noise)
  end
end

function Base.show(io::IO, ch::BasebandReplayChannel)
  print(io, "BasebandReplayChannel($(size(ch.h,2)) × $(round(size(ch.h,3)/ch.fs*ch.step; digits=1)) s, $(ch.fc) Hz, $(ch.fs) Sa/s)")
end

"""
    BasebandReplayChannel(h, θ, fs, fc, step=1; noise=nothing)
    BasebandReplayChannel(h, fs, fc, step=1; noise=nothing)

Construct a baseband replay channel with impulse responses `h` and optional
phase estimates `θ` (theta_hat mode). `fs` is the sampling frequency in Sa/s,
`fc` is the carrier frequency in Hz, and `step` is the decimation rate for the
time axis of `h`. The effective sampling frequency of the impulse responses is
`fs ÷ step` impulse responses per second.

An additive noise model may be optionally specified as `noise`. If specified,
it is used to corrupt the received signals.
"""
function BasebandReplayChannel(h, θ, fs, fc, step::Int=1; noise=nothing)
  fs = in_units(u"Hz", fs)
  fc = in_units(u"Hz", fc)
  φ = Matrix{Float64}(undef, 0, 0)
  BasebandReplayChannel(h, θ, φ, fs, fc, step; noise)
end

function BasebandReplayChannel(h, fs, fc, step::Int=1; noise=nothing)
  fs = in_units(u"Hz", fs)
  fc = in_units(u"Hz", fc)
  θ = Matrix{Float64}(undef, size(h,2), 0)
  φ = Matrix{Float64}(undef, 0, 0)
  BasebandReplayChannel(h, θ, φ, fs, fc, step; noise)
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
    all(["version", "h_hat", "params"] .∈ Ref(keys(data))) || error("Bad channel file format")
    data["version"] == 1.0 || @warn "Unsupported channel file version"
    h = reverse(data["h_hat"]; dims=1)
    M = size(h, 2)
    rxs === (:) && (rxs = 1:M)
    ndims(rxs) == 0 && (rxs = [rxs])
    h = h[:,rxs,:]
    θ = Matrix{Float64}(undef, 0, 0)
    φ = Matrix{Float64}(undef, 0, 0)
    if haskey(data, "phi_hat")
      φ_data = data["phi_hat"]
      size(φ_data, 1) == M || error("Invalid phi_hat size")
      φ = Float64.(transpose(φ_data[rxs,:]))
    elseif haskey(data, "theta_hat")
      θ_data = data["theta_hat"]
      size(θ_data, 1) == M || error("Invalid theta_hat size")
      θ = Float64.(transpose(θ_data[rxs,:]))
    end
    fs = data["params"]["fs_delay"]
    fs_time = data["params"]["fs_time"]
    fc = data["params"]["fc"]
    f_resamp = haskey(data, "f_resamp") ? Float64(only(data["f_resamp"])) : 1.0
    if upsample && fs != fs_time
      step = 1
      h = resample(h, fs / fs_time; dims=3)
    else
      step = round(Int, fs / fs_time)
    end
    # spec: size(phase, 2)/fs_delay == size(h_hat, 3)/fs_time  (phase spans the IR duration)
    let nphase = size(φ, 1) > 0 ? size(φ, 1) : size(θ, 1)
      if nphase > 0
        dur_phase = nphase / fs
        dur_h = size(data["h_hat"], 3) / fs_time
        isapprox(dur_phase, dur_h; rtol=1e-3) ||
          error("Phase/IR duration mismatch: phase spans $(round(dur_phase;digits=3))s, " *
                "h_hat spans $(round(dur_h;digits=3))s (spec requires equal durations)")
      end
    end
    return BasebandReplayChannel(h, θ, φ, fs, fc, step, f_resamp; noise)
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
sampled at the channel's sampling rate. If the channel specifies a passband
resampling factor (`f_resamp`), the output is resampled by that factor to
reproduce the nominal Doppler offset.

If `abstime` is `true`, the returned signals begin at the start of transmission.
Otherwise, the result is relative to the earliest arrival time of the signal
at any receiver. If `noisy` is `true` and the channel has a noise model
associated with it, the received signal is corrupted by additive noise.

If `start` is specified, it specifies the starting time index in the replay channel.
If not specified, a random start time is chosen.
"""
function transmit(ch::BasebandReplayChannel, x; txs=:, rxs=:, abstime=false, noisy=true, fs=nothing, start=nothing)
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
  ȳ = similar(x̄, nframes(x̄) + L - 1, length(rxs))
  h = @view ch.h[:,rxs,start:start+Treq]
  _apply_tvir!(ȳ, x̄, ch.step == 1 ? h : resample(h, ch.step; dims=3))
  if size(ch.φ, 2) > 0
    # phi_hat: apply phase then re-interpolate at time-shifted grid to insert delay drift
    i = (start - 1) * ch.step + 1
    φ_seg = @view(ch.φ[i:i+nframes(ȳ)-1, rxs])
    ȳ .*= cis.(φ_seg)
    t = range(0.0, step=1.0/ch.fs, length=nframes(ȳ))
    for (j, _) ∈ enumerate(rxs)
      drift = Float64.(φ_seg[:, j] ./(2π * ch.fc))
      itp = extrapolate(scale(interpolate(ȳ[:, j], BSpline(Cubic(Line(OnGrid())))), t), Line())
      ȳ[:, j] .= itp.(t .+ drift)
    end
  elseif size(ch.θ, 2) > 0
    # theta_hat: phase only, no delay interpolation
    i = (start - 1) * ch.step + 1
    ȳ .*= cis.(@view(ch.θ[i:i+nframes(ȳ)-1, rxs]))
  end
  # resample to original sampling rate and upconvert to passband
  y = resample(ȳ, fs/ch.fs; dims=1)
  y .*= cispi.(2 * ch.fc * (0:nframes(y)-1) ./ fs)
  # resample in passband to reproduce the nominal Doppler offset, if needed
  isone(ch.f_resamp) || (y = resample(y, Float64(ch.f_resamp); dims=1))
  input_was_analytic || (y = real(y))
  # normalize total mean power across receivers to the number of receivers
  y .*= sqrt(length(rxs) * size(y,1) / sum(abs2, y))
  y = signal(y, fs)
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