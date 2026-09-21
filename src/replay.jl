import SignalAnalysis: duration, nchannels, SampledSignal, samples, signal
import SignalAnalysis: framerate, nframes, resample, isanalytic, analytic, padded
import Interpolations: interpolate, BSpline, Cubic, Line, OnGrid, scale, extrapolate

export BasebandReplayChannel

# fs, fc and doppler stay Float64 even when T1 is Float32, since errors in them accumulate over time
struct BasebandReplayChannel{T1,T2} <: AbstractChannelModel
  h::Array{Complex{T1},3}   # channel impulse responses (delay × rx × time)
  θ::Matrix{T1}             # theta_hat phase estimates (time × rx), or 0×0 if unused
  φ::Matrix{T1}             # phi_hat phase estimates (time × rx), or 0×0 if unused
  fs::Float64               # delay-axis sampling rate, fs_delay (Sa/s)
  fc::Float64               # carrier frequency (Hz)
  step::Int                 # step size for h time axis (fs ÷ step IRs/s)
  doppler::Float64          # passband resampling factor (f_resamp)
  noise::T2
  function BasebandReplayChannel(h, θ::AbstractMatrix, φ::AbstractMatrix, fs::Number, fc::Number, step::Int=1, doppler::Real=1.0; noise=nothing)
    fs = in_units(u"Hz", fs)
    fc = in_units(u"Hz", fc)
    T1 = float(real(eltype(h)))
    new{T1,typeof(noise)}(Complex{T1}.(h), T1.(θ), T1.(φ), Float64(fs), Float64(fc), step, Float64(doppler), noise)
  end
end

"""
    Float32(ch::BasebandReplayChannel)
    Float64(ch::BasebandReplayChannel)

Convert the impulse responses and phase estimates of a replay channel to the
given precision. `Float32(ch)` halves the memory used by a channel loaded from
a file.
"""
(::Type{T})(ch::BasebandReplayChannel) where {T<:AbstractFloat} =
  BasebandReplayChannel(Complex{T}.(ch.h), ch.θ, ch.φ, ch.fs, ch.fc, ch.step, ch.doppler; noise=ch.noise)

function Base.show(io::IO, ch::BasebandReplayChannel)
  print(io, "BasebandReplayChannel($(size(ch.h,2)) × $(round(size(ch.h,3)/ch.fs*ch.step; digits=1)) s, $(ch.fc) Hz, $(ch.fs) Sa/s)")
end

"""
    BasebandReplayChannel(h, θ, φ, fs, fc, step=1, doppler=1.0; noise=nothing)
    BasebandReplayChannel(h, θ, fs, fc, step=1; noise=nothing)
    BasebandReplayChannel(h, fs, fc, step=1; noise=nothing)

Construct a baseband replay channel with impulse responses `h` and optional
phase estimates `θ` (theta_hat, phase tracking only) or `φ` (phi_hat, delay
tracking). `fs` is the sampling frequency in Sa/s, `fc` is the carrier frequency
in Hz, and `step` is the decimation rate for the time axis of `h`. The effective
sampling frequency of the impulse responses is `fs ÷ step` impulse responses per
second. `doppler` is a time-invariant passband resampling factor. The impulse
responses and phase estimates are stored at the precision of `h`.

To use `φ` without `θ`, pass `zeros(0, 0)` for `θ`. If both are given, `φ`
takes precedence.

An additive noise model may be optionally specified as `noise`. If specified,
it is used to corrupt the received signals.
"""
function BasebandReplayChannel(h, θ::AbstractMatrix, fs::Number, fc::Number, step::Int=1; noise=nothing)
  φ = Matrix{Float64}(undef, 0, 0)
  BasebandReplayChannel(h, θ, φ, fs, fc, step; noise)
end

function BasebandReplayChannel(h, fs::Number, fc::Number, step::Int=1; noise=nothing)
  θ = Matrix{Float64}(undef, 0, 0)
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
  See https://github.com/uwa-channels/ for details. Loading `.mat` files
  requires the `MAT` package to be loaded (`using MAT`).
"""
function BasebandReplayChannel(filename::AbstractString; upsample=false, rxs=:, noise=nothing)
  endswith(filename, ".mat") || error("Unsupported file format")
  applicable(_load_mat_replay_channel, filename, upsample, rxs, noise) ||
    error("Loading .mat replay channels requires the MAT package; run `using MAT` first")
  _load_mat_replay_channel(filename, upsample, rxs, noise)
end

# implemented in MATExt
function _load_mat_replay_channel end

"""
    transmit(ch::BasebandReplayChannel, x; txs=:, rxs=:, abstime=false, noisy=true, fs=nothing, start=nothing)

Simulate the transmission of passband signal `x` through the channel model `ch`.
If `txs` is specified, it specifies the indices of the sources active in the
simulation. The number of sources must match the number of channels in the
input signal. If `rxs` is specified, it specifies the indices of the
receivers active in the simulation. Returns the received signal at the
specified (or all) receivers.

`fs` specifies the sampling rate of the input signal. The output signal is
sampled at the same rate. If `fs` is not specified but `x` is a `SampledSignal`,
the sampling rate of `x` is used; otherwise an error is raised. If the channel
has a passband resampling factor (`doppler`), the output is also resampled by
that factor to reproduce the nominal Doppler offset.

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
  maxtime = ((T - 2) * ch.step - L + 1) / ch.fs
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
  x̄ = samples(resample(x .* cispi.(-2 * ch.fc * (0:nframes(x)-1) ./ fs), ch.fs/fs))
  Treq = ceil(Int, (nframes(x̄) + L - 1) / ch.step) + 1
  Treq < T || error("Signal duration ($(round(duration(x); digits=3)) s) exceeds maximum replayable duration ($(floor(maxtime; digits=3)) s)")
  start = something(start, rand(1:T-Treq))
  1 ≤ start ≤ T - Treq || error("Invalid start index ($start ∉ 1:$(T-Treq))")
  ȳ = similar(x̄, nframes(x̄) + L - 1, length(rxs))
  # only the spline in _interp_ir needs the pad; step == 1 must use the exact window
  pad = ch.step == 1 ? 0 : 2
  lo = max(1, start - pad)
  hi = min(T, start + Treq + pad)
  h = @view ch.h[:,rxs,lo:hi]
  _apply_tvir!(ȳ, x̄, ch.step == 1 ? h : _interp_ir(h, ch.step, nframes(ȳ), (start - lo) * ch.step))
  if size(ch.φ, 2) > 0
    i = (start - 1) * ch.step + 1
    φ_seg = @view(ch.φ[i:i+nframes(ȳ)-1, rxs])
    ȳ .*= cis.(φ_seg)
    t = range(0.0, step=1.0/ch.fs, length=nframes(ȳ))
    for (j, _) ∈ enumerate(rxs)
      drift = φ_seg[:, j] ./(2π * ch.fc)
      itp = extrapolate(scale(interpolate(@view(ȳ[:, j]), BSpline(Cubic(Line(OnGrid())))), t), 0.0)
      ȳ[:, j] .= itp.(t .+ drift)
    end
  elseif size(ch.θ, 2) > 0
    # theta_hat: phase only; the delay drift is already in h
    i = (start - 1) * ch.step + 1
    ȳ .*= cis.(@view(ch.θ[i:i+nframes(ȳ)-1, rxs]))
  end
  y = resample(ȳ, fs/ch.fs; dims=1)
  y .*= cispi.(2 * ch.fc * (0:nframes(y)-1) ./ fs)
  isone(ch.doppler) || (y = resample(y, ch.doppler; dims=1))
  input_was_analytic || (y = real(y) .* √2) # undo analytic()'s 1/√2 scaling
  y = signal(y, fs)
  if noisy && ch.noise !== nothing
    if input_was_analytic
      y .+= analytic(rand(ch.noise, size(y); fs))
    else
      y .+= rand(ch.noise, size(y); fs)
    end
  end
  y
end

function _apply_tvir!(y, x, h)
  L = size(h, 1)
  x = padded(x, L - 1)
  for i ∈ 1:size(y,1)
    y[i,:] .= @views transpose(h[:,:,i]) * x[i-L+1:i]
  end
  y
end

function _interp_ir(h, step, n, offset=0)
  L, M, T = size(h)
  out = similar(h, L, M, n)
  ts = range(0.0, step=float(step), length=T)
  for m ∈ 1:M, l ∈ 1:L
    itp = extrapolate(scale(interpolate(@view(h[l, m, :]), BSpline(Cubic(Line(OnGrid())))), ts), 0.0)
    for i ∈ 1:n
      out[l, m, i] = itp(float(i - 1 + offset))
    end
  end
  out
end
