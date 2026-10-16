using TestItems

@testsnippet ReplaySetup begin
  using SignalAnalysis
  using StableRNGs

  const FS_IN, FC, FS_DELAY, STEP = 96_000.0, 12_000.0, 24_000.0, 20
  const RATIO = FS_IN / FS_DELAY
  const L, M, T = 150, 1, 100
  const TAPS = [(30, 1.0), (90, 0.7)]

  function make_h(taps)
    h = zeros(ComplexF64, L, M, T)
    for (idx, g) in taps
      h[idx, :, :] .= g
    end
    h
  end

  function make_h_tv(taps)
    h = zeros(ComplexF64, L, M, T)
    for (idx, g) in taps
      for t in 1:T
        h[idx, :, t] .= g * (1 + 0.3 * sin(2π * t / T))
      end
    end
    h
  end

  # asymmetric in (i,j), so that a transposed beta changes the covariance
  function make_beta(m, nlags)
    β = zeros(Float64, m, m, nlags)
    for i ∈ 1:m, j ∈ 1:m, k ∈ 1:nlags
      β[i,j,k] = cos(1.3i + 0.7j + 0.5k) / k
    end
    β
  end

  # StableRNG so the probe is identical across Julia versions and platforms.
  function make_probe(seed=42)
    rng = StableRNG(seed)
    nsym, rate = 240, 4800.0
    ups = round(Int, FS_IN / rate)
    bb = repeat(Float64.(rand(rng, (-1.0, 1.0), nsym)); inner=ups)
    bb .* cos.(2π .* FC .* (0:length(bb)-1) ./ FS_IN)
  end

  function arrival_mag(y_m, probe)
    n = length(y_m)
    yv = y_m .* exp.(-im .* 2π .* FC .* (0:n-1) ./ FS_IN)
    xv = probe .* exp.(-im .* 2π .* FC .* (0:length(probe)-1) ./ FS_IN)
    maxlag = ceil(Int, (L + 50) * RATIO)
    mag = zeros(maxlag + 1)
    for lag in 0:maxlag
      s = 0.0im
      for k in (lag+1):min(n, length(xv) + lag)
        s += yv[k] * conj(xv[k-lag])
      end
      mag[lag+1] = abs(s)
    end
    mag
  end

  function second_arrival(mag, p1)
    w = round(Int, 0.0003 * FS_IN)
    m2 = copy(mag)
    m2[max(1, p1+1-w):min(end, p1+1+w)] .= 0
    argmax(m2) - 1
  end
end

@testitem "replay physics" setup=[ReplaySetup] begin
  probe = make_probe()
  ch = BasebandReplayChannel(make_h(TAPS), FS_DELAY, FC, STEP)
  y = collect(transmit(ch, signal(probe, FS_IN); start=1, noisy=false))
  mag = arrival_mag(y[:, 1], probe)
  p1 = argmax(mag) - 1
  p2 = second_arrival(mag, p1)
  gap_true = abs(TAPS[1][1] - TAPS[2][1]) * RATIO
  @test abs(p2 - p1) ≈ gap_true atol=3
end

@testitem "replay phi=0 identity" setup=[ReplaySetup] begin
  # the phi branch still runs the drift interpolator, so agreement is to round-off, not exact
  h = make_h(TAPS)
  x = signal(make_probe(), FS_IN)
  ch_none = BasebandReplayChannel(h, FS_DELAY, FC, STEP)
  ch_phi = BasebandReplayChannel(h, Matrix{Float64}(undef, 0, 0),
                                 zeros(Float64, T * STEP, M), FS_DELAY, FC, STEP)
  y_none = collect(transmit(ch_none, x; start=1, noisy=false))
  y_phi = collect(transmit(ch_phi, x; start=1, noisy=false))
  reldiff = maximum(abs.(y_none .- y_phi)) / maximum(abs.(y_none))
  @test reldiff < 1e-6
end

@testitem "replay constant phase" setup=[ReplaySetup] begin
  # the phase also shifts the delay by φ0/2πfc ≈ 0.22 samples at FS_DELAY, which
  # biases the estimate by φ0·⟨f⟩/fc ≈ 0.02 rad, ⟨f⟩ being the probe's mean
  # baseband frequency; hence the 0.05 rad tolerance
  φ0 = 0.7
  h = make_h(TAPS)
  x = analytic(signal(make_probe(), FS_IN))
  ch_none = BasebandReplayChannel(h, FS_DELAY, FC, STEP)
  ch_phi = BasebandReplayChannel(h, Matrix{Float64}(undef, 0, 0),
                                 fill(φ0, T * STEP, M), FS_DELAY, FC, STEP)
  yn = collect(transmit(ch_none, x; start=1, noisy=false))[:, 1]
  yp = collect(transmit(ch_phi, x; start=1, noisy=false))[:, 1]
  φ_est = angle(sum(yp .* conj(yn)))
  @test abs(rem(φ_est - φ0, 2π, RoundNearest)) < 0.05
end

@testitem "replay theta=0 identity" setup=[ReplaySetup] begin
  h = make_h(TAPS)
  x = signal(make_probe(), FS_IN)
  ch_none = BasebandReplayChannel(h, FS_DELAY, FC, STEP)
  ch_theta = BasebandReplayChannel(h, zeros(Float64, T * STEP, M), FS_DELAY, FC, STEP)
  y_none = collect(transmit(ch_none, x; start=1, noisy=false))
  y_theta = collect(transmit(ch_theta, x; start=1, noisy=false))
  reldiff = maximum(abs.(y_none .- y_theta)) / maximum(abs.(y_none))
  @test reldiff < 1e-10
end

@testitem "replay theta phase" setup=[ReplaySetup] begin
  θ0 = 0.7
  h = make_h(TAPS)
  x = analytic(signal(make_probe(), FS_IN))
  ch_none = BasebandReplayChannel(h, FS_DELAY, FC, STEP)
  ch_theta = BasebandReplayChannel(h, fill(θ0, T * STEP, M), FS_DELAY, FC, STEP)
  yn = collect(transmit(ch_none, x; start=1, noisy=false))[:, 1]
  yt = collect(transmit(ch_theta, x; start=1, noisy=false))[:, 1]
  φ_est = angle(sum(yt .* conj(yn)))
  @test abs(rem(φ_est - θ0, 2π, RoundNearest)) < 1e-6
end

@testitem "replay time-varying h" setup=[ReplaySetup] begin
  # exercises _interp_ir, which constant-in-time channels don't meaningfully test
  probe = make_probe()
  ch = BasebandReplayChannel(make_h_tv(TAPS), FS_DELAY, FC, STEP)
  y = collect(transmit(ch, signal(probe, FS_IN); start=1, noisy=false))
  @test all(isfinite, y)
  mag = arrival_mag(y[:, 1], probe)
  p1 = argmax(mag) - 1
  p2 = second_arrival(mag, p1)
  @test abs(p2 - p1) ≈ abs(TAPS[1][1] - TAPS[2][1]) * RATIO atol=3
end

@testitem "replay multi-receiver phases" setup=[ReplaySetup] begin
  # a different phase per receiver catches column-indexing bugs that M=1 cannot
  Lm, Mm, Tm = 150, 3, 100
  h = zeros(ComplexF64, Lm, Mm, Tm)
  for (idx, g) in TAPS; h[idx, :, :] .= g; end
  φ0s = [0.3, 0.7, 1.1]
  φ = repeat(reshape(φ0s, 1, Mm), Tm * STEP, 1)
  x = analytic(signal(make_probe(), FS_IN))
  ch_none = BasebandReplayChannel(h, FS_DELAY, FC, STEP)
  ch_phi = BasebandReplayChannel(h, Matrix{Float64}(undef, 0, 0), φ, FS_DELAY, FC, STEP)
  yn = collect(transmit(ch_none, x; start=1, noisy=false))
  yp = collect(transmit(ch_phi, x; start=1, noisy=false))
  for m in 1:Mm
    φ_est = angle(sum(yp[:, m] .* conj(yn[:, m])))
    @test abs(rem(φ_est - φ0s[m], 2π, RoundNearest)) < 0.05
  end
end

@testitem "replay receiver subset" setup=[ReplaySetup] begin
  Lm, Mm, Tm = 150, 3, 100
  h = zeros(ComplexF64, Lm, Mm, Tm)
  for (idx, g) in TAPS; h[idx, :, :] .= g; end
  ch = BasebandReplayChannel(h, FS_DELAY, FC, STEP)
  x = signal(make_probe(), FS_IN)
  y_all = collect(transmit(ch, x; start=1, noisy=false))
  y_sub = collect(transmit(ch, x; rxs=[2], start=1, noisy=false))
  @test size(y_sub, 2) == 1
  reldiff = maximum(abs.(y_sub[:, 1] .- y_all[:, 2])) / maximum(abs.(y_all[:, 2]))
  @test reldiff < 1e-10
end

@testitem "replay from file" setup=[ReplaySetup] begin
  using MAT: matwrite
  Lf, Mf, Tf = 150, 2, 100
  h_file = zeros(ComplexF64, Lf, Mf, Tf)
  for (idx, g) in TAPS; h_file[idx, :, :] .= g; end
  phi_file = zeros(Float64, Mf, Tf * STEP)            # [rx, time], length = T*step (spec)

  tmp = joinpath(tempdir(), "uacr_roundtrip_test.mat")
  matwrite(tmp, Dict(
    "version" => 1.0,
    "h_hat" => h_file,
    "phi_hat" => phi_file,
    "params" => Dict("fs_delay" => FS_DELAY, "fs_time" => FS_DELAY / STEP, "fc" => FC),
  ))

  ch_file = BasebandReplayChannel(tmp)
  @test size(ch_file.h) == (Lf, Mf, Tf)
  @test size(ch_file.φ, 2) > 0
  @test size(ch_file.θ, 2) == 0

  # the loader reverses the delay axis, so the direct channel gets the reversed array
  ch_direct = BasebandReplayChannel(reverse(h_file; dims=1), FS_DELAY, FC, STEP)
  x = signal(make_probe(), FS_IN)
  y_file = collect(transmit(ch_file, x; start=1, noisy=false))
  y_direct = collect(transmit(ch_direct, x; start=1, noisy=false))
  @test maximum(abs.(y_file .- y_direct)) / maximum(abs.(y_direct)) < 1e-10

  rm(tmp; force=true)
end

@testitem "replay vs python reference" setup=[ReplaySetup] begin
  # Compares against stored outputs of the Python reference
  # (github.com/uwa-channels/python). The channel is closed-form and rebuilt
  # here; only the outputs are committed. Regenerate them with
  # test/data/gen_references.py.
  #
  # The residual (~2.6e-3 relative, amplitude ratio ~0.9988, the same in every
  # mode) comes from the two resampling steps: DSP.jl's resample and scipy's
  # resample_poly use different anti-aliasing filters, whose gains across the
  # probe band differ by ~0.05-0.07% per step. Given identical input, the
  # channel steps match the reference to float precision, apart from spline end
  # effects (≤ 6e-4, tv case). The tolerances below allow for that residual.
  using MAT: matwrite

  datadir = joinpath(@__DIR__, "data")
  NL, NM = 16, 2

  # must match make_h in gen_references.py
  function ref_h(nl, nm, nt, step, fd_scale)
    h = zeros(ComplexF64, nl, nm, nt)
    for i ∈ 0:nl-1, j ∈ 0:nm-1, k ∈ 0:nt-1
      g = 0.6^i * cis(π * (i + 3j) / 7)
      fd = fd_scale * sin(1.7i + 0.9j)
      h[i+1, j+1, k+1] = g * cis(2π * fd * k * step / FS_DELAY)
    end
    h
  end

  # must match make_phase in gen_references.py
  function ref_phase(nm, nphase)
    p = zeros(Float64, nm, nphase)
    for m ∈ 0:nm-1, q ∈ 0:nphase-1
      p[m+1, q+1] = 0.4 * sin(2π * 1.3q / FS_DELAY + 0.7m) + 1.5q / FS_DELAY
    end
    p
  end

  # must match make_probe in gen_references.py
  function ref_probe(fs; D=0.008, f0=9000.0, f1=15000.0, tau=0.002)
    n = 0:round(Int, D * fs)-1
    t = n ./ fs
    x = cos.(2π .* (f0 .* t .+ 0.5 * (f1 - f0) / D .* t .^ 2))
    w = ones(length(t))
    for (idx, tt) ∈ enumerate(t)
      tt < tau && (w[idx] = 0.5 * (1 - cos(π * tt / tau)))
      tt > D - tau && (w[idx] = 0.5 * (1 - cos(π * (D - tt) / tau)))
    end
    x .* w
  end

  function read_ref(path)
    rows = Vector{Vector{Float64}}()
    for line ∈ eachline(path)
      s = strip(line)
      (isempty(s) || startswith(s, "#")) && continue
      push!(rows, parse.(Float64, split(s)))
    end
    reduce(vcat, transpose.(rows))
  end

  probe = ref_probe(FS_IN)

  #  name    step   T   mode         Doppler (Hz)
  cases = [("none",   1, 240, nothing,     3.0),
           ("theta",  1, 240, "theta_hat", 3.0),
           ("phi",    1, 240, "phi_hat",   3.0),
           ("tv",    20, 120, "phi_hat",  80.0)]

  for (name, step, nt, mode, fd) ∈ cases
    data = Dict{String,Any}(
      "version" => 1.0,
      "h_hat" => ref_h(NL, NM, nt, step, fd),
      "params" => Dict("fs_delay" => FS_DELAY, "fs_time" => FS_DELAY / step, "fc" => FC))
    mode === nothing || (data[mode] = ref_phase(NM, nt * step))

    tmp = joinpath(tempdir(), "uacr_ref_$(name).mat")
    matwrite(tmp, data)
    ch = BasebandReplayChannel(tmp)
    y = collect(transmit(ch, signal(probe, FS_IN); start=1, noisy=false))
    rm(tmp; force=true)

    y_ref = read_ref(joinpath(datadir, "y_$(name).txt"))
    @test size(y, 2) == size(y_ref, 2)
    n = min(size(y, 1), size(y_ref, 1))
    for m ∈ 1:size(y_ref, 2)
      a, b = y[1:n, m], y_ref[1:n, m]
      @test maximum(abs.(a .- b)) / maximum(abs.(b)) < 5e-3
      @test sqrt(sum(abs2, a)) / sqrt(sum(abs2, b)) ≈ 1 atol=0.01
    end
  end
end

@testitem "replay bounds checking" setup=[ReplaySetup] begin
  h = make_h(TAPS)
  ch = BasebandReplayChannel(h, FS_DELAY, FC, STEP)
  x = signal(make_probe(), FS_IN)

  # mirrors the Treq computation in transmit()
  Treq = ceil(Int, (round(Int, nframes(x) * FS_DELAY / FS_IN) + L - 1) / STEP) + 1
  maxstart = T - Treq
  @test maxstart ≥ 1

  @test_throws ErrorException transmit(ch, x; start=0, noisy=false)
  @test_throws ErrorException transmit(ch, x; start=maxstart+1, noisy=false)
  @test size(collect(transmit(ch, x; start=maxstart, noisy=false)), 2) == M

  # fills the raw channel duration, leaving no room for the impulse response tail
  nlong = round(Int, T * STEP * FS_IN / FS_DELAY)
  xlong = signal(zeros(nlong), FS_IN)
  @test_throws ErrorException transmit(ch, xlong; noisy=false)
end

@testitem "replay noise statistics" setup=[ReplaySetup] begin
  using Statistics: cov
  Mn, nlags, fs = 3, 4, 48_000.0
  β = make_beta(Mn, nlags)
  x = samples(rand(StableRNG(11), ReplayNoise(β, fs), 200_000, Mn; fs))
  @test cov(x) ≈ sum(β[:,:,k] * transpose(β[:,:,k]) for k ∈ 1:nlags) rtol=0.05
  # the lag-1 covariance pins the direction of time, which the zero-lag one cannot
  R1 = transpose(x[1:end-1,:]) * x[2:end,:] ./ (size(x, 1) - 1)
  @test R1 ≈ sum(β[:,:,k+1] * transpose(β[:,:,k]) for k ∈ 1:nlags-1) rtol=0.05
end

@testitem "replay noise receivers" setup=[ReplaySetup] begin
  Mn, fs = 4, 48_000.0
  β = make_beta(Mn, 3)
  x_all = samples(rand(StableRNG(3), ReplayNoise(β, fs), 500, Mn; fs))
  x_sub = samples(rand(StableRNG(3), ReplayNoise(β, fs; rxs=[3,1]), 500, 2; fs))
  @test x_sub ≈ x_all[:,[3,1]]
  @test samples(rand(StableRNG(3), ReplayNoise(β, fs; σ=0.5), 500, Mn; fs)) ≈ 0.5 .* x_all
  @test ndims(samples(rand(ReplayNoise(β, fs; rxs=2), 500; fs))) == 1
  @test_throws ErrorException rand(ReplayNoise(β, fs), 500, 2; fs)
  @test_throws ErrorException ReplayNoise(β, fs; rxs=[1,1])
  @test_throws ErrorException ReplayNoise(β, fs; rxs=[Mn+1])
  @test_throws ErrorException ReplayNoise(β, fs, 2.5)
  @test_throws ErrorException ReplayNoise(β, fs, 0.05)
  @test_throws ErrorException ReplayNoise(zeros(2, 3, 4), fs)
end

@testitem "replay noise resampling" setup=[ReplaySetup] begin
  β = make_beta(2, 3)
  noise = ReplayNoise(β, 48_000.0)
  for fs ∈ (24_000.0, 48_000.0, 96_000.0)
    x = rand(StableRNG(5), noise, 1000, 2; fs)
    @test size(x) == (1000, 2)
    @test framerate(x) == fs
    @test all(isfinite, samples(x))
  end
  # lag-1 autocorrelation is 0.48 at the measured rate, 0.83 at 96 kHz and 0.26 at 24 kHz
  r1(x) = sum(x[1:end-1] .* x[2:end]) / sum(abs2, x)
  @test r1(samples(rand(StableRNG(5), noise, 100_000, 2; fs=96_000.0))[:,1]) > 0.7
  @test r1(samples(rand(StableRNG(5), noise, 100_000, 2; fs=24_000.0))[:,1]) < 0.35
  @test eltype(samples(rand(StableRNG(5), ReplayNoise(Float32.(β), 48_000.0), 1000, 2; fs=96_000.0))) === Float32
end

@testitem "replay noise impulsive" setup=[ReplaySetup] begin
  using Statistics: mean
  fs = 48_000.0
  β = zeros(2, 2, 1)
  β[1,1,1] = β[2,2,1] = 1.0
  x = samples(rand(StableRNG(9), ReplayNoise(β, fs, 1.5), 100_000, 2; fs))
  # identity mixing gives SαS noise with scale 1/√2, whose characteristic function is
  # exp(-|t/√2|^α); the estimate has std < 0.002
  for t ∈ (1.0, 2.0)
    @test mean(cos.(t .* x)) ≈ exp(-abs(t / √2)^1.5) atol=0.01
  end
  # with Float32 draws, these 5M would include u = -π/2
  @test all(isfinite, UnderwaterAcoustics._sαsrand(StableRNG(1), Float32, 5_000_000, 1, 1.5))
  # α = 0.1 variates overflow Float32 but not Float64
  @test_throws ErrorException rand(StableRNG(1), ReplayNoise(ones(Float32, 1, 1, 1), fs, 0.1), 100_000; fs)
  @test all(isfinite, samples(rand(StableRNG(1), ReplayNoise(ones(1, 1, 1), fs, 0.1), 100_000; fs)))
end

@testitem "replay noise from file" setup=[ReplaySetup] begin
  using MAT: matwrite
  β = make_beta(3, 4)
  fs = 48_000.0
  dir = mktempdir()
  tmp = joinpath(dir, "noise.mat")
  matwrite(tmp, Dict("version" => 1.0, "Fs" => fs, "alpha" => 1.5, "beta" => β))
  noise = ReplayNoise(tmp)
  @test noise.β == β
  @test noise.fs == fs
  @test noise.α == 1.5
  @test noise.rxs == 1:3
  @test ReplayNoise(tmp; rxs=[2]).rxs == [2]
  @test ReplayNoise(tmp; σ=2u"µPa").β == 2 .* β
  bad = joinpath(dir, "bad.mat")
  matwrite(bad, Dict("version" => 1.0, "Fs" => fs))
  @test_throws ErrorException ReplayNoise(bad)
  @test_throws ErrorException ReplayNoise("noise.txt")
end

@testitem "replay noise in channel" setup=[ReplaySetup] begin
  using Statistics: mean
  # the noise is at twice the signal rate, so transmit() resamples it
  ch = BasebandReplayChannel(make_h(TAPS), FS_DELAY, FC, STEP; noise=ReplayNoise(make_beta(M, 3), 2 * FS_IN))
  x = signal(make_probe(), FS_IN)
  y_clean = collect(transmit(ch, x; start=1, noisy=false))
  y_noisy = collect(transmit(ch, x; start=1))
  @test size(y_noisy) == size(y_clean)
  @test y_noisy != y_clean
  @test all(isfinite, y_noisy)

  # noise powers 1, 100 and 10⁴ identify the receivers
  Mm = 3
  h = zeros(ComplexF64, L, Mm, T)
  for (idx, g) in TAPS; h[idx, :, :] .= g; end
  β = zeros(Mm, Mm, 1)
  for i ∈ 1:Mm; β[i,i,1] = 10.0^(i-1); end
  ch3 = BasebandReplayChannel(h, FS_DELAY, FC, STEP; noise=ReplayNoise(β, FS_IN))
  n2 = collect(transmit(ch3, x; rxs=[2], start=1)) .- collect(transmit(ch3, x; rxs=[2], start=1, noisy=false))
  @test 50 < mean(abs2, n2) < 200
  ch4 = BasebandReplayChannel(h, FS_DELAY, FC, STEP; noise=ReplayNoise(make_beta(4, 3), FS_IN))
  @test_throws ErrorException transmit(ch4, x; start=1)
end

@testitem "replay noise with channel file" setup=[ReplaySetup] begin
  using MAT: matwrite
  using Statistics: mean
  Mf, fs = 3, 48_000.0
  h_file = zeros(ComplexF64, L, Mf, T)
  for (idx, g) in TAPS; h_file[idx, :, :] .= g; end
  tmp = joinpath(mktempdir(), "channel.mat")
  matwrite(tmp, Dict(
    "version" => 1.0,
    "h_hat" => h_file,
    "params" => Dict("fs_delay" => FS_DELAY, "fs_time" => FS_DELAY / STEP, "fc" => FC),
  ))
  # noise powers 1, 100 and 10⁴ identify the receivers; analytic input covers transmit()'s analytic branch
  βd = zeros(Mf, Mf, 1)
  for i ∈ 1:Mf; βd[i,i,1] = 10.0^(i-1); end
  ch = BasebandReplayChannel(tmp; rxs=[3,1], noise=ReplayNoise(βd, FS_IN))
  x = analytic(signal(make_probe(), FS_IN))
  n = collect(transmit(ch, x; start=1)) .- collect(transmit(ch, x; start=1, noisy=false))
  @test 5e3 < mean(abs2, n[:,1]) < 2e4
  @test 0.5 < mean(abs2, n[:,2]) < 2
  n1 = collect(transmit(ch, x; rxs=[2], start=1)) .- collect(transmit(ch, x; rxs=[2], start=1, noisy=false))
  @test 0.5 < mean(abs2, n1) < 2
  β = make_beta(Mf, 3)
  @test BasebandReplayChannel(tmp; rxs=[3,1], noise=ReplayNoise(β, fs; rxs=[3,1])).noise.rxs == [3,1]
  @test_throws ErrorException BasebandReplayChannel(tmp; rxs=[3,1], noise=ReplayNoise(β, fs; rxs=[1,2]))
  @test_throws ErrorException BasebandReplayChannel(tmp; noise=ReplayNoise(make_beta(2, 3), fs))
end

@testitem "replay pink noise" setup=[ReplaySetup] begin
  using Statistics: mean, var, cov, cor
  fs = 48_000.0
  @test sum(abs2, UnderwaterAcoustics._pink_noise_filter()) ≈ 1
  x = samples(rand(StableRNG(13), PinkGaussianNoise(1.0), 500_000, 4; fs))
  # power at low frequencies leaves few independent samples: variances 0.92-1.12, correlations up to 0.02
  @test mean(var(x; dims=1)) ≈ 1 rtol=0.1
  @test maximum(abs(cor(x)[i,j]) for i ∈ 1:4, j ∈ 1:4 if i != j) < 0.1
  p = welch_pgram(x[:,1], 8192; fs, window=hanning)
  f, P = freq(p), power(p)
  m = (f .>= 100) .& (f .<= 10_000)
  lf, lp = log10.(f[m]), 10 .* log10.(P[m])
  @test cov(lf, lp) / var(lf) ≈ -17 rtol=0.15  # measured -16.95
  y = samples(rand(StableRNG(13), PinkGaussianNoise(2.0), 1000, 2; fs))
  @test y ≈ 2 .* samples(rand(StableRNG(13), PinkGaussianNoise(1.0), 1000, 2; fs))
  # short draws keep the full variance (0.99 measured), where a "same" convolution gives 0.80
  z = samples(rand(StableRNG(14), PinkGaussianNoise(1.0), 1000, 2000; fs))
  @test 0.9 < mean(abs2, z) < 1.1
  @test ndims(samples(rand(PinkGaussianNoise(1.0), 1000; fs))) == 1
  @test eltype(samples(rand(PinkGaussianNoise(1f0), 1000, 2; fs))) === Float32
end

@testitem "replay storage types" setup=[ReplaySetup] begin
  # 24000.1 isn't exact in Float32, so this catches fs passing through Float32
  ch = BasebandReplayChannel(make_h(TAPS), zeros(T * STEP, M), 24000.1, 12000.3, STEP)
  @test ch.h isa Array{ComplexF64,3}
  @test ch.θ isa Matrix{Float64}
  @test ch.fs === 24000.1
  @test ch.fc === 12000.3
  @test ch.doppler === 1.0
  ch32 = Float32(ch)
  @test ch32.h isa Array{ComplexF32,3}
  @test ch32.θ isa Matrix{Float32}
  @test ch32.fs === 24000.1
  @test BasebandReplayChannel(ComplexF32.(make_h(TAPS)), FS_DELAY, FC, STEP).h isa Array{ComplexF32,3}
end
