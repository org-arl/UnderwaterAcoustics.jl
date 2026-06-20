using TestItems

@testsnippet ReplaySetup begin
  using SignalAnalysis

  # number of input samples per delay-rate sample for the params below
  const FS_IN, FC, FS_DELAY, STEP = 96_000.0, 12_000.0, 24_000.0, 20
  const RATIO = FS_IN / FS_DELAY
  const L, M, T = 150, 1, 100        # delay taps, receivers, time snapshots
  const TAPS = [(30, 1.0), (90, 0.7)]

  # (L, M, T) impulse-response cube with constant-in-time taps
  function make_h(taps)
    h = zeros(ComplexF64, L, M, T)
    for (idx, g) in taps
      h[idx, :, :] .= g
    end
    h
  end

  # random BPSK-style probe, upsampled and modulated to the carrier
  function make_probe()
    nsym, rate = 240, 4800.0
    ups = round(Int, FS_IN / rate)
    bb = repeat(Float64.(rand((-1.0, 1.0), nsym)); inner=ups)
    bb .* cos.(2π .* FC .* (0:length(bb)-1) ./ FS_IN)
  end

  # matched-filter magnitude vs lag on a received signal
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
end

@testitem "replay physics" setup=[ReplaySetup] begin
  # two echoes a known number of taps apart must arrive that far apart
  probe = make_probe()
  ch = BasebandReplayChannel(make_h(TAPS), FS_DELAY, FC, STEP)
  y = collect(transmit(ch, signal(probe, FS_IN); start=1, noisy=false))
  mag = arrival_mag(y[:, 1], probe)
  p1 = argmax(mag) - 1
  w = round(Int, 0.0003 * FS_IN)
  m2 = copy(mag); m2[max(1, p1+1-w):min(end, p1+1+w)] .= 0
  p2 = argmax(m2) - 1
  gap_true = abs(TAPS[1][1] - TAPS[2][1]) * RATIO
  @test abs(abs(p2 - p1) - gap_true) ≤ 3
end

@testitem "replay phi=0 identity" setup=[ReplaySetup] begin
  # a zero-phase phi channel must equal plain convolution to machine precision
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
  # a known constant phase must be re-inserted at the correct magnitude
  φ0 = 0.7
  h = make_h(TAPS)
  x = analytic(signal(make_probe(), FS_IN))
  ch_none = BasebandReplayChannel(h, FS_DELAY, FC, STEP)
  ch_phi = BasebandReplayChannel(h, Matrix{Float64}(undef, 0, 0),
                                 fill(φ0, T * STEP, M), FS_DELAY, FC, STEP)
  yn = collect(transmit(ch_none, x; start=1, noisy=false))[:, 1]
  yp = collect(transmit(ch_phi, x; start=1, noisy=false))[:, 1]
  φ_est = angle(sum(yp .* conj(yn)))
  @test abs(rem(φ_est - φ0, 2π, RoundNearest)) < 0.12
end