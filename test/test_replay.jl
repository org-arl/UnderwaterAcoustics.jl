using TestItems

@testsnippet ReplaySetup begin
  using SignalAnalysis
  using StableRNGs

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

  # time-varying impulse response: taps whose gains drift across snapshots
  function make_h_tv(taps)
    h = zeros(ComplexF64, L, M, T)
    for (idx, g) in taps
      for t in 1:T
        h[idx, :, t] .= g * (1 + 0.3 * sin(2π * t / T))   # smooth time variation
      end
    end
    h
  end

  # BPSK-style probe, upsampled and modulated to the carrier. Seeded with a
  # StableRNG so the probe is identical across Julia versions and platforms.
  function make_probe(seed=42)
    rng = StableRNG(seed)
    nsym, rate = 240, 4800.0
    ups = round(Int, FS_IN / rate)
    bb = repeat(Float64.(rand(rng, (-1.0, 1.0), nsym)); inner=ups)
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

  # index of the second strongest arrival, excluding a window around the first
  function second_arrival(mag, p1)
    w = round(Int, 0.0003 * FS_IN)
    m2 = copy(mag)
    m2[max(1, p1+1-w):min(end, p1+1+w)] .= 0
    argmax(m2) - 1
  end
end

@testitem "replay physics" setup=[ReplaySetup] begin
  # two echoes a known number of taps apart must arrive that far apart
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
  # a zero-phase phi channel must reduce to plain convolution. The phi branch
  # still runs the drift interpolator (on a zero drift), so agreement is to
  # interpolation round-off rather than to machine precision.
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
  # a known constant phase φ0 must be recovered from the phi_hat output.
  # tolerance of 0.05 rad accounts for the sub-sample delay drift the phase
  # itself induces (Δτ = φ0/2πfc ≈ 0.22 samples), which biases the estimate by
  # ~φ0·⟨f⟩/fc ≈ 0.02 rad for this probe — a test-design artifact, not a bug.
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
  # a zero-phase theta channel must equal plain convolution exactly: theta is a
  # pure phase multiply with no interpolation, so cis(0) == 1 is exact
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
  # theta is phase-only (no delay drift), so a known constant phase is recovered
  # essentially exactly — far tighter than the phi case, which carries a
  # sub-sample drift bias. The contrast is the point of having both.
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
  # a channel whose IR varies across snapshots exercises the _interp_ir cubic
  # spline path, which constant-in-time channels never meaningfully trigger
  probe = make_probe()
  ch = BasebandReplayChannel(make_h_tv(TAPS), FS_DELAY, FC, STEP)
  y = collect(transmit(ch, signal(probe, FS_IN); start=1, noisy=false))
  @test all(isfinite, y)
  # the two echoes must still land at the correct separation despite the drift
  mag = arrival_mag(y[:, 1], probe)
  p1 = argmax(mag) - 1
  p2 = second_arrival(mag, p1)
  @test abs(p2 - p1) ≈ abs(TAPS[1][1] - TAPS[2][1]) * RATIO atol=3
end

@testitem "replay multi-receiver phases" setup=[ReplaySetup] begin
  # M=3 receivers each with a DIFFERENT constant phi must each recover their own
  # phase — catches per-receiver column-indexing bugs invisible with M=1
  Lm, Mm, Tm = 150, 3, 100
  h = zeros(ComplexF64, Lm, Mm, Tm)
  for (idx, g) in TAPS; h[idx, :, :] .= g; end
  φ0s = [0.3, 0.7, 1.1]
  φ = repeat(reshape(φ0s, 1, Mm), Tm * STEP, 1)   # (time × rx), per-rx phase
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
  # selecting a subset of receivers must return exactly those receivers'
  # signals. The output no longer depends on how many receivers were selected
  # (the per-transmit power normalisation was removed to match the reference),
  # so the subset must match the corresponding column of the full transmit.
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
  # round-trip a channel through the real .mat loader (phi mode), exercising the
  # file reader (format checks, phi selection, duration validation, delay-axis
  # reversal) that the in-memory constructors bypass
  using MAT: matwrite
  Lf, Mf, Tf = 150, 2, 100
  h_file = zeros(ComplexF64, Lf, Mf, Tf)              # [delay, rx, time], UACR layout
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
  @test size(ch_file.h) == (Lf, Mf, Tf)               # loaded with correct shape
  @test size(ch_file.φ, 2) > 0                        # phi mode selected
  @test size(ch_file.θ, 2) == 0                       # theta correctly absent

  # the loaded channel must transmit identically to a direct-constructor channel
  # built from the same data. The loader reverses the delay axis on read, so the
  # direct constructor is given the reversed array to match. Phi is zero here,
  # so both reduce to plain convolution.
  ch_direct = BasebandReplayChannel(reverse(h_file; dims=1), FS_DELAY, FC, STEP)
  x = signal(make_probe(), FS_IN)
  y_file = collect(transmit(ch_file, x; start=1, noisy=false))
  y_direct = collect(transmit(ch_direct, x; start=1, noisy=false))
  @test maximum(abs.(y_file .- y_direct)) / maximum(abs.(y_direct)) < 1e-10

  rm(tmp; force=true)
end

@testitem "replay vs python reference" setup=[ReplaySetup] begin
  # Cross-check the full replayed signal against reference outputs from the
  # Python implementation (github.com/uwa-channels/python). The channel is
  # deterministic and closed-form, so it is rebuilt here bit-for-bit and written
  # to a temporary .mat; only the reference output is committed, as plain text.
  # Regenerate with test/data/gen_references.py, which downloads the reference
  # implementation at a pinned commit.
  #
  # The residual (~2.6e-3 relative, amplitude ratio ~0.9988, the same in every
  # mode) comes from the two resampling steps: DSP.jl's resample and scipy's
  # resample_poly use different anti-aliasing filters, whose gains across the
  # probe band differ by ~0.05-0.07% per step. Given identical input, the
  # channel stage itself matches the reference to ~1e-7, apart from spline end
  # effects at the first snapshot (≤ 6e-4, tv case). The tolerances below
  # allow for that residual.
  using MAT: matwrite

  datadir = joinpath(@__DIR__, "data")
  NL, NM = 16, 2

  # closed-form impulse response, [delay, rx, time]; must match gen_references.py
  function ref_h(nl, nm, nt, step, fd_scale)
    h = zeros(ComplexF64, nl, nm, nt)
    for i ∈ 0:nl-1, j ∈ 0:nm-1, k ∈ 0:nt-1
      g = 0.6^i * cis(π * (i + 3j) / 7)
      fd = fd_scale * sin(1.7i + 0.9j)
      h[i+1, j+1, k+1] = g * cis(2π * fd * k * step / FS_DELAY)
    end
    h
  end

  # closed-form phase vector, [rx, time] at FS_DELAY; must match gen_references.py
  function ref_phase(nm, nphase)
    p = zeros(Float64, nm, nphase)
    for m ∈ 0:nm-1, q ∈ 0:nphase-1
      p[m+1, q+1] = 0.4 * sin(2π * 1.3q / FS_DELAY + 0.7m) + 1.5q / FS_DELAY
    end
    p
  end

  # the probe used to generate the references
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

  # read a reference output (whitespace-separated, '#' comment lines)
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
  # the usable signal length is shorter than the raw channel duration by
  # roughly the impulse response length, and an explicitly supplied start
  # must lie within the range the channel can accommodate
  h = make_h(TAPS)
  ch = BasebandReplayChannel(h, FS_DELAY, FC, STEP)

  # a signal short enough to replay, for the start-index cases
  x = signal(make_probe(), FS_IN)

  # T - Treq is the largest valid start for this signal
  Treq = ceil(Int, (round(Int, nframes(x) * FS_DELAY / FS_IN) + L - 1) / STEP) + 1
  maxstart = T - Treq
  @test maxstart ≥ 1                                    # sanity: probe fits

  # start below range
  @test_throws ErrorException transmit(ch, x; start=0, noisy=false)
  # start above range
  @test_throws ErrorException transmit(ch, x; start=maxstart+1, noisy=false)
  # a valid start still works
  @test size(collect(transmit(ch, x; start=maxstart, noisy=false)), 2) == M

  # a signal too long to replay: fills the whole channel duration, which the
  # old duration-only check allowed but the convolution cannot accommodate
  nlong = round(Int, T * STEP * FS_IN / FS_DELAY)
  xlong = signal(zeros(nlong), FS_IN)
  @test_throws ErrorException transmit(ch, xlong; noisy=false)
end

@testitem "replay storage types" setup=[ReplaySetup] begin
  # sampled data keeps the precision it is given; scalars are always Float64 and
  # must not pass through Float32 (24000.1 would become 24000.099609375)
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
