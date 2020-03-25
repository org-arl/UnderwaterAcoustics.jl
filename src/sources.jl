export generate, ShipNoise, Pinger

function generate(signal::AbstractArray, nsamples, fs)
  length(signal) >= nsamples && return signal[1:nsamples]
  padded(signal, (0,nsamples-length(signal)))
end

Base.@kwdef struct Pinger
  spl = db2amp(180.0)
  frequency = 38000.0
  duration = 0.02
  start = 0.0
  interval = 1.0
  window = nothing
end

function generate(pinger::Pinger, nsamples, fs)
  x = zeros(nsamples)
  ping = real(cw(pinger.frequency, pinger.duration, fs; window=pinger.window).data)
  t = pinger.start
  while t < nsamples/fs
    ndx = round(Int, t*fs) + 1
    n = length(ping)
    ndx+n-1 > nsamples && (n = 1+nsamples-ndx)
    x[ndx:ndx+n-1] .= ping[1:n]
    t += pinger.interval
  end
  pinger.spl * x
end

Base.@kwdef struct ShipNoise
  spl = db2amp(172.0)
  spec = [(50.0, 1.0)]
end

function generate(ship::ShipNoise, nsamples, fs)
  x = rand(PinkGaussian(nsamples))
  y = copy(x)
  for s in ship.spec
    y += s[2]*real(cw(s[1], nsamples/fs, fs).data)[1:nsamples] .* x
  end
  lpf = digitalfilter(Lowpass(3000.0, fs=fs), FIRWindow(hanning(15)))
  y = filtfilt(lpf, y)
  ship.spl * y/rms(y)
end
