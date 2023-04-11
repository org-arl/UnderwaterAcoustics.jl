using RecipesBase
using Colors

@recipe function plot(env::UnderwaterEnvironment;
    sources=[], receivers=[], transmissionloss=[], rays=[], xmargin=0.0, dynamicrange=42.0,
    colors=range(colorant"darkred", colorant"deepskyblue"; length=256)
  )
  if length(transmissionloss) > 0
    size(transmissionloss) == size(receivers) || throw(ArgumentError("Mismatched receivers and transmissionloss"))
    receivers isa AcousticReceiverGrid2D || throw(ArgumentError("Receivers must be an instance of AcousticReceiverGrid2D"))
    minloss = minimum(transmissionloss)
    clims --> (-minloss-dynamicrange, -minloss)
    colorbar --> true
    cguide --> "dB"
    @series begin
      seriestype := :heatmap
      receivers.xrange, receivers.zrange, -transmissionloss'
    end
  end
  xmin = Inf64
  xmax = -Inf64
  length(rays) > 0 && (xmin = min(xmin, minimum(minimum(r1[1] for r1 ∈ r.raypath) for r ∈ rays)))
  length(rays) > 0 && (xmax = max(xmax, maximum(maximum(r1[1] for r1 ∈ r.raypath) for r ∈ rays)))
  length(sources) > 0 && (xmin = min(xmin, minimum(p[1] for p ∈ location.(sources))))
  length(sources) > 0 && (xmax = max(xmax, maximum(p[1] for p ∈ location.(sources))))
  length(receivers) > 0 && (xmin = min(xmin, minimum(p[1] for p ∈ location.(receivers))))
  length(receivers) > 0 && (xmax = max(xmax, maximum(p[1] for p ∈ location.(receivers))))
  isinf(xmin) && (xmin = 0.0)
  isinf(xmax) && (xmax = 1000.0)
  xmin -= xmargin
  xmax += xmargin
  xrange = range(xmin, xmax; length=1000)
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "z (m) "
  z = altimetry(env)
  @series begin
    seriestype := :line
    linecolor := :royalblue
    xrange, [altitude(z, x, 0.0) for x ∈ xrange]
  end
  z = bathymetry(env)
  @series begin
    seriestype := :line
    linecolor := :brown
    xrange, [-depth(z, x, 0.0) for x ∈ xrange]
  end
  if length(rays) > 0
    reverse!(rays)
    ampl = [r.phasor === missing ? -(r.surface + 3*r.bottom) : amp2db(abs.(r.phasor)) for r ∈ rays]
    ampl .-= minimum(ampl)
    if maximum(ampl) > 0.0
      cndx = round.(Int, 1 .+ (length(colors)-1) .* ampl ./ maximum(ampl))
    else
      cndx = ones(Int, length(ampl)) .* length(colors)
    end
    for (j, eigenray) ∈ enumerate(rays)
      r = eigenray.raypath
      @series begin
        seriestype := :line
        linecolor := colors[cndx[j]]
        [r[i][1] for i ∈ 1:length(r)], [r[i][3] for i ∈ 1:length(r)]
      end
    end
  end
  if length(sources) > 0
    @series begin
      seriestype := :scatter
      marker := :star
      seriescolor := :red
      [p[1] for p ∈ location.(sources)], [p[3] for p ∈ location.(sources)]
    end
  end
  if length(receivers) > 0 && length(transmissionloss) == 0
    @series begin
      seriestype := :scatter
      seriescolor := :blue
      [p[1] for p ∈ location.(receivers)], [p[3] for p ∈ location.(receivers)]
    end
  end
end

@recipe function plot(ssp::SoundSpeedProfile; maxdepth=missing, x=0.0)
  D = 10.0
  if maxdepth === missing
    ssp isa SampledSSP && (D = maximum(-ssp.z))
  else
    D = maxdepth
  end
  d = 0.0:-0.1:-D
  c = [soundspeed(ssp, x, 0.0, d1) for d1 ∈ d]
  clim = extrema(c)
  if clim[2] - clim[1] > 50.0
    clim = (round(clim[1]-5.0), round(clim[2]+5.0))
  else
    μ = (clim[2] + clim[1]) / 2.0
    clim = (round(μ - 30.0), round(μ + 30.0))
  end
  ticks --> :native
  xlims --> clim
  xguide --> "soundspeed (m/s)"
  yguide --> "z (m) "
  @series begin
    c, d
  end
end
