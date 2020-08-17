using Plots
using Colors

@recipe function plot(env::UnderwaterEnvironment;
    sources::Vector{AcousticSource}=[],
    receivers::Vector{AcousticReceiver}=[],
    eigenrays::Vector{Arrival}=[],
    colors=range(colorant"darkred", colorant"deepskyblue"; length=256),
    xmargin=5.0
  )
  xmin = Inf64
  xmax = -Inf64
  length(eigenrays) > 0 && (xmin = min(xmin, minimum(minimum(r1[1] for r1 ∈ r.raypath) for r ∈ eigenrays)))
  length(eigenrays) > 0 && (xmax = max(xmax, maximum(maximum(r1[1] for r1 ∈ r.raypath) for r ∈ eigenrays)))
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
  if length(sources) > 0
    @series begin
      seriestype := :scatter
      marker := :star
      color := :red
      [p[1] for p ∈ location.(sources)], [p[3] for p ∈ location.(sources)]
    end
  end
  if length(receivers) > 0
    @series begin
      seriestype := :scatter
      color := :blue
      [p[1] for p ∈ location.(receivers)], [p[3] for p ∈ location.(receivers)]
    end
  end
  if length(eigenrays) > 0
    ampl = [amp2db(abs.(r.phasor)) for r ∈ eigenrays]
    ampl .-= minimum(ampl)
    if maximum(ampl) > 0.0
      cndx = round.(Int, (length(colors)-1) .* ampl ./ maximum(ampl))
    else
      cndx = zeros(Int, length(ampl))
    end
    for (j, eigenray) ∈ enumerate(eigenrays)
      r = eigenray.raypath
      @series begin
        seriestype := :line
        linecolor := colors[1+cndx[j]]
        [r[i][1] for i ∈ 1:length(r)], [r[i][3] for i ∈ 1:length(r)]
      end
    end
  end
end
