module UnderwaterAcousticsPlots

using RecipesBase
using Colors
using UnderwaterAcoustics
import UnderwaterAcoustics: AbstractAcousticSource, AbstractAcousticReceiver, RayArrival, value
import DSP: amp2db

@recipe function plot(env::UnderwaterEnvironment)
  z1 = env.altimetry
  z2 = env.bathymetry
  xmin, xmax = get(plotattributes, :xlims, (0.0, 10*value(z2, 0)))
  xrange = range(xmin, xmax; length=1000)
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "z (m)"
  @series begin
    seriestype := :line
    linecolor := :royalblue
    xrange, [value(z1, (x, 0.0)) for x ∈ xrange]
  end
  @series begin
    seriestype := :line
    linecolor := :brown
    xrange, [-value(z2, (x, 0.0)) for x ∈ xrange]
  end
end

@recipe function plot(tx::Union{AbstractAcousticSource,AbstractVector{<:AbstractAcousticSource}})
  tx isa AbstractVector || (tx = [tx])
  pos = location.(vec(tx))
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "z (m)"
  @series begin
    seriestype := :scatter
    marker := :star
    seriescolor := :red
    [p.x for p ∈ pos], [p.z for p ∈ pos]
  end
end

@recipe function plot(rx::Union{AbstractAcousticReceiver,AbstractVector{<:AbstractAcousticReceiver}})
  rx isa AbstractVector || (rx = [rx])
  pos = location.(vec(rx))
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "z (m)"
  @series begin
    seriestype := :scatter
    markersize := 3
    seriescolor := :blue
    [p.x for p ∈ pos], [p.z for p ∈ pos]
  end
end

@recipe function plot(rx::Union{AcousticReceiverGrid2D,AcousticReceiverGrid3D})
  pos = location.(vec(rx))
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "z (m)"
  @series begin
    seriestype := :scatter
    markersize := 3
    seriescolor := :blue
    [p.x for p ∈ pos], [p.z for p ∈ pos]
  end
end

@recipe function plot(rx::AcousticReceiverGrid2D, x::AbstractMatrix; crange=42.0)
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "z (m)"
  colorbar_title --> "Transmission loss (dB)"
  cmin = minimum(x)
  clims --> (cmin, cmin + crange)
  colorbar --> true
  color --> :YlGnBu
  cguide --> "dB"
  @series begin
    seriestype := :heatmap
    rx.xrange, rx.zrange, x'
  end
end

@recipe function plot(rays::AbstractVector{<:RayArrival};
  colors=range(colorant"lightgray", colorant"darkblue"; length=256)
)
  rays = filter(r -> r.path isa AbstractVector, rays)
  sort!(rays; by=r->r.t, rev=true)
  ampl = [r.ϕ == 0 ? -(r.surface + 3*r.bottom) : amp2db(abs.(r.ϕ)) for r ∈ rays]
  ampl .-= minimum(ampl)
  if maximum(ampl) > 0.0
    cndx = round.(Int, 1 .+ (length(colors)-1) .* ampl ./ maximum(ampl))
  else
    cndx = ones(Int, length(ampl)) .* length(colors)
  end
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "z (m)"
  for (j, r) ∈ enumerate(rays)
    r = r.path
    @series begin
      seriestype := :line
      linecolor := colors[cndx[j]]
      [r[i][1] for i ∈ 1:length(r)], [r[i][3] for i ∈ 1:length(r)]
    end
  end
end

end # module
