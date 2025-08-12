module PlotsExt

using RecipesBase
using Colors
using UnderwaterAcoustics
import UnderwaterAcoustics: AbstractAcousticSource, AbstractAcousticReceiver, RayArrival, value
import UnderwaterAcoustics: SampledFieldZ, SampledFieldX, SampledFieldXZ, SampledFieldXY, ModeArrival
import UnderwaterAcoustics: BasebandReplayChannel, DepthDependent, PositionDependent
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
    seriestype := :path
    linecolor := :royalblue
    xrange, [value(z1, (x, 0.0)) for x ∈ xrange]
  end
  @series begin
    seriestype := :path
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
  X = copy(x')
  X[isinf.(X)] .= prevfloat(typemax(eltype(X)))
  @series begin
    seriestype := :heatmap
    rx.xrange, rx.zrange, X
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
      seriestype := :path
      linecolor := colors[cndx[j]]
      [r[i][1] for i ∈ 1:length(r)], [r[i][3] for i ∈ 1:length(r)]
    end
  end
end

@recipe function plot(fld::SampledFieldZ; npts=1000)
  zr = range(minimum(fld.zrange), maximum(fld.zrange); length=npts)
  ticks --> :native
  legend --> false
  yguide --> "z (m)"
  @series begin
    seriestype := :path
    [fld(z) for z ∈ zr], zr
  end
end

@recipe function plot(fld::DepthDependent, zrange; npts=1000)
  zr = range(minimum(zrange), maximum(zrange); length=npts)
  ticks --> :native
  legend --> false
  yguide --> "z (m)"
  @series begin
    seriestype := :path
    [value(fld, z) for z ∈ zr], zr
  end
end

@recipe function plot(fld::DepthDependent, zmin, zmax; npts=1000)
  zr = range(zmin, zmax; length=npts)
  ticks --> :native
  legend --> false
  yguide --> "z (m)"
  @series begin
    seriestype := :path
    [value(fld, z) for z ∈ zr], zr
  end
end

@recipe function plot(fld::SampledFieldX; npts=1000)
  xr = range(minimum(fld.xrange), maximum(fld.xrange); length=npts)
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  @series begin
    seriestype := :path
    xr, [fld(x) for x ∈ xr]
  end
end

@recipe function plot(fld::PositionDependent, xrange; npts=1000)
  xr = range(minimum(xrange), maximum(xrange); length=npts)
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  @series begin
    seriestype := :path
    [value(fld, x) for x ∈ xr], xr
  end
end

@recipe function plot(fld::PositionDependent, xmin, xmax; npts=1000)
  xr = range(xmin, xmax; length=npts)
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  @series begin
    seriestype := :path
    [value(fld, x) for x ∈ xr], xr
  end
end

@recipe function plot(fld::SampledFieldXY; npts=1000)
  xr = range(minimum(fld.xrange), maximum(fld.xrange); length=npts)
  yr = range(minimum(fld.yrange), maximum(fld.yrange); length=npts)
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "y (m)"
  colorbar --> true
  @series begin
    seriestype := :heatmap
    xr, yr, [fld(x, y) for x ∈ xr, y ∈ yr]
  end
end

@recipe function plot(fld::SampledFieldXZ; npts=1000)
  xr = range(minimum(fld.xrange), maximum(fld.xrange); length=npts)
  zr = range(minimum(fld.zrange), maximum(fld.zrange); length=npts)
  ticks --> :native
  legend --> false
  xguide --> "x (m)"
  yguide --> "z (m)"
  colorbar --> true
  @series begin
    seriestype := :heatmap
    xr, zr, [fld(x, z) for x ∈ xr, z ∈ zr]
  end
end

@recipe function plot(m::ModeArrival, D=nothing; npts=1000)
  D = something(D, -minimum(m.ψ.zrange))
  zr = range(0, -D; length=npts)
  ticks --> :native
  legend --> false
  xguide --> "Mode #"
  yguide --> "z (m)"
  xlims --> (0, 2)
  @series begin
    seriestype := :line
    color := :lightgray
    linestyle := :dot
    [1, 1], [0, -D]
  end
  x = m.ψ.(zr)
  s = max(maximum(abs, real(x)), maximum(abs, imag(x)))
  s > 0 && (s = 0.25/s)
  @series begin
    seriestype := :path
    1 .+ s * real(x), zr
  end
  @series begin
    seriestype := :path
    linestyle := :dash
    1 .+ s * imag(x), zr
  end
end

@recipe function plot(m::AbstractVector{<:ModeArrival}, D=nothing; npts=1000)
  n = length(m)
  D === nothing && (D = mapreduce(m1 -> -minimum(m1.ψ.zrange), min, m))
  zr = range(0, -D; length=npts)
  ticks --> :native
  legend --> false
  xguide --> "Mode #"
  yguide --> "z (m)"
  xlims --> (0, n+1)
  s = mapreduce(min, m) do m1
    x = m1.ψ.(zr)
    s = max(maximum(abs, real(x)), maximum(abs, imag(x)))
    s > 0 && (s = 0.25/s)
    s
  end
  for i ∈ 1:n
    @series begin
      seriestype := :line
      color := :lightgray
      linestyle := :dot
      [i, i], [0, -D]
    end
    x = m[i].ψ.(zr)
    @series begin
      seriestype := :path
      i .+ s * real(x), zr
    end
    @series begin
      seriestype := :path
      linestyle := :dash
      i .+ s * imag(x), zr
    end
  end
end

@recipe function plot(ch::BasebandReplayChannel, i=1; npts=1000)
  T = size(ch.h, 3)
  Tstep = ceil(Int, T / npts)
  h = 20 * log10.(reverse(abs.(ch.h[:,i,1:Tstep:T]'); dims=2))
  hmax = ceil(Int, maximum(h))
  n, m = size(h)
  t = (0:n-1) / ch.fs * ch.step * Tstep
  τ = (0:m-1) / ch.fs
  ticks --> :native
  legend --> false
  xguide --> "delay (ms)"
  yguide --> "time (s)"
  yflip --> true
  clims --> (hmax-30, hmax)
  colorbar --> true
  @series begin
    seriestype := :heatmap
    1000τ, t, h
  end
end

end # module
