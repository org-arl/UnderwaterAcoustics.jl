import Accessors: @set
export AdiabaticExt

"""
    AdiabaticExt(model, env; dx=0.0, dz=0.0, reciprocal=false, kwargs...)

A 2D adiabatic mode propagation model based on a range-independent modal
propagation `model`. The adiabatic mode model supports range-dependent
bathymetry. Any `kwargs` passed in are transferred to the underlying `model`.

`dx` is the step size in the range direction for integration. `dz` is the mesh
size in the depth direction to compute modes. If `dx` and/or `dz` is set to
zero, it is automatically determined (~10 points per wavelength).

Adiabatic mode models are usually not reciprocal, i.e., exchanging a source and
receiver changes the answer. This is because the bathymetry is assumed to have
azimuthal symmetry around the source. In `reciprocal` mode, the adiabatic
computation is modified to ensure reciprocity, but azimuthal symmetry is lost.
"""
struct AdiabaticExt{T1,T2,T3} <: AbstractModePropagationModel
  env::T1
  kwargs::T2
  dx::Float64
  dz::Float64
  reciprocal::Bool
  function AdiabaticExt(model::Type{<:AbstractModePropagationModel}, env; dx=0.0, dz=0.0, reciprocal=false, kwargs...)
    new{typeof(env),typeof(kwargs),model}(env, kwargs, dx, dz, reciprocal)
  end
end

function Base.show(io::IO, pm::AdiabaticExt)
  print(io, "AdiabaticExt($(_model(pm)))")
end

## interface methods

function arrivals(pm::AdiabaticExt, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver)
  min_depth = minimum(pm.env.bathymetry)
  max_depth = maximum(pm.env.bathymetry)
  if min_depth == max_depth
    env = let env = pm.env
      @set env.bathymetry = min_depth
    end
    pm1 = _model(pm)(env; pm.kwargs...)
    return arrivals(pm1, tx, rx)
  end
  # only 2D models in the x-z plane supported at present
  location(tx).y == 0 || error("Source must be in the x-z plane")
  location(rx).y == 0 || error("All receivers must be in the x-z plane")
  # dz is allowable error in z for cached mode, dx is step size in x for integration
  dz = pm.dz > 0 ? pm.dz : minimum(pm.env.soundspeed) / (frequency(tx) * 10)
  dx = pm.dx > 0 ? pm.dx : minimum(pm.env.soundspeed) / (frequency(tx) * 10)
  # precompute modes and cache them
  tx_depth = value(pm.env.bathymetry, location(tx))
  cache = _precompute_arrivals(pm, tx, tx_depth, rx, min_depth, max_depth, dz)
  # compute all modes
  xs = range(location(tx).x, location(rx).x; step = location(tx).x ≤ location(rx).x ? dx : -dx)
  modes = map(x -> _cached_arrivals(cache, -value(pm.env.bathymetry, (x=x, y=0.0, z=0.0)), dz), xs)
  rx_modes = _cached_arrivals(cache, -value(pm.env.bathymetry, (x=location(rx).x, y=0.0, z=0.0)), dz)
  nmodes = min(mapreduce(length, min, modes), length(rx_modes))
  ω = 2π * frequency(tx)
  map(1:nmodes) do m
    v = length(modes) / sum(mode -> 1 / mode[m].v, modes)
    ModeArrival(m, rx_modes[m].kᵣ, rx_modes[m].ψ, v, ω / real(rx_modes[m].kᵣ))
  end
end

function acoustic_field(pm::AdiabaticExt, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; mode=:coherent)
  only(acoustic_field(pm, tx, [rx]; mode))
end

function acoustic_field(pm::AdiabaticExt, tx::AbstractAcousticSource, rxs::AbstractArray{<:AbstractAcousticReceiver}; mode=:coherent)
  min_depth = minimum(pm.env.bathymetry)
  max_depth = maximum(pm.env.bathymetry)
  if min_depth == max_depth
    env = let env = pm.env
      @set env.bathymetry = min_depth
    end
    pm1 = _model(pm)(env; pm.kwargs...)
    return acoustic_field(pm1, tx, rxs; mode)
  end
  # only 2D models in the x-z plane supported at present
  location(tx).y == 0 || error("Source must be in the x-z plane")
  all(map(rx -> location(rx).y == 0, rxs)) || error("All receivers must be in the x-z plane")
  # dz is allowable error in z for cached mode, dx is step size in x for integration
  dz = pm.dz > 0 ? pm.dz : minimum(pm.env.soundspeed) / (frequency(tx) * 10)
  dx = pm.dx > 0 ? pm.dx : minimum(pm.env.soundspeed) / (frequency(tx) * 10)
  # precompute modes and cache them using a virtual rx as far as possible just
  # to help with choice of solver parameters for some models
  tx_depth = value(pm.env.bathymetry, location(tx))
  max_r = maximum(rx -> abs(location(rx).x - location(tx).x), rxs)
  virtual_rx = AcousticReceiver(location(tx).x + max_r, 0, tx_depth)
  cache = _precompute_arrivals(pm, tx, tx_depth, virtual_rx, min_depth, max_depth, dz)
  # get x positions of interest
  xs = unique!(sort!(map(rx -> location(rx).x, vec(rxs))))
  ndx = map(x -> findall(rx -> location(rx).x == x, rxs), xs)
  x0 = location(tx).x
  z0 = -value(pm.env.bathymetry, (x=x0, y=0.0, z=0.0))
  i0 = something(findfirst(≥(x0), xs), length(xs)+1)
  fld = zeros(ComplexF64, size(rxs))
  tx_modes = _cached_arrivals(cache, z0, dz)
  for m ∈ eachindex(tx_modes)
    kri = zeros(ComplexF64, size(xs))
    # compute kr-integral for all receivers on the right of the source
    x = x0
    for i ∈ i0:length(xs)
      i > i0 && (kri[i] = kri[i-1])
      m1 = _cached_arrivals(cache, -value(pm.env.bathymetry, (x=x, y=0.0, z=0.0)), dz)
      while x < xs[i] - dx
        x += dx
        m2 = _cached_arrivals(cache, -value(pm.env.bathymetry, (x=x, y=0.0, z=0.0)), dz)
        min(length(m1), length(m2)) < m && break
        kri[i] += 0.5 * (m1[m].kᵣ + m2[m].kᵣ) * dx
        m1 = m2
      end
      if x < xs[i]
        m2 = _cached_arrivals(cache, -value(pm.env.bathymetry, (x=xs[i], y=0.0, z=0.0)), dz)
        if x < xs[i] - dx || min(length(m1), length(m2)) < m
          kri[i] = 0
          break
        end
        kri[i] += 0.5 * (m1[m].kᵣ + m2[m].kᵣ) * (xs[i] - x)
        x = xs[i]
      end
    end
    # compute kr-integral for all receivers on the left of the source
    x = x0
    for i ∈ i0-1:-1:1
      i < i0-1 && (kri[i] = kri[i+1])
      m1 = _cached_arrivals(cache, -value(pm.env.bathymetry, (x=x, y=0.0, z=0.0)), dz)
      while x > xs[i] + dx
        x -= dx
        m2 = _cached_arrivals(cache, -value(pm.env.bathymetry, (x=x, y=0.0, z=0.0)), dz)
        min(length(m1), length(m2)) < m && break
        kri[i] += 0.5 * (m1[m].kᵣ + m2[m].kᵣ) * dx
        m1 = m2
      end
      if x > xs[i]
        m2 = _cached_arrivals(cache, -value(pm.env.bathymetry, (x=xs[i], y=0.0, z=0.0)), dz)
        if x > xs[i] + dx || min(length(m1), length(m2)) < m
          kri[i] = 0
          break
        end
        kri[i] += 0.5 * (m1[m].kᵣ + m2[m].kᵣ) * (x - xs[i])
        x = xs[i]
      end
    end
    # compute field contribution
    for (x1, kri1, ndx1) ∈ zip(xs, kri, ndx)
      if x1 == x0
        fld[ndx1] .= NaN
        continue
      end
      kri1 == 0 && break
      z = -value(pm.env.bathymetry, (x=x1, y=0.0, z=0.0))
      rx_modes = _cached_arrivals(cache, z, dz)
      if m ≤ length(rx_modes)
        R = abs(x1 - x0)
        A = tx_modes[m].ψ(location(tx).z) * cis(-kri1) /
            sqrt(pm.reciprocal ? kri1 : rx_modes[m].kᵣ * R)
        for j ∈ ndx1
          zᵣ = location(rxs[j]).z
          if zᵣ ≥ z
            fld[j] += mode === :coherent ? A * rx_modes[m].ψ(zᵣ) : abs2(A * rx_modes[m].ψ(zᵣ))
          end
        end
      end
    end
  end
  mode === :coherent || (fld .= sqrt.(fld))
  fld .*= sqrt(2π) * cispi(0.25) * db2amp(spl(tx))
  fld
end

## private methods

_model(pm::AdiabaticExt{T1,T2,T3}) where {T1,T2,T3} = T3

struct ModeCache{T}
  zs::Vector{Float64}
  modes::Vector{Vector{T}}
  ModeCache(T) = new{T}(Float64[], Vector{T}[])
end

function _cached_arrivals(cache, z, dz)
  i = searchsortedfirst(cache.zs, z)
  i ≤ length(cache.zs) && cache.zs[i] == z && return cache.modes[i]
  dz1 = i > 1 ? z - cache.zs[i-1] : Inf
  dz2 = i ≤ length(cache.zs) ? cache.zs[i] - z : Inf
  if dz1 < dz2
    dz1 > 1.01dz && error("Cache miss for required mode")
    cache.modes[i-1]
  else
    dz2 > 1.01dz && error("Cache miss for required mode")
    cache.modes[i]
  end
end

function _cache_arrivals(cache, z, v)
  i = searchsortedfirst(cache.zs, z)
  if i > length(cache.zs)
    push!(cache.zs, z)
    push!(cache.modes, v)
  elseif cache.zs[i] == z
    cache.modes[i] = v
  else
    insert!(cache.zs, i, z)
    insert!(cache.modes, i, v)
  end
  nothing
end

function _precompute_arrivals(pm, tx, tx_depth, rx, min_depth, max_depth, dz)
  n = ceil(Int, (max_depth - min_depth) / 2dz) + 1
  arr = _compute_arrivals(pm, tx, -tx_depth, rx)
  cache = ModeCache(eltype(arr))
  _cache_arrivals(cache, -tx_depth, arr)
  for z ∈ range(-max_depth, -min_depth; length=n)
    i = searchsortedlast(cache.zs, z)
    if i == 0 || cache.zs[i] != z
      arr = _compute_arrivals(pm, tx, z, rx)
      _align_arrivals!(arr, cache.modes[max(1,i)])
      _cache_arrivals(cache, z, arr)
    end
  end
  cache
end

function _compute_arrivals(pm, tx, z, rx)
  # compute arrivals for a water depth of z
  env = let env = pm.env
    @set env.bathymetry = -z
  end
  pm1 = _model(pm)(env; pm.kwargs...)
  arr = arrivals(pm1, tx, rx)
  convert(Vector{Union{Nothing,eltype(arr)}}, arr)
end

function _align_arrivals!(arr1, arr0; candidates=3, corr_threshold=0.5, npts=101)
  arr2 = similar(arr0)
  arr2 .= nothing
  for m1 ∈ arr1
    zs = range(extrema(m1.ψ.zrange)...; length=npts)
    m1ψ = m1.ψ.(zs)
    ndx = sortperm(map(zip(arr0, arr2)) do (m0, m2)
      m0 === nothing || m2 !== nothing ? Inf : abs(m0.kᵣ - m1.kᵣ)
    end)
    resize!(ndx, min(candidates, length(ndx)))
    v, i = findmax(map(m0 -> abs(m1ψ' * m0.ψ.(zs)), arr0[ndx]))
    v / abs(m1ψ' * m1ψ) > corr_threshold && (arr2[ndx[i]] = m1)
  end
  arr2
end
