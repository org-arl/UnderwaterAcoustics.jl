export AdiabaticExt

struct AdiabaticExt{T1,T2,T3} <: AbstractModePropagationModel
  env::T1
  kwargs::T2
  nmesh_per_λ::Int
  max_modes::Int
  reciprocal::Bool
  zs::Vector{Float64}
  cache::Vector{Vector{ModeArrival}}
  function AdiabaticExt(model, env; nmesh_per_λ=10, max_modes=0, reciprocal=false, kwargs...)
    new{typeof(env),typeof(kwargs),model}(env, kwargs, nmesh_per_λ, max_modes, reciprocal, Float64[], Vector{ModeArrival}[])
  end
end

function Base.show(io::IO, pm::AdiabaticExt)
  print(io, "AdiabaticExt($(_model(pm)))")
end

# TODO
# function arrivals(pm::AdiabaticExt, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver)
# end

function acoustic_field(pm::AdiabaticExt, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; mode=:coherent)
  only(acoustic_field(pm, tx, [rx]; mode))
end

function acoustic_field(pm::AdiabaticExt, tx::AbstractAcousticSource, rxs::AbstractArray{<:AbstractAcousticReceiver}; mode=:coherent)
  min_depth = minimum(pm.env.bathymetry)
  max_depth = maximum(pm.env.bathymetry)
  if min_depth == max_depth
    env = UnderwaterEnvironment(
      bathymetry = min_depth,
      altimetry = pm.env.altimetry,
      temperature = pm.env.temperature,
      salinity = pm.env.salinity,
      soundspeed = pm.env.soundspeed,
      density = pm.env.density,
      seabed = pm.env.seabed,
      surface = pm.env.surface
    )
    pm1 = _model(pm)(env; pm.kwargs...)
    return acoustic_field(pm1, tx, rxs; mode)
  end
  # only 2D models in the x-z plane supported at present
  location(tx).y == 0 || error("Source must be in the x-z plane")
  all(map(rx -> location(rx).y == 0, rxs)) || error("All receivers must be in the x-z plane")
  # dz is allowable error in z for cached mode, dx is step size in x for integration
  dz = minimum(pm.env.soundspeed) / (frequency(tx) * pm.nmesh_per_λ)
  dx = dz
  # precompute modes and cache them
  n = ceil(Int, (max_depth - min_depth) / 2dz) + 1
  for z ∈ range(-max_depth, -min_depth; length=n)
    v = _cached(pm, z, dz)
    if v === nothing
      env = UnderwaterEnvironment(
        bathymetry = -z,
        altimetry = pm.env.altimetry,
        temperature = pm.env.temperature,
        salinity = pm.env.salinity,
        soundspeed = pm.env.soundspeed,
        density = pm.env.density,
        seabed = pm.env.seabed,
        surface = pm.env.surface
      )
      pm1 = _model(pm)(env; pm.kwargs...)
      rx1 = AcousticReceiver(location(tx))
      arr = arrivals(pm1, tx, rx1)
      _cache(pm, z, arr)
    end
  end
  # get x positions of interest
  xs = unique!(sort!(map(rx -> location(rx).x, vec(rxs))))
  ndx = map(x -> findall(rx -> location(rx).x == x, rxs), xs)
  x0 = location(tx).x
  z0 = -value(pm.env.bathymetry, (x=x0, y=0.0, z=0.0))
  i0 = something(findfirst(≥(x0), xs), length(xs)+1)
  fld = zeros(ComplexF64, size(rxs))
  tx_modes = _cached(pm, z0, dz)
  nmodes = length(tx_modes)
  pm.max_modes > 0 && nmodes < pm.max_modes && (nmodes = pm.max_modes)
  a = absorption(frequency(tx), 1.0, pm.env.salinity, pm.env.temperature, min_depth / 2)  # nominal absorption
  for m ∈ 1:nmodes
    kri = zeros(ComplexF64, size(xs))
    # compute kr-integral for all receivers on the right of the source
    x = x0
    for i ∈ i0:length(xs)
      i > i0 && (kri[i] = kri[i-1])
      m1 = _cached(pm, -value(pm.env.bathymetry, (x=x, y=0.0, z=0.0)), dz)
      while x < xs[i] - dx
        x += dx
        m2 = _cached(pm, -value(pm.env.bathymetry, (x=x, y=0.0, z=0.0)), dz)
        min(length(m1), length(m2)) < m && break
        kri[i] += 0.5 * (m1[m].kᵣ + m2[m].kᵣ) * dx
        m1 = m2
      end
      if x < xs[i]
        m2 = _cached(pm, -value(pm.env.bathymetry, (x=xs[i], y=0.0, z=0.0)), dz)
        if x < xs[i] - dx || min(length(m1), length(m2)) < m
          kri[i] = 0
          break
        end
        kri[i] += 0.5 * (m1[m].kᵣ + m2[m].kᵣ) * (xs[i] - x)
        x = xs[i]
      end
    end
    # compute kr-integral for all receivers on the left of the source
    for i ∈ i0-1:-1:1
      # TODO
    end
    # compute field contribution
    for (x1, kri1, ndx1) ∈ zip(xs, kri, ndx)
      kri1 == 0 && break
      z = -value(pm.env.bathymetry, (x=x1, y=0.0, z=0.0))
      rx_modes = _cached(pm, z, dz)
      if m ≤ length(rx_modes)
        A = tx_modes[m].ψ(location(tx).z) * cis(kri1)
        R = x1 - x0
        A /= sqrt(pm.reciprocal ? kri1 : rx_modes[m].kᵣ * R)
        A *= a ^ R
        for j ∈ ndx1
          zᵣ = location(rxs[j]).z
          if zᵣ ≥ z
            fld[j] += mode === :coherent ? A * rx_modes[m].ψ(zᵣ) : abs2(A * rx_modes[m].ψ(zᵣ))
          end
        end
      end
    end
  end
  #ρ = value(pm.env.density, location(tx))
  #fld .*= cispi(0.25) / (ρ * √(8π)) * db2amp(spl(tx))
  mode === :coherent || (fld .= sqrt.(fld))
  fld .*= sqrt(2π) * cispi(0.25) * db2amp(spl(tx))
  fld
end

## TODO make impulse_response() work with this

_model(pm::AdiabaticExt{T1,T2,T3}) where {T1,T2,T3} = T3

function _cached(pm, z, dz)
  i = searchsortedfirst(pm.zs, z)
  i ≤ length(pm.zs) && pm.zs[i] == z && return pm.cache[i]
  dz1 = i > 1 ? z - pm.zs[i-1] : Inf
  dz2 = i ≤ length(pm.zs) ? pm.zs[i] - z : Inf
  if dz1 < dz2
    dz1 ≤ dz ? pm.cache[i-1] : nothing
  else
    dz2 ≤ dz ? pm.cache[i] : nothing
  end
end

function _cache(pm, z, v)
  i = searchsortedfirst(pm.zs, z)
  if i > length(pm.zs)
    push!(pm.zs, z)
    push!(pm.cache, v)
  elseif pm.zs[i] == z
    pm.cache[i] = v
  else
    insert!(pm.zs, i, z)
    insert!(pm.cache, i, v)
  end
  nothing
end
