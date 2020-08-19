export PekerisRayModel

struct PekerisRayModel{T} <: PropagationModel{T}
  env::T
  rays::Int
  function PekerisRayModel(env, rays)
    rays > 0 || throw(ArgumentError("Number of rays should be more than 0"))
    new{typeof(env)}(check(PekerisRayModel, env), rays)
  end
end

### interface functions

function check(::Type{PekerisRayModel}, env::Union{<:UnderwaterEnvironment,Missing})
  if env !== missing
    altimetry(env) isa FlatSurface || throw(ArgumentError("PekerisRayModel only supports environments with flat sea surface"))
    bathymetry(env) isa ConstantDepth || throw(ArgumentError("PekerisRayModel only supports constant depth environments"))
    ssp(env) isa IsoSSP || throw(ArgumentError("PekerisRayModel only supports isovelocity environments"))
  end
  env
end

function arrivals(model::PekerisRayModel, tx1::AcousticSource, rx1::AcousticReceiver)
  # based on Chitre (2007)
  c = soundspeed(ssp(model.env), 0.0, 0.0, 0.0)
  h = depth(bathymetry(model.env), 0.0, 0.0)
  f = nominalfrequency(tx1)
  p1 = location(tx1)
  p2 = location(rx1)
  R² = abs2(p1[1] - p2[1]) + abs2(p1[2] - p2[2])
  d1 = -p1[3]
  d2 = -p2[3]
  [arrival(j, model, R², d1, d2, h, c, f) for j ∈ 1:model.rays]
end

function eigenrays(model::PekerisRayModel, tx1::AcousticSource, rx1::AcousticReceiver)
  # based on Chitre (2007)
  c = soundspeed(ssp(model.env), 0.0, 0.0, 0.0)
  h = depth(bathymetry(model.env), 0.0, 0.0)
  f = nominalfrequency(tx1)
  k = c / (2π * f)
  p1 = location(tx1)
  p2 = location(rx1)
  R² = abs2(p1[1] - p2[1]) + abs2(p1[2] - p2[2])
  d1 = -p1[3]
  d2 = -p2[3]
  [arrival(j, model, R², d1, d2, h, c, f, p1, p2) for j ∈ 1:model.rays]
end

function transfercoef(model::PekerisRayModel, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
  arr = arrivals(model, tx1, rx1)
  if mode === :coherent
    f = nominalfrequency(tx1)
    tc = sum(a.phasor * cis(2π * a.time * f) for a ∈ arr)
  elseif mode === :incoherent
    tc = Complex(√sum(abs2(a.phasor) for a ∈ arr), 0)
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  tc
end

### helper functions

# memoized version of absorption
const cached_absorption = Ref{Tuple{Float64,Float64,Float64}}()
function fast_absorption(f::Float64, D, S::Float64)
  if cached_absorption[][1:2] == (f, S)
    db2amp(cached_absorption[][3] * D)
  else
    dBperm = amp2db(absorption(f, 1.0, S))
    cached_absorption[] = (f, S, dBperm)
    db2amp(dBperm * D)
  end
end

# fallback
fast_absorption(f, D, S) = absorption(f, D, S)

# complex ForwardDiff friendly version of x^n
ipow(x, n::Int) = prod(x for _ ∈ 1:n)

function arrival(j, model, R², d1, d2, h, c, f, p1=missing, p2=missing)
  upward = iseven(j)
  s1 = 2*upward - 1
  n = div(j, 2)
  s = div(n + upward, 2)
  b = div(n + (1-upward), 2)
  s2 = 2*iseven(n) - 1
  dz = 2*b*h + s1*d1 - s1*s2*d2
  D = √(R² + abs2(dz))
  R = R² == 0 ? R² : √R²   # ForwardDiff compatible version of √R²
  θ = atan(R/dz)
  t = D/c
  A = Complex(1.0, 0.0) / D * fast_absorption(f, D, salinity(model.env))
  #A = cis(2π*t*f) / D * fast_absorption(f, D, salinity(model.env))
  s > 0 && (A *= ipow(reflectioncoef(seasurface(model.env), f, θ), s))
  b > 0 && (A *= ipow(reflectioncoef(seabed(model.env), f, θ), b))
  λ = π/2 - θ
  if typeof(p1) === Missing
    Arrival(t, conj(A), s, b, s1*λ, -s1*s2*λ)    # conj(A) needed to match with Bellhop
  else
    raypath = Array{typeof(p1)}(undef, 2+s+b)
    raypath[1] = p1
    if s + b > 0
      dx = p2[1] - p1[1]
      dy = p2[2] - p1[2]
      z = (1-upward) * h
      r = abs(z-d1) * tan(θ)
      raypath[2] = (p1[1] + r/R * dx, p1[2] + r/R * dy, -z)
      for i ∈ 3:length(raypath)-1
        r += h * tan(θ)
        z = h - z
        raypath[i] = (p1[1] + r/R * dx, p1[2] + r/R * dy, -z)
      end
    end
    raypath[end] = p2
    Arrival(t, conj(A), s, b, s1*λ, -s1*s2*λ, raypath)
  end
end
