export PekerisRayModel

struct PekerisRayModel{T} <: PropagationModel{T}
  env::T
  rays::Int
  PekerisRayModel(env, rays) = new{typeof(env)}(checkenvironment(PekerisRayModel, env), rays)
end

function checkenvironment(::Type{PekerisRayModel}, env::UnderwaterEnvironment)
  altimetry(env) isa FlatSurface || throw(ArgumentError("PekerisRayModel only supports environments with flat sea surface"))
  bathymetry(env) isa ConstantDepth || throw(ArgumentError("PekerisRayModel only supports constant depth environments"))
  ssp(env) isa IsoSSP || throw(ArgumentError("PekerisRayModel only supports isovelocity environments"))
  env
end

environment(model::PekerisRayModel) = model.env

function arrivals(model::PekerisRayModel, tx1::AcousticSource, rx1::AcousticReceiver)
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
  D = √(R² + abs2(d1 - d2))
  t = D/c
  A = cis(D/k) / D * absorption(f, D, salinity(model.env))
  arr = Array{Arrival{typeof(t),typeof(A),Missing}}(undef, model.rays)
  arr[1] = Arrival(t, A, 0, 0)
  for j ∈ 2:model.rays
    upward = iseven(j)
    s1 = 2*upward - 1
    n = div(j, 2)
    s = div(n + upward, 2)
    b = div(n + (1-upward), 2)
    s2 = 2*iseven(n) - 1
    dz = 2*b*h + s1*d1 - s1*s2*d2
    D = √(R² + abs2(dz))
    θ = atan(√R²/dz)
    t = D/c
    A = cis(D/k) / D * absorption(f, D, salinity(model.env))
    s > 0 && (A *= reflectioncoef(seasurface(model.env), f, θ)^s)
    b > 0 && (A *= reflectioncoef(seabed(model.env), f, θ)^b)
    arr[j] = Arrival(t, A, s, b)
  end
  arr
end

function transfercoef(model::PekerisRayModel, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
  arr = arrivals(model, tx1, rx1)
  if mode === :coherent
    tc = sum(a.phasor for a ∈ arr)
  elseif mode === :incoherent
    tc = √sum(abs2(a.phasor) for a ∈ arr)
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  tc
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
  D = √(R² + abs2(d1 - d2))
  t = D/c
  A = cis(D/k) / D * absorption(f, D, salinity(model.env))
  arr = Array{Arrival{typeof(t),typeof(A),eltype(p1)}}(undef, model.rays)
  arr[1] = Arrival(t, A, 0, 0, [p1,p2])
  for j ∈ 2:model.rays
    upward = iseven(j)
    s1 = 2*upward - 1
    n = div(j, 2)
    s = div(n + upward, 2)
    b = div(n + (1-upward), 2)
    s2 = 2*iseven(n) - 1
    dz = 2*b*h + s1*d1 - s1*s2*d2
    D = √(R² + abs2(dz))
    R = √R²
    θ = atan(R/dz)
    t = D/c
    A = cis(D/k) / D * absorption(f, D, salinity(model.env))
    s > 0 && (A *= reflectioncoef(seasurface(model.env), f, θ)^s)
    b > 0 && (A *= reflectioncoef(seabed(model.env), f, θ)^b)
    raypath = Array{typeof(p1)}(undef, 2+s+b)
    raypath[1] = p1
    dx = p2[1] - p1[1]
    dy = p2[2] - p1[2]
    z = -(1-upward) * h
    r = abs(z-d1) * tan(θ)
    raypath[2] = (p1[1] + r/R * dx, p1[2] + r/R * dy, z)
    for i ∈ 3:length(raypath)-1
      r += h * tan(θ)
      z = -h - z
      raypath[i] = (p1[1] + r/R * dx, p1[2] + r/R * dy, z)
    end
    raypath[end] = p2
    arr[j] = Arrival(t, A, s, b, raypath)
  end
  arr
end

function record(model::PekerisRayModel, tx::AbstractArray{AcousticSource}, rx::AbstractArray{AcousticReceiver}, duration, fs; start=0.0)

end
