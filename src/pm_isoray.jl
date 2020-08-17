export IsoRayModel

struct IsoRayModel{T} <: PropagationModel{T}
  env::T
  rays::Int
  IsoRayModel(env, rays) = new(checkenv(IsoRayModel, env), rays)
end

function checkenv(::Type{IsoRayModel}, env::UnderwaterEnvironment)
  altimetry(env) isa FlatSurface || throw(ArgumentError("IsoRayModel only supports environments with flat sea surface"))
  bathymetry(env) isa ConstantDepth || throw(ArgumentError("IsoRayModel only supports constant depth environments"))
  ssp(env) isa IsoSSP || throw(ArgumentError("IsoRayModel only supports isovelocity environments"))
  env
end

function arrivals(model::IsoRayModel, tx1::AcousticSource, rx1::AcousticReceiver)
  # based on Chitre (2007)
  c = soundspeed(ssp(model.env), 0.0, 0.0, 0.0)
  h = depth(bathymetry(model.env), 0.0, 0.0)
  f = nominalfrequency(tx1)
  p1 = location(tx1)
  p2 = location(rx1)
  R² = abs2(p1[1] - p2[1]) + abs2(p1[2] - p2[2])
  d1 = -p1[3]
  d2 = -p2[3]
  D = √(R² + abs2(d1 - d2))
  t = D/c
  A = 1/D * absorption(f, D, salinity(model.env))
  arr = Array{Tuple{typeof(t),typeof(A)}}(undef, model.rays)
  arr[1] = (t, A)
  for j ∈ 2:model.rays
    je = iseven(j)
    s1 = 2*je - 1
    n = div(j, 2)
    s = div(n + je, 2)
    b = div(n + (1-je), 2)
    s2 = 2*iseven(n) - 1
    dz = 2*b*h + s1*d1 - s1*s2*d2
    D = √(R² + abs2(dz))
    θ = atan(R/dz)
    t = D/c
    A = 1/D * absorption(f, D, salinity(model.env))
    s > 0 && (A *= reflectioncoef(seasurface(model.env), f, θ)^s)
    b > 0 && (A *= reflectioncoef(seabed(model.env), f, θ)^b)
    arr[j] = (t, A)
  end
  arr
end

function transfercoef(model::IsoRayModel, tx1::AcousticSource, rx::AbstractArray{AcousticReceiver}; mode=:coherent)

end

function eigenrays(model::IsoRayModel, tx1::AcousticSource, rx1::AcousticReceiver)

end

function record(model::IsoRayModel, tx::AbstractArray{AcousticSource}, rx::AbstractArray{AcousticReceiver}, duration, fs; start=0.0)

end
