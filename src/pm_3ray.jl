export transmit, ThreeRayModel, BasicSeabed

Base.@kwdef mutable struct ThreeRayModel <: PropagationModel
  sources::Vector{AcousticSource} = []
  receivers::Vector{AcousticReceiver} = []
  waterdepth = 25.0
  soundspeed = 1500.0
  density = 1023.0
  seabed = BasicSeabed(1600.0, 1600.0, 0.1)
  surfaceloss = db2amp(1.0)
  bottomloss = db2amp(1.0)
  noisemodel = RedGaussian
  noiselevel = 3000
  fs = 96000.0
end

Base.@kwdef struct BasicSeabed
  soundspeed = 1600.0
  density = 1600.0
  absorption = 0.1
end

function transmit(model::ThreeRayModel, x::AbstractVector, spos, rpos)
  rpos2 = copy(rpos)
  rpos2[end] = -rpos2[end]
  rpos3 = copy(rpos)
  rpos3[end] = -2*model.waterdepth - rpos3[end]
  d1 = norm(rpos - spos)
  d2 = norm(rpos2 - spos)
  d3 = norm(rpos3 - spos)
  t1 = round(Int, d1/model.soundspeed*model.fs)
  t2 = round(Int, d2/model.soundspeed*model.fs)
  t3 = round(Int, d3/model.soundspeed*model.fs)
  breflect = abs(reflectioncoef(
    acos(-(rpos[end]+spos[end])/d2),
    model.seabed.density, model.seabed.soundspeed, model.seabed.absorption,
    model.density, model.soundspeed))
  y = zeros(length(x))
  y[t1:end] .+=  x[1:end-t1+1]/d1
  y[t2:end] .+= -x[1:end-t2+1]/d2 / model.surfaceloss
  y[t3:end] .+=  x[1:end-t3+1]/d3 / model.bottomloss * breflect
  return y
end
