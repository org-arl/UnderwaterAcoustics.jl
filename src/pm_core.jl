export PropagationModel, source!, receiver!, record

abstract type PropagationModel end

Base.@kwdef struct AcousticSource
  pos = []
  signal = nothing
end

Base.@kwdef struct AcousticReceiver
  pos = []
end

source!(model::PropagationModel, pos, signal=nothing) = push!(model.sources, AcousticSource(pos, signal))
receiver!(model::PropagationModel, pos) = push!(model.receivers, AcousticReceiver(pos))

source!(model::PropagationModel, x::AcousticSource) = push!(model.sources, x)
receiver!(model::PropagationModel, x::AcousticReceiver) = push!(model.receivers, x)

function record(model::PropagationModel, duration; start=0.0)
  if start isa Symbol
    distances = [norm(model.sources[j].pos-model.receivers[k].pos)
      for j in 1:length(model.sources), k in 1:length(model.receivers)]
    if start == :first
      start = minimum(distances)/minimum(model.soundspeed)
    elseif start == :last
      start = maximum(distances)/minimum(model.soundspeed)
    else
      throw(ArgumentError("Unknown start option"))
    end
  end
  nsamples = round(Int, (start+duration)*model.fs)
  x = zeros(nsamples, length(model.receivers))
  for j = 1:length(model.sources)
    sig = generate(model.sources[j].signal, nsamples, model.fs)
    for k = 1:length(model.receivers)
      x[:,k] .+= transmit(model, sig, model.sources[j].pos, model.receivers[k].pos)
    end
  end
  for k = 1:length(model.receivers)
    x[:,k] .+= rand(model.noisemodel(nsamples)) * model.noiselevel * √model.fs
  end
  x[end-round(Int,duration*model.fs)+1:end,:]
end
