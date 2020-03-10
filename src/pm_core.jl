export PropagationModel, source!, receiver!, record

abstract type PropagationModel end

Base.@kwdef struct AcousticSource
  pos = []
  signal = nothing
end

Base.@kwdef struct AcousticReceiver
  pos = []
end

source!(model::PropagationModel, pos, signal=nothing) = append!(model.sources, AcousticSource(pos, signal))
receiver!(model::PropagationModel, pos) = append!(model.receivers, AcousticReceiver(pos))

function record(model::PropagationModel, duration)
  nsamples = round(Int, duration*model.fs)
  x = zeros(nsamples, length(model.receivers))
  for j = 1:length(model.sources)
    sig = generate(model.sources[j].signal, nsamples, model.fs)
    for k = 1:length(model.receivers)
      x[:,k] .+= transmit(model, sig, model.sources[j].pos, model.receivers[k].pos)
    end
  end
  for k = 1:length(model.receivers)
    x[:,k] .+= rand(model.noise, nsamples) * √model.fs
  end
  return x
end
