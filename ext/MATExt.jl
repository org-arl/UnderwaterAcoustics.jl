module MATExt

using UnderwaterAcoustics
import MAT: matread
import SignalAnalysis: resample
import UnderwaterAcoustics: BasebandReplayChannel

function UnderwaterAcoustics._load_mat_replay_channel(filename, upsample, rxs, noise)
  # TODO: support UACR noise models
  data = matread(filename)
  all(["version", "h_hat", "theta_hat", "params"] .∈ Ref(keys(data))) || error("Bad channel file format")
  data["version"] == 1.0 || @warn "Unsupported channel file version"
  h = reverse(data["h_hat"]; dims=1)
  M = size(h, 2)
  rxs === (:) && (rxs = 1:M)
  ndims(rxs) == 0 && (rxs = [rxs])
  h = h[:,rxs,:]
  θ = data["theta_hat"]
  size(θ, 1) == M || error("Invalid phase estimates size")
  θ = transpose(θ[rxs,:])
  fs = data["params"]["fs_delay"]
  fc = data["params"]["fc"]
  if upsample && fs != data["params"]["fs_time"]
    step = 1
    h = resample(h, fs / data["params"]["fs_time"]; dims=3)
  else
    step = round(Int, fs / data["params"]["fs_time"])
  end
  BasebandReplayChannel(h, θ, fs, fc, step; noise)
end

end # module
