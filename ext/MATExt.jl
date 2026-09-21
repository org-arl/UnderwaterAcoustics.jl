module MATExt

using UnderwaterAcoustics
import MAT: matread
import SignalAnalysis: resample
import UnderwaterAcoustics: BasebandReplayChannel

function UnderwaterAcoustics._load_mat_replay_channel(filename, upsample, rxs, noise)
  # TODO: support UACR noise models
  data = matread(filename)
  all(["version", "h_hat", "params"] .∈ Ref(keys(data))) || error("Bad channel file format")
  data["version"] >= 1.0 || @warn "Unsupported channel file version"
  # the file stores taps in forward-delay order; _apply_tvir! expects them reversed
  h = reverse(data["h_hat"]; dims=1)
  M = size(h, 2)
  rxs === (:) && (rxs = 1:M)
  ndims(rxs) == 0 && (rxs = [rxs])
  h = h[:,rxs,:]
  θ = Matrix{Float64}(undef, 0, 0)
  φ = Matrix{Float64}(undef, 0, 0)
  if haskey(data, "phi_hat")
    φ_data = data["phi_hat"]
    size(φ_data, 1) == M || error("Invalid phi_hat size")
    φ = transpose(φ_data[rxs,:])
  elseif haskey(data, "theta_hat")
    θ_data = data["theta_hat"]
    size(θ_data, 1) == M || error("Invalid theta_hat size")
    θ = transpose(θ_data[rxs,:])
  end
  fs = data["params"]["fs_delay"]
  fs_time = data["params"]["fs_time"]
  fc = data["params"]["fc"]
  doppler = haskey(data, "f_resamp") ? only(data["f_resamp"]) : 1.0
  ratio = fs / fs_time
  step = round(Int, ratio)
  isapprox(ratio, step; rtol=1e-9) || error("fs_delay/fs_time must be an integer ratio (got $ratio)")
  if upsample && step != 1
    h = UnderwaterAcoustics._interp_ir(h, step, (size(h, 3) - 1) * step + 1)
    step = 1
  end
  let nphase = size(φ, 1) > 0 ? size(φ, 1) : size(θ, 1)
    if nphase > 0
      dur_phase = nphase / fs
      dur_h = size(data["h_hat"], 3) / fs_time
      isapprox(dur_phase, dur_h; rtol=1e-3) ||
        error("Phase/IR duration mismatch: phase spans $(round(dur_phase;digits=3))s, " *
              "h_hat spans $(round(dur_h;digits=3))s (spec requires equal durations)")
    end
  end
  BasebandReplayChannel(h, θ, φ, fs, fc, step, doppler; noise)
end

end
