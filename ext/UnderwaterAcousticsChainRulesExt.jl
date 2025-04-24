module UnderwaterAcousticsChainRulesExt

using ChainRulesCore
import UnderwaterAcoustics: _arr2ir

function ChainRulesCore.rrule(::typeof(_arr2ir), ts, ϕs; T, t0, fs, n)
  ir = _arr2ir(ts, ϕs; T, t0, fs, n)
  function _arr2ir_pullback(dx)
    dx = unthunk(dx)
    T1 = promote_type(eltype(ts), eltype(dx))
    T2 = promote_type(eltype(ϕs), eltype(dx))
    dt = zeros(T1, length(ts))
    dϕ = zeros(T2, length(ϕs))
    for i ∈ eachindex(dt)
      t = (ts[i] - t0) * fs + 1
      t̄ = floor(Int, t)
      α, β = sincospi(0.5 * (t - t̄))
      if t̄ ≤ n
        dϕ[i] += β * dx[t̄]
        dt[i] -= 0.5π * fs * α * ϕs[i] * dx[t̄]
      end
      if t̄ < n
        dϕ[i] += α * dx[t̄+1]
        dt[i] += 0.5π * fs  * β * ϕs[i] * dx[t̄+1]
      end
    end
    return NoTangent(), dt, dϕ
  end
  return ir, _arr2ir_pullback
end

end # module
