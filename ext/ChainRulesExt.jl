module ChainRulesExt

using ChainRulesCore
import UnderwaterAcoustics: _arr2ir, tmap, _blackman, _bisect, _mode_cache_get, _mode_cache_put!

# the bisection root-find is value-only; derivatives are recovered by the
# Newton polish step in _solve_mode_γ (implicit function theorem)
ChainRulesCore.@non_differentiable _bisect(f, lo, hi)
ChainRulesCore.@non_differentiable _mode_cache_get(pm, ω)
ChainRulesCore.@non_differentiable _mode_cache_put!(pm, ω, γ)

function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(tmap), f, X::AbstractArray)
  cache = tmap(X) do x
    rrule_via_ad(config, f, x)
  end
  Y = map(first, cache)
  function tmap_pullback(dY_raw)
    dY = unthunk(dY_raw)
    backevals = map(cache, dY) do (y, back), dy
      back(dy)
    end
    df = ProjectTo(f)(sum(first, backevals))
    dX = map(last, backevals)
    return (NoTangent(), df, dX)
  end
  return Y, tmap_pullback
end

function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(tmap), f, T::Type, X::AbstractArray)
  cache = tmap(X) do x
    rrule_via_ad(config, f, x)
  end
  Y = map(T∘first, cache)
  function tmap_pullback(dY_raw)
    dY = unthunk(dY_raw)
    backevals = map(cache, dY) do (y, back), dy
      back(dy)
    end
    df = ProjectTo(f)(sum(first, backevals))
    dX = map(last, backevals)
    return (NoTangent(), df, NoTangent(), dX)
  end
  return Y, tmap_pullback
end

end # module
