export models

const allmodels = [
  PekerisRayModel,
  RaySolver,
  Bellhop
]

"""
    models()
    models(env::UnderwaterEnvironment)

List available models, optionally for a given environment.
"""
function models(env=missing)
  mlist = []
  for m ∈ allmodels
    try
      check(m, env)
      push!(mlist, m)
    catch ex
      # don't add to list
    end
  end
  mlist
end

"""
$(SIGNATURES)
Register model in the list of available models.
"""
addmodel!(mtype) = push!(allmodels, mtype)
