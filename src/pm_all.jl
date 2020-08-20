export models

const allmodels = [
  PekerisRayModel,
  RaySolver,
  Bellhop
]

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

addmodel!(mtype) = push!(allmodels, mtype)
