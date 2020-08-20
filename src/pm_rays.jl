using DifferentialEquations
using ForwardDiff

export RaySolver

struct RaySolver{T} <: PropagationModel{T}
  env::T
end

### interface functions

function check(::Type{RaySolver}, env::Union{<:UnderwaterEnvironment,Missing})
  if env !== missing
    altimetry(env) isa FlatSurface || throw(ArgumentError("RaySolver only supports environments with flat sea surface"))
    bathymetry(env) isa ConstantDepth || throw(ArgumentError("RaySolver only supports constant depth environments"))
  end
  env
end

function arrivals(model::RaySolver, tx1::AcousticSource, rx1::AcousticReceiver)
  # TODO
end

function eigenrays(model::RaySolver, tx1::AcousticSource, rx1::AcousticReceiver)
  # TODO
end

function rays(model::RaySolver, tx1::AcousticSource, θ::Real, rmax)
  θ₀ = θ
  f = nominalfrequency(tx1)
  zmax = -maxdepth(bathymetry(model.env))
  p = location(tx1)
  raypath = Array{typeof(p)}(undef, 0)
  a = Complex(1.0, 0.0)
  s = 0
  b = 0
  t = 0.0
  D = 0.0
  while true
    dD, r, z, θ, dt, u = traceray1(model, p[1], p[3], θ, rmax)
    for u1 ∈ u
      push!(raypath, (u1[1], 0.0, u1[2]))
    end
    t += dt
    D += dD
    if isapprox(z, 0.0; atol=1e-3)        # FIXME: assumes flat altimetry
      s += 1
      a *= reflectioncoef(seasurface(model.env), f, θ)
    elseif isapprox(z, zmax; atol=1e-3)   # FIXME: assumes flat bathymetry
      b += 1
      a *= reflectioncoef(seabed(model.env), f, θ)
    else
      break
    end
    p = (r, 0.0, z)
    θ = -θ                                # FIXME: assumes flat altimetry/bathymetry
  end
  RayArrival(t, a/D, s, b, θ₀, θ, raypath)
end

function transfercoef(model::RaySolver, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
  # TODO
end

### helper functions

function rayeqns!(du, u, params, s)
  r, z, ξ, ζ, t = u
  c, ∂c = params
  cᵥ = c(z)
  cᵥ² = cᵥ * cᵥ
  du[1] = cᵥ * ξ
  du[2] = cᵥ * ζ
  du[3] = 0               # TODO: support range-dependent soundspeed
  du[4] = -∂c(z) / cᵥ²
  du[5] = 1 / cᵥ
end

function checkray!(out, u, s, integrator, a::FlatSurface, b::ConstantDepth)
  bz = u[4] ≥ 0 ? 0.0 : -maxdepth(b)
  out[1] = u[2] - bz
  out[2] = u[3]
end

function traceray1(model, r0, z0, θ, rmax; ds=1.0)
  a = altimetry(model.env)
  b = bathymetry(model.env)
  c = z -> soundspeed(ssp(model.env), 0.0, 0.0, z)
  ∂c = z -> ForwardDiff.derivative(c, z)
  cᵥ = c(z0)
  u0 = [r0, z0, cos(θ)/cᵥ, sin(θ)/cᵥ, 0.0]
  prob = ODEProblem(rayeqns!, u0, (0.0, (rmax-r0)/cos(θ)), (c, ∂c))
  cb = VectorContinuousCallback(
    (out, u, s, i) -> checkray!(out, u, s, i, a, b),
    (i, ndx) -> terminate!(i), 2; rootfind=true)
  soln = solve(prob; saveat=ds, callback=cb)
  s2 = soln[end]
  soln.t[end], s2[1], s2[2], atan(s2[4], s2[3]), s2[5], soln.u
end
