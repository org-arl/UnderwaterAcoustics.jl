using DifferentialEquations
using ForwardDiff
using Optim

export RaySolver

Base.@kwdef struct RaySolver{T} <: PropagationModel{T}
  env::T
  nbeams::Int = 0
  minangle::Float64 = -80°
  maxangle::Float64 = +80°
  atol::Float64 = 1e-3
  rugocity::Float64 = 1.5
  athreshold::Float64 = 1e-5
  function RaySolver(env, nbeams, minangle, maxangle, atol, rugocity, athreshold)
    -π/2 ≤ minangle ≤ π/2 || throw(ArgumentError("minangle should be between -π/2 and π/2"))
    -π/2 ≤ maxangle ≤ π/2 || throw(ArgumentError("maxangle should be between -π/2 and π/2"))
    minangle < maxangle || throw(ArgumentError("maxangle should be more than minangle"))
    nbeams ≤ 0 && (nbeams = round(Int, (maxangle - minangle)/(1°)))
    new{typeof(env)}(check(Bellhop, env), nbeams, minangle, maxangle, atol, rugocity, athreshold)
  end
end

RaySolver(env; kwargs...) = RaySolver(; env=env, kwargs...)

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
  p = location(rx1)
  θ = range(model.minangle, model.maxangle; length=model.nbeams)
  r = tmap(θ1 -> rays(model, tx1, θ1, p[1]), θ)
  ie = Array{Tuple{Int,Float64}}(undef, 0)
  for i ∈ 1:length(θ)
    r1 = r[i].raypath[end]
    if isapprox(r1[1], p[1]; atol=model.atol) && isapprox(r1[2], p[2]; atol=model.atol)
      push!(ie, (i, r1[3] - p[3]))
    end
  end
  n = length(ie)
  er = Array{eltype(r)}(undef, 0)
  for i ∈ 1:n
    if isapprox(ie[i][2], 0.0; atol=model.atol)
      push!(er, r[ie[i][1]])
    elseif i > 1 && sign(ie[i-1][2]) != sign(ie[i][2])
      a, b = ordered(θ[ie[i-1][1]], θ[ie[i][1]])
      soln = optimize(x -> abs(rays(model, tx1, x, p[1]).raypath[end][3] - p[3]), a, b; abs_tol=model.atol)
      push!(er, rays(model, tx1, soln.minimizer, p[1]))
    elseif i > 2 && isnearzero(ie[i-2][2], ie[i-1][2], ie[i][2])
      a, b = ordered(θ[ie[i-2][1]], θ[ie[i][1]])
      soln = optimize(x -> abs(rays(model, tx1, x, p[1]).raypath[end][3] - p[3]), a, b; abs_tol=model.atol)
      push!(er, rays(model, tx1, soln.minimizer, p[1]))
    end
  end
  er
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
    if abs(a)/D < model.athreshold
      break
    end
    p = (r, 0.0, z)
    θ = -θ                                # FIXME: assumes flat altimetry/bathymetry
  end
  a /= D
  a *= fast_absorption(f, D, salinity(model.env))
  RayArrival(t, a, s, b, θ₀, θ, raypath)
end

function transfercoef(model::RaySolver, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
  # TODO
end

### helper functions

ordered(a, b) = a < b ? (a, b) : (b, a)

function isnearzero(a, b, c)
  sign(a) != sign(b) && return false
  sign(a) != sign(c) && return false
  abs(a) < abs(b) && return false
  abs(c) < abs(b) && return false
  return true
end

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

function checkray!(out, u, s, integrator, a::Altimetry, b::Bathymetry, rmax)
  bz = u[4] ≥ 0 ? altitude(a, u[1], 0.0) : -depth(b, u[1], 0.0)
  out[1] = u[2] - bz
  out[2] = u[3]
  out[3] = u[1] - rmax
end

# FIXME: type inference fails for prob and soln, but will be fixed in PR570 for DiffEqBase soon
function traceray1(model, r0, z0, θ, rmax; ds=1.0)
  a = altimetry(model.env)
  b = bathymetry(model.env)
  c = z -> soundspeed(ssp(model.env), 0.0, 0.0, z)
  ∂c = z -> ForwardDiff.derivative(c, z)
  cᵥ = c(z0)
  u0 = [r0, z0, cos(θ)/cᵥ, sin(θ)/cᵥ, 0.0]
  prob = ODEProblem(rayeqns!, u0, (0.0, model.rugocity * (rmax-r0)/cos(θ)), (c, ∂c))
  cb = VectorContinuousCallback(
    (out, u, s, i) -> checkray!(out, u, s, i, a, b, rmax),
    (i, ndx) -> terminate!(i), 3; rootfind=true)
  soln = solve(prob; saveat=ds, callback=cb)
  s2 = soln[end]
  soln.t[end], s2[1], s2[2], atan(s2[4], s2[3]), s2[5], soln.u
end
