using LinearAlgebra
using DifferentialEquations
using ForwardDiff
using Optim

export RaySolver

"""
$(TYPEDEF)
A pure Julia implementation of a ray/Gaussian beam propapagation model. The model supports complex
environments, but retains differentiability.
"""
Base.@kwdef struct RaySolver{T1,T2} <: PropagationModel{T1}
  env::T1
  nbeams::Int = 0
  minangle::Float64 = -80°
  maxangle::Float64 = +80°
  ds::Float64 = 1.0
  atol::Float64 = 1e-4
  rugocity::Float64 = 1.5
  athreshold::Float64 = 1e-5
  solver::T2 = missing
  solvertol::Float64 = 1e-4
  function RaySolver(env, nbeams, minangle, maxangle, ds, atol, rugocity, athreshold, solver, solvertol)
    nbeams < 0 && (nbeams = 0)
    -π/2 ≤ minangle ≤ π/2 || throw(ArgumentError("minangle should be between -π/2 and π/2"))
    -π/2 ≤ maxangle ≤ π/2 || throw(ArgumentError("maxangle should be between -π/2 and π/2"))
    minangle < maxangle || throw(ArgumentError("maxangle should be more than minangle"))
    if solver === missing
      if ssp(env) isa IsoSSP || (ssp(env) isa SampledSSP && length(ssp(env).z) ≤ 2)
        solver = Tsit5()
      else
        solver = Rodas5()
      end
    end
    new{typeof(env),typeof(solver)}(check(RaySolver, env), nbeams, minangle, maxangle, ds, atol, rugocity, athreshold, solver, solvertol)
  end
end

"""
    RaySolver(env; nbeams, minangle, maxangle, ds, atol, rugocity, athreshold, solver, solvertol)

Create a RaySolver propagation model.
"""
RaySolver(env; kwargs...) = RaySolver(; env=env, kwargs...)

### interface functions

function check(::Type{RaySolver}, env::Union{<:UnderwaterEnvironment,Missing})
  env
end

function arrivals(model::RaySolver, tx1::AcousticSource, rx1::AcousticReceiver)
  check2d([tx1], [rx1])
  a = [RayArrival(r.time, r.phasor, r.surface, r.bottom, r.launchangle, r.arrivalangle) for r ∈ eigenrays(model, tx1, rx1; ds=0.0)]
  sort(a; by = a1 -> a1.time)
end

function eigenrays(model::RaySolver, tx1::AcousticSource, rx1::AcousticReceiver; ds=model.ds)
  check2d([tx1], [rx1])
  p2 = location(rx1)
  nbeams = model.nbeams
  if nbeams == 0
    p1 = location(tx1)
    R = abs(p2[1] - p1[1])
    h = maxdepth(bathymetry(model.env))
    nbeams = ceil(Int, 16 * (model.maxangle - model.minangle) / atan(h, R))
  end
  θ = range(model.minangle, model.maxangle; length=nbeams)
  n = length(θ)
  rmax = ForwardDiff.value(p2[1])
  err = fill(NaN64, n)
  Threads.@threads for i ∈ 1:n
    p3 = traceray(model, tx1, θ[i], rmax).raypath[end]
    if isapprox(p3[1], p2[1]; atol=model.atol) && isapprox(p3[2], p2[2]; atol=model.atol)
      err[i] = ForwardDiff.value(p3[3] - p2[3])
    end
  end
  T = promote_type(envrealtype(model.env), eltype(location(tx1)), typeof(nominalfrequency(tx1)), eltype(location(rx1)))
  erays = Array{RayArrival{T,T},1}(undef, 0)
  for i ∈ 1:n
    if isapprox(err[i], 0.0; atol=model.atol)
      eray = traceray(model, tx1, convert(T, θ[i]), rmax, ds)
      push!(erays, eray)
    elseif i > 1 && !isnan(err[i-1]) && !isnan(err[i]) && sign(err[i-1]) * sign(err[i]) == -1
      a, b = ordered(θ[i-1], θ[i])
      soln = optimize(ϕ -> abs2(traceray(model, tx1, ϕ, rmax).raypath[end][3] - p2[3]), a, b; abs_tol=1e-4)
      if Optim.converged(soln) && soln.minimum ≤ 1e-2
        eray = traceray(model, tx1, convert(T, soln.minimizer), rmax, ds)
        push!(erays, eray)
      end
    elseif i > 2 && isnearzero(err[i-2], err[i-1], err[i])
      a, b = ordered(θ[i-2], θ[i])
      soln = optimize(ϕ -> abs2(traceray(model, tx1, ϕ, rmax).raypath[end][3] - p2[3]), a, b; abs_tol=1e-4)
      if Optim.converged(soln) && soln.minimum ≤ 1e-2
        eray = traceray(model, tx1, convert(T, soln.minimizer), rmax, ds)
        push!(erays, eray)
      end
    end
  end
  erays
end

function rays(model::RaySolver, tx1::AcousticSource, θ::Real, rmax)
  check2d([tx1], [])
  traceray(model, tx1, θ, rmax, model.ds)
end

function transfercoef(model::RaySolver, tx1::AcousticSource, rx::AcousticReceiverGrid2D{T1}; mode=:coherent) where {T1}
  check2d([tx1], rx)
  mode === :coherent || mode === :incoherent || throw(ArgumentError("Unknown mode :" * string(mode)))
  # implementation primarily based on ideas from COA (Computational Ocean Acoustics, 2nd ed., ch. 3)
  f = nominalfrequency(tx1)
  T = promote_type(envrealtype(model.env), eltype(location(tx1)), typeof(f), T1, ComplexF64)
  tc = zeros(T, size(rx,1), size(rx,2), Threads.nthreads())
  rmax = maximum(rx.xrange) + 0.1
  nbeams = model.nbeams
  if nbeams == 0
    p1 = location(tx1)
    R = abs(rmax - p1[1])
    h = maxdepth(bathymetry(model.env))
    nbeams = clamp(ceil(Int, 20 * (model.maxangle - model.minangle) / atan(h, R)), 100, 1000)
  end
  θ = range(model.minangle, model.maxangle; length=nbeams)
  δθ = Float64(θ.step)
  c₀ = soundspeed(ssp(model.env), location(tx1)...)
  ω = 2π * f
  G = 1 / (2π)^(1/4)
  Threads.@threads for θ1 ∈ θ
    β = (mode === :incoherent ? 2 : 1) * cos(θ1)/c₀
    traceray(model, tx1, θ1, rmax, 1.0; cb = (s1, u1, s2, u2, A₀, D₀, t₀, cₛ, kmah1, kmah2) -> begin
      r1, z1, ξ1, ζ1, t1, _, q1 = u1
      r2, z2, ξ2, ζ2, t2, _, q2 = u2
      ndx = findall(r -> r1 ≤ r < r2, rx.xrange)
      rz1 = (r1, z1)
      vlen = norm((r2-r1, z2-z1))
      rvlen = norm((ξ1, ζ1))
      tᵥ = (ξ1, ζ1) ./ rvlen
      nᵥ = (ζ1, -ξ1) ./ rvlen
      D = D₀ + (s1 + s2) / 2
      γ = fast_absorption(f, D, salinity(model.env))
      Wmax = max(q1, q2) * δθ
      ndx2 = findall(z -> min(z1, z2) - 4*Wmax ≤ z ≤ max(z1, z2) + 4*Wmax, rx.zrange)
      for j ∈ ndx
        for i ∈ ndx2
          rz = (rx.xrange[j], rx.zrange[i])
          v = rz .- rz1
          s = dot(v, tᵥ)
          n = dot(v, nᵥ)
          α = s / vlen
          t = t₀ + t1 + (t2 - t1) * α
          q = q1 + (q2 - q1) * α
          A = A₀ * γ * G * √abs(β * cₛ / (rz[1] * q)) # COA (3.76)
          kmah = round(Int, kmah1 + (kmah2 - kmah1) * α)
          A *= cis(-π/2 * kmah)                       # COA section 3.4.1 (KMAH correction)
          W = abs(q * δθ)                             # COA (3.74)
          A *= exp(-(n / W)^2)
          tc[j, i, Threads.threadid()] += mode === :coherent ? conj(A) * cis(-ω * t) : Complex(abs2(A), 0.0)
        end
      end
    end)
  end
  rv = dropdims(sum(tc; dims=3); dims=3)
  mode === :incoherent && (rv = sqrt.(rv))
  rv
end

### helper functions

ordered(a, b) = a < b ? (a, b) : (b, a)

function isnearzero(a, b, c)
  (isnan(a) || isnan(b) || isnan(c)) && return false
  sign(a) != sign(b) && return false
  sign(a) != sign(c) && return false
  abs(a) < abs(b) && return false
  abs(c) < abs(b) && return false
  return true
end

function check2d(tx, rx)
  # TODO: remove restrictions on coordinates
  all(location(tx1)[1] == 0.0 for tx1 ∈ tx) || throw(ArgumentError("RaySolver requires transmitters at (0, 0, z)"))
  all(location(tx1)[2] == 0.0 for tx1 ∈ tx) || throw(ArgumentError("RaySolver requires transmitters in the x-z plane"))
  all(location(rx1)[1] >= 0.0 for rx1 ∈ rx) || throw(ArgumentError("RaySolver requires receivers to be in the +x halfspace"))
  all(location(rx1)[2] == 0.0 for rx1 ∈ rx) || throw(ArgumentError("RaySolver requires receivers in the x-z plane"))
end

function rayeqns!(du, u, params, s)
  # implementation based on COA (3.161-164, 3.58-63)
  # assumes range-independent soundspeed
  r, z, ξ, ζ, t, p, q = u
  c, ∂c, ∂²c = params
  cᵥ = c(z)
  cᵥ² = cᵥ * cᵥ
  c̄ = ∂²c(z) * ξ * ξ
  du[1] = cᵥ * ξ
  du[2] = cᵥ * ζ
  du[3] = 0
  du[4] = -∂c(z) / cᵥ²
  du[5] = 1 / cᵥ
  du[6] = -c̄ * q
  du[7] = cᵥ * p
end

function checkray!(out, u, s, integrator, a::Altimetry, b::Bathymetry, rmax)
  out[1] = altitude(a, u[1], 0.0) - u[2]   # surface reflection
  out[2] = u[2] + depth(b, u[1], 0.0)      # bottom reflection
  out[3] = rmax - u[1]                     # maximum range
  out[4] = u[3]                            # ray turned back
end

function traceray1(T, model, r0, z0, θ, rmax, ds, p0, q0)
  a = altimetry(model.env)
  b = bathymetry(model.env)
  c = z -> soundspeed(ssp(model.env), 0.0, 0.0, z)
  ∂c = z -> ForwardDiff.derivative(c, z)
  ∂²c = z -> ForwardDiff.derivative(∂c, z)
  cᵥ = c(z0)
  u0 = [convert(T, r0), z0, cos(θ)/cᵥ, sin(θ)/cᵥ, zero(T), p0, q0]
  prob = ODEProblem{true}(rayeqns!, u0, (zero(T), one(T) * model.rugocity * (rmax-r0)/cos(θ)), (c, ∂c, ∂²c))
  cb = VectorContinuousCallback(
    (out, u, s, i) -> checkray!(out, u, s, i, a, b, rmax),
    (i, ndx) -> terminate!(i), 4; rootfind=true)
  if ds ≤ 0
    soln = solve(prob, model.solver; abstol=model.solvertol, save_everystep=false, callback=cb)
  else
    soln = solve(prob, model.solver; abstol=model.solvertol, saveat=ds, callback=cb)
  end
  s2 = soln[end]
  soln.t[end], s2[1], s2[2], atan(s2[4], s2[3]), s2[5], s2[6], s2[7], soln.u, soln.t
end

function traceray(model::RaySolver, tx1::AcousticSource, θ::Real, rmax, ds=0.0; cb=nothing)
  θ₀ = θ
  f = nominalfrequency(tx1)
  zmin = -maxdepth(bathymetry(model.env))
  ϵ = one(zmin) * model.solvertol
  p = location(tx1)
  c₀ = soundspeed(ssp(model.env), p...)
  T = promote_type(envrealtype(model.env), eltype(p), typeof(f), typeof(θ), typeof(rmax))
  raypath = Array{NTuple{3,T}}(undef, 0)
  A = one(Complex{T})   # phasor
  s = 0                 # surface bounces
  b = 0                 # bottom bounces
  t = zero(T)           # time along ray
  D = zero(T)           # distance along ray
  q = zero(T)           # spreading factor
  qp = one(T) / c₀      # spreading rate
  while true
    dD, r, z, θ, dt, qp, q, u, svec = traceray1(T, model, p[1], p[3], θ, rmax, ds, qp, q)
    oq = u[1][7]
    kmah = 0
    kmahhist = zeros(length(u))
    for i ∈ 1:length(u)
      sign(oq) * sign(u[i][7]) == -1 && (kmah += 1)
      kmahhist[i] = kmah
      u[i][7] != 0.0 && (oq = u[i][7])
      push!(raypath, (u[i][1], 0.0, u[i][2]))
    end
    if cb !== nothing
      for i ∈ 2:length(u)
        cₛ = soundspeed(ssp(model.env), (u[i-1][1]+u[i][1])/2, 0.0, (u[i-1][2]+u[i][2])/2)
        cb(svec[i-1], u[i-1], svec[i], u[i], A, D, t, cₛ, kmahhist[i-1], kmahhist[i])
      end
    end
    t += dt
    D += dD
    A *= cis(-π/2 * kmah)                 # COA section 3.4.1 (KMAH correction)
    zmin = -depth(bathymetry(model.env), r, 0.0)
    zmax = altitude(altimetry(model.env), r, 0.0)
    p = (r, 0.0, clamp(z, zmin+ϵ, zmax-ϵ))
    if r ≥ rmax - 1e-3
      break
    elseif isapprox(z, zmax; atol=1e-3)   # hit the surface
      s += 1
      A *= reflectioncoef(seasurface(model.env), f, π/2 - θ)
      α = atan(ForwardDiff.derivative(x -> altitude(altimetry(model.env), x, 0.0), r))
      θ = -θ + 2α
    elseif isapprox(z, zmin; atol=1e-3)   # hit the bottom
      b += 1
      A *= reflectioncoef(seabed(model.env), f, π/2 + θ)
      α = atan(ForwardDiff.derivative(x -> depth(bathymetry(model.env), x, 0.0), r))
      θ = -θ - 2α
    else
      break
    end
    if abs(θ) ≥ π/2
      break
    end
    if abs(A)/q < model.athreshold
      break
    end
  end
  cₛ = soundspeed(ssp(model.env), p...)
  A *= √abs(cₛ * cos(θ₀) / (p[1] * c₀ * q))        # COA (3.65)
  A *= fast_absorption(f, D, salinity(model.env))
  RayArrival{T,T}(t, conj(A), s, b, θ₀, -θ, raypath)    # conj(A) needed to match with Bellhop
end
