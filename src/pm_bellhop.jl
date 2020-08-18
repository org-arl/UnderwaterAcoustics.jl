using Statistics
using Printf

export Bellhop

struct Bellhop{T} <: PropagationModel{T}
  env::T
  nbeams::Int
  minangle::Float32   # degrees
  maxangle::Float32   # degrees
  function Bellhop(env, nbeams, minangle, maxangle)
    nbeams > 0 || throw(ArgumentError("nbeams should be more than 0"))
    -90 ≤ minangle ≤ 90 || throw(ArgumentError("minangle should be between -90 and 90"))
    -90 ≤ maxangle ≤ 90 || throw(ArgumentError("maxangle should be between -90 and 90"))
    minangle < maxangle || throw(ArgumentError("maxangle should be more than minangle"))
    new{typeof(env)}(check(Bellhop, env), nbeams, Float32(minangle), Float32(maxangle))
  end
end

Bellhop(env) = Bellhop(env, 0, -80, 80)

### interface functions

function check(::Type{Bellhop}, env::Union{<:UnderwaterEnvironment,Missing})
  if env === missing
    mktempdir(prefix="bellhop_") do dirname
      bellhop(dirname)
    end
  else
    # TODO
  end
  env
end

function arrivals(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver)

end

function transfercoef(model::Bellhop, tx1::AcousticSource, rx::AbstractArray{<:AcousticReceiver}; mode=:coherent)
  if mode === :coherent
    taskcode = "C"
  elseif mode === :incoherent
    taskcode = "I"
  elseif mode === :semicoherent
    taskcode = "S"
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  mktempdir(prefix="bellhop_") do dirname
    writeenv(model, [tx1], rx, taskcode, dirname)
    bellhop(dirname)

  end
end

transfercoef(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent) = transfercoef(model, tx1, [rx1]; mode=mode)

function eigenrays(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver)

end

### helper functions

function bellhop(dirname)
  try
    infilename = joinpath(dirname, "input.env")
    outfilename = joinpath(dirname, "output.txt")
    run(pipeline(ignorestatus(`bellhop.exe $infilename`); stdout=outfilename, stderr=outfilename))
  catch
    error("Unable to execute bellhop.exe")
  end
end

function writeenv(model::Bellhop, tx::Vector{<:AcousticSource}, rx::Vector{<:AcousticReceiver}, taskcode, dirname)
  env = model.env
  name = split(basename(dirname), "_")[end]
  filename = joinpath(dirname, "input.env")
  open(filename, "w") do io
    println(io, "'", name, "'")
    flist = [nominalfrequency(tx1) for tx1 ∈ tx]
    f = mean(flist)
    maximum(abs.(flist .- f))/f > 0.2 && @warn("Source frequency varies by more than 20% from nominal frequency")
    @printf(io, "%0.6f\n", f)
    println(io, "1")
    ss = ssp(env)
    sspi = "C"
    ss isa SampledSSP && ss.interp == :cubic && (sspi = "S")
    print(io, "'", sspi, "VWT")
    altimetry(env) isa FlatSurface || throw(ArgumentError("Non-flat altimetry not yet supported"))  # TODO: support altimetry
    println(io, "'")
    bathy = bathymetry(env)
    depth = maxdepth(bathy)
    @printf(io, "1 0.0 %0.6f\n", depth)
    # TODO: support range-dependent soundspeed
    if ss isa IsoSSP
      @printf(io, "0.0 %0.6f /\n", soundspeed(ss, 0.0, 0.0, 0.0), )
      @printf(io, "%0.6f %0.6f /\n", depth, soundspeed(ss, 0.0, 0.0, 0.0))
    elseif ss isa SampledSSP
      for i ∈ 1:length(ss.z)
        @printf(io, "%0.6f %0.6f /\n", -ss.z[i], ss.c[i])
      end
    else
      for d ∈ 0.0:1.0:depth
        @printf(io, "%0.6f %0.6f /\n", d, soundspeed(ss, 0.0, 0.0, -d))
      end
      floor(depth) != depth && @printf(io, "%0.6f %0.6f /\n", depth, soundspeed(ss, 0.0, 0.0, -depth))
    end
    # TODO: support bathymetry
    # TODO: support bottom roughness
    println(io, "'A' 0.0")
    bed = seabed(env)
    if bed isa Rayleigh
      # FIXME: handle relative density / soundspeed better
      @printf(io, "%0.6f %0.6f 0.0 %0.6f %0.6f /\n", depth, 1500 * bed.cᵣ, bed.ρᵣ, bed.δ)  # TODO: check units for absorption
    else
      throw(ArgumentError("Seabed type not supported"))
    end
    for i ∈ 1:length(tx)
      p = location(tx[i])
      (p[1] == 0.0 && p[2] == 0.0) || throw(ArgumentError("Transmitters must be located at (0, 0)"))
    end
    printarray(io, [-location(tx1)[3] for tx1 ∈ tx])
    printarray(io, [-location(rx1)[3] for rx1 ∈ rx])
    r = [sqrt(sum(abs2, location(rx1)[1:2]))/1000.0 for rx1 ∈ rx]
    printarray(io, r)
    # TODO: support source directionality
    println(io, "'", taskcode, "'")
    @printf(io, "%d\n", model.nbeams)
    @printf(io, "%0.6f %0.6f /\n", model.minangle, model.maxangle)
    @printf(io, "0.0 %0.6f %0.6f\n", 1.01*depth, 1.01*maximum(r))
  end
end

function printarray(io, a::AbstractVector)
  println(io, length(a))
  for a1 ∈ a
    @printf(io, "%0.6f ", a1)
  end
  println(io, "/")
end
