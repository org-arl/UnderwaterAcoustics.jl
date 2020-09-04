export Bellhop

"""
$(TYPEDEF)
A propagation model based on an external FORTRAN Bellhop executable.
"""
struct Bellhop{T} <: PropagationModel{T}
  env::T
  nbeams::Int
  minangle::Float32
  maxangle::Float32
  debug::Bool
  function Bellhop(env, nbeams, minangle, maxangle, debug)
    nbeams < 0 && (nbeams = 0)
    -π/2 ≤ minangle ≤ π/2 || throw(ArgumentError("minangle should be between -π/2 and π/2"))
    -π/2 ≤ maxangle ≤ π/2 || throw(ArgumentError("maxangle should be between -π/2 and π/2"))
    minangle < maxangle || throw(ArgumentError("maxangle should be more than minangle"))
    new{typeof(env)}(check(Bellhop, env), nbeams, Float32(minangle), Float32(maxangle), debug)
  end
end

"""
    Bellhop(env; debug=false)
    Bellhop(env, nbeams, minangle, maxangle, debug)

Create a Bellhop propagation model.
"""
Bellhop(env; debug=false) = Bellhop(env, 0, -80°, 80°, debug)

### interface functions

function check(::Type{Bellhop}, env::Union{<:UnderwaterEnvironment,Missing})
  if env === missing
    mktempdir(prefix="bellhop_") do dirname
      try
        bellhop(dirname, false)
      catch e
        e isa BellhopError && e.details == ["Unable to execute bellhop.exe"] && throw(e)
      end
    end
  else
    seabed(env) isa Rayleigh || throw(ArgumentError("Seabed type not supported"))
    seasurface(env) === Vacuum || throw(ArgumentError("Only vacuum seasurface supported"))
  end
  env
end

function arrivals(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver)
  mktempdir(prefix="bellhop_") do dirname
    writeenv(model, [tx1], [rx1], "A", dirname)
    bellhop(dirname, model.debug)
    readarrivals(joinpath(dirname, "model.arr"))
  end
end

function transfercoef(model::Bellhop, tx1::AcousticSource, rx::AcousticReceiverGrid2D; mode=:coherent)
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
    bellhop(dirname, model.debug)
    readshd(joinpath(dirname, "model.shd"))
  end
end

function transfercoef(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
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
    writeenv(model, [tx1], [rx1], taskcode, dirname)
    bellhop(dirname, model.debug)
    readshd(joinpath(dirname, "model.shd"))[1]
  end
end

function eigenrays(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver)
  mktempdir(prefix="bellhop_") do dirname
    writeenv(model, [tx1], [rx1], "E", dirname)
    bellhop(dirname, model.debug)
    readrays(joinpath(dirname, "model.ray"))
  end
end

function rays(model::Bellhop, tx1::AcousticSource, θ::AbstractArray, rmax)
  θ isa AbstractRange || length(θ) == 1 || throw(ArgumentError("Bellhop only supports uniformly spaced angles"))
  all(-π/2 .< θ .< π/2) || throw(ArgumentError("θ must be between -π/2 and π/2"))
  mktempdir(prefix="bellhop_") do dirname
    writeenv(model, [tx1], [AcousticReceiver(rmax, 0.0)], "R", dirname; minangle=minimum(θ), maxangle=maximum(θ), nbeams=length(θ))
    bellhop(dirname, model.debug)
    readrays(joinpath(dirname, "model.ray"))
  end
end

rays(model::Bellhop, tx1::AcousticSource, θ, rmax) = rays(model, tx1, [θ], rmax)[1]

### helper functions

struct BellhopError <: Exception
  details::Vector{String}
end

function Base.show(io::IO, e::BellhopError)
  if length(e.details) == 1
    println(io, e.details[1])
  else
    println(io, "Bellhop said:")
    for s ∈ e.details
      println(io, "  ", s)
    end
  end
end

function bellhop(dirname, debug)
  infilebase = joinpath(dirname, "model")
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(ignorestatus(`bellhop.exe $infilebase`); stdout=outfilename, stderr=outfilename))
    if debug
      @info "Bellhop run completed in $dirname, press ENTER to delete intermediate files..."
      readline()
    end
  catch
    throw(BellhopError(["Unable to execute bellhop.exe"]))
  end
  err = String[]
  checkerr!(err, outfilename)
  checkerr!(err, joinpath(dirname, "model.prt"))
  if length(err) > 0
    throw(BellhopError(err))
  end
end

function checkerr!(err, filename)
  output = false
  open(filename) do f
    for s in eachline(f)
      if output || occursin("ERROR", uppercase(s))
        push!(err, s)
        output = true
      end
    end
  end
end

function writeenv(model::Bellhop, tx::Vector{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver}, taskcode, dirname; minangle=model.minangle, maxangle=model.maxangle, nbeams=model.nbeams)
  all(location(tx1)[1] == 0.0 for tx1 ∈ tx) || throw(ArgumentError("Bellhop requires transmitters at (0, 0, z)"))
  all(location(tx1)[2] == 0.0 for tx1 ∈ tx) || throw(ArgumentError("Bellhop 2D requires transmitters in the x-z plane"))
  all(location(rx1)[1] >= 0.0 for rx1 ∈ rx) || throw(ArgumentError("Bellhop requires receivers to be in the +x halfspace"))
  all(location(rx1)[2] == 0.0 for rx1 ∈ rx) || throw(ArgumentError("Bellhop 2D requires receivers in the x-z plane"))
  env = model.env
  name = split(basename(dirname), "_")[end]
  filename = joinpath(dirname, "model.env")
  open(filename, "w") do io
    println(io, "'", name, "'")
    flist = [nominalfrequency(tx1) for tx1 ∈ tx]
    f = sum(flist) / length(flist)
    maximum(abs.(flist .- f))/f > 0.2 && @warn("Source frequency varies by more than 20% from nominal frequency")
    @printf(io, "%0.6f\n", f)
    println(io, "1")
    if length(rx) == 1
      maxr = location(rx[1])[1]
    elseif rx isa AcousticReceiverGrid2D
      maxr = maximum(rx.xrange)
    else
      throw(ArgumentError("Receivers must be on a 2D grid"))
    end
    ss = ssp(env)
    sspi = "S"
    ss isa SampledSSP && ss.interp === :linear && (sspi = "C")
    print(io, "'", sspi, "VWT")
    alt = altimetry(env)
    if !(alt isa FlatSurface)
      print(io, "*")
      createadfile(joinpath(dirname, "model.ati"), alt, (p...) -> -altitude(p...), maxr, f)
    end
    println(io, "'")
    bathy = bathymetry(env)
    waterdepth = maxdepth(bathy)
    @printf(io, "1 0.0 %0.6f\n", waterdepth)
    if ss isa IsoSSP
      @printf(io, "0.0 %0.6f /\n", soundspeed(ss, 0.0, 0.0, 0.0), )
      @printf(io, "%0.6f %0.6f /\n", waterdepth, soundspeed(ss, 0.0, 0.0, 0.0))
    elseif ss isa SampledSSP
      for i ∈ 1:length(ss.z)
        @printf(io, "%0.6f %0.6f /\n", -ss.z[i], ss.c[i])
      end
    else
      for d ∈ range(0.0, waterdepth; length=recommendlength(waterdepth, f))
        @printf(io, "%0.6f %0.6f /\n", d, soundspeed(ss, 0.0, 0.0, -d))
      end
      floor(waterdepth) != waterdepth && @printf(io, "%0.6f %0.6f /\n", waterdepth, soundspeed(ss, 0.0, 0.0, -waterdepth))
    end
    print(io, "'A")
    if !(bathy isa ConstantDepth)
      print(io, "*")
      createadfile(joinpath(dirname, "model.bty"), bathy, depth, maxr, f)
    end
    println(io, "' 0.0") # bottom roughness = 0
    bed = seabed(env)
    c2 = soundspeed(ss, 0.0, 0.0, -waterdepth) * bed.cᵣ
    α = bed.δ * 40π / log(10)      # based on APL-UW TR 9407 (1994), IV-9 equation (4)
    @printf(io, "%0.6f %0.6f 0.0 %0.6f %0.6f /\n", waterdepth, c2, bed.ρᵣ, α)
    printarray(io, [-location(tx1)[3] for tx1 ∈ tx])
    if length(rx) == 1
      printarray(io, [-location(rx[1])[3]])
      printarray(io, [maxr / 1000.0])
    elseif rx isa AcousticReceiverGrid2D
      printarray(io, -rx.zrange)
      printarray(io, rx.xrange ./ 1000.0)
    end
    println(io, "'", taskcode, "'")
    @printf(io, "%d\n", nbeams)
    @printf(io, "%0.6f %0.6f /\n", rad2deg(minangle), rad2deg(maxangle))
    @printf(io, "0.0 %0.6f %0.6f\n", 1.01*waterdepth, 1.01*maxr)
  end
end

function printarray(io, a::AbstractVector)
  println(io, length(a))
  for a1 ∈ a
    @printf(io, "%0.6f ", a1)
  end
  println(io, "/")
end

function recommendlength(x, f)
  # recommendation based on nominal half-wavelength spacing
  λ = 1500.0 / f
  clamp(round(Int, 2x / λ) + 1, 25, 1000)
end

function createadfile(filename, data, func, maxr, f)
  open(filename, "w") do io
    interp = "L"
    if data isa SampledDepth || data isa SampledAltitude
      x = data.x
      data.interp !== :linear && (interp = "C")
    else
      x = range(0.0, maxr; length=recommendlength(maxr, f))
    end
    println(io, "'", interp, "'")
    println(io, length(x))
    for i ∈ 1:length(x)
      @printf(io, "%0.6f %0.6f\n", x[i]/1000.0, func(data, x[i], 0.0))
    end
  end
end

function readrays(filename)
  rays = RayArrival{Float64,Float64}[]
  open(filename, "r") do io
    [readline(io) for i ∈ 1:7]
    while !eof(io)
      s = strip(readline(io))
      length(s) == 0 && break
      aod = parse(Float64, s)
      pts, sb, bb = parse.(Int, split(strip(readline(io)) ,r" +"))
      raypath = Array{NTuple{3,Float64}}(undef, pts)
      for k ∈ 1:pts
        x, d = parse.(Float64, split(strip(readline(io)) ,r" +"))
        raypath[k] = (x, 0.0, -d)
      end
      push!(rays, RayArrival(NaN64, NaN64, sb, bb, -deg2rad(aod), NaN64, raypath))
    end
  end
  rays
end

function readarrivals(filename)
  arrivals = RayArrival{Float64,Missing}[]
  open(filename, "r") do io
    s = strip(readline(io))
    if occursin("2D", s)
      f = parse(Float64, strip(readline(io)))
      v = split(strip(readline(io)) ,r" +")
      n = parse(Int, v[1])
      txdepth = parse.(Float64, v[2:end])
      n == length(txdepth) || error("Wrong number of txdepth entries in arrivals")
      v = split(strip(readline(io)) ,r" +")
      n = parse(Int, v[1])
      rxdepth = parse.(Float64, v[2:end])
      n == length(rxdepth) || error("Wrong number of rxdepth entries in arrivals")
      v = split(strip(readline(io)) ,r" +")
      n = parse(Int, v[1])
      rxrange = parse.(Float64, v[2:end])
      n == length(rxrange) || error("Wrong number of rxrange entries in arrivals")
    else
      v = split(s ,r" +")
      f = parse(Float64, v[1])
      n1, n2, n3 = parse.(Int, v[2:4])
      txdepth = parse.(Float64, split(strip(readline(io)) ,r" +"))
      rxdepth = parse.(Float64, split(strip(readline(io)) ,r" +"))
      rxrange = parse.(Float64, split(strip(readline(io)) ,r" +"))
      n1 == length(txdepth) || error("Wrong number of txdepth entries in arrivals")
      n2 == length(rxdepth) || error("Wrong number of rxdepth entries in arrivals")
      n3 == length(rxrange) || error("Wrong number of rxrange entries in arrivals")
    end
    for j ∈ 1:length(txdepth)
      readline(io)
      for k ∈ 1:length(rxdepth)
        for m ∈ 1:length(rxrange)
          count = parse(Int, strip(readline(io)))
          for n ∈ 1:count
            v = split(strip(readline(io)) ,r" +")
            length(v) == 8 || error("Wrong number of data entries in arrivals")
            A, ph, t, _, aod, aoa = parse.(Float64, v[1:6])
            sb, bb = parse.(Int, v[7:8])
            push!(arrivals, RayArrival(t, A * cis(deg2rad(ph)), sb, bb, -deg2rad(aod), deg2rad(aoa)))
          end
        end
      end
    end
  end
  sort(arrivals; by = a -> a.time)
end

function readshd(filename)
  open(filename, "r") do io
    r = read(io, UInt32)
    seek(io, 4r)
    b = Array{UInt8}(undef, 10)
    read!(io, b)
    strip(String(b)) == "rectilin" || error("Bad shd file format: incorrect ptype")
    seek(io, 8r)
    nfreq = read(io, UInt32)
    nfreq == 1 || error("Bad shd file format: incorrect nfreq")
    nθ = read(io, UInt32)
    nθ == 1 || error("Bad shd file format: incorrect nθ")
    nsx = read(io, UInt32)
    nsy = read(io, UInt32)
    nsd = read(io, UInt32)
    nsd == 1 || error("Bad shd file format: incorrect nsd")
    nrd = read(io, UInt32)
    nrr = read(io, UInt32)
    pressure = Array{ComplexF32}(undef, nrr, nrd)
    for ird ∈ 0:nrd-1
      recnum = 10 + ird
      seek(io, recnum * 4r)
      temp = Array{ComplexF32}(undef, nrr)
      read!(io, temp)
      pressure[:,ird+1] .= -temp    # negative because Bellhop seems to have a 180° phase inversion
    end
    pressure
  end
end
