using Statistics

export Bellhop

struct Bellhop{T} <: PropagationModel{T}
  env::T
  nbeams::Int
  minangle::Float32   # degrees
  maxangle::Float32   # degrees
  function Bellhop(env, nbeams, minangle, maxangle)
    nbeams < 0 && (nbeams = 0)
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
    altimetry(env) isa FlatSurface || throw(ArgumentError("Non-flat altimetry not yet supported"))
    seabed(env) isa Rayleigh || throw(ArgumentError("Seabed type not supported"))
    seasurface(env) === Vacuum || throw(ArgumentError("Only vacuum seasurface supported"))
  end
  env
end

function arrivals(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver)
  mktempdir(prefix="bellhop_") do dirname
    writeenv(model, [tx1], [rx1], "A", dirname)
    bellhop(dirname)
    readarrivals(joinpath(dirname, "model.arr"))
  end
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
    readshd(joinpath(dirname, "model.shd"))
  end
end

transfercoef(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent) = transfercoef(model, tx1, [rx1]; mode=mode)

function eigenrays(model::Bellhop, tx1::AcousticSource, rx1::AcousticReceiver)
  mktempdir(prefix="bellhop_") do dirname
    writeenv(model, [tx1], [rx1], "E", dirname)
    bellhop(dirname)
    readrays(joinpath(dirname, "model.ray"))
  end
end

### helper functions

function bellhop(dirname)
  try
    infilebase = joinpath(dirname, "model")
    outfilename = joinpath(dirname, "output.txt")
    run(pipeline(ignorestatus(`bellhop.exe $infilebase`); stdout=outfilename, stderr=outfilename))
  catch
    error("Unable to execute bellhop.exe")
  end
end

function writeenv(model::Bellhop, tx::Vector{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver}, taskcode, dirname)
  env = model.env
  name = split(basename(dirname), "_")[end]
  filename = joinpath(dirname, "model.env")
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
    # TODO: support altimetry
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
    c2 = soundspeed(ss, 0.0, 0.0, -depth) * bed.cᵣ
    α = bed.δ / (c2/(f/1000)) * 40π / log(10)       # based on APL-UW TR 9407 (1994), IV-8 equation (4)
    @printf(io, "%0.6f %0.6f 0.0 %0.6f %0.6f /\n", depth, c2, bed.ρᵣ, α)
    for i ∈ 1:length(tx)
      p = location(tx[i])
      (p[1] == 0.0 && p[2] == 0.0) || throw(ArgumentError("Transmitters must be located at (0, 0)"))
    end
    printarray(io, [-location(tx1)[3] for tx1 ∈ tx])
    if length(rx) == 1
      printarray(io, [-location(rx[1])[3]])
      maxr = sqrt(sum(abs2, location(rx[1])[1:2])) / 1000.0
      printarray(io, [maxr])
    elseif rx isa AcousticReceiverGrid2D
      printarray(io, -rx.zrange)
      printarray(io, rx.xrange ./ 1000.0)
      maxr = maximum(rx.xrange) / 1000.0
    else
      throw(ArgumentError("Receivers must be on a 2D grid"))
    end
    # TODO: support source directionality
    println(io, "'", taskcode, "'")
    @printf(io, "%d\n", model.nbeams)
    @printf(io, "%0.6f %0.6f /\n", model.minangle, model.maxangle)
    @printf(io, "0.0 %0.6f %0.6f\n", 1.01*depth, 1.01*maxr)
  end
end

function printarray(io, a::AbstractVector)
  println(io, length(a))
  for a1 ∈ a
    @printf(io, "%0.6f ", a1)
  end
  println(io, "/")
end

function readrays(filename)
  rays = Arrival{Missing,Missing,Float64,Float64}[]
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
      push!(rays, Arrival(missing, missing, sb, bb, -deg2rad(aod), NaN64, raypath))
    end
  end
  rays
end

function readarrivals(filename)
  arrivals = Arrival{Float64,ComplexF64,Float64,Missing}[]
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
            push!(arrivals, Arrival(t, A * cis(deg2rad(ph)), sb, bb, -deg2rad(aod), deg2rad(aoa)))
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
      pressure[:,ird+1] .= temp
    end
    pressure
  end
end
