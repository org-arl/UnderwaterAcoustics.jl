export Kraken

"""
$(TYPEDEF)
A propagation model based on an external FORTRAN Kraken executable.
"""
struct Kraken{T} <: PropagationModel{T}
  env::T
  nmodes::Int
  nmedia::Int    # number of media (<20)
  nmesh::Int     # set to 0 to allow Kraken to automatically calculate
  clow::Float64  # lower phase speed limit
  chigh::Float64 # upper phase speed limit
  debug::Bool
  function Kraken(env, nmodes, nmedia, nmesh, clow, chigh,  debug)
    nmodes > 1 || throw(ArgumentError("number of modes should be a postive integer"))
    nmedia < 20 || throw(ArgumentError("number of media is recommended to be smaller than 20"))
    new{typeof(env)}(check(Kraken, env), nmodes, nmedia, nmesh, clow, chigh, debug)
  end
end

"""
    Kraken(env; debug=false)
    Kraken(env, chigh; debug=false)
    Kraken(env, nmodes, nmedia, nmesh, clow, chigh, debug)

Create a Kraken propagation model.
"""
Kraken(env; debug=false) = Kraken(env, 9999, 1, 0,  0.0, 1600.0, debug)
Kraken(env, chigh; debug=false) = Kraken(env, 9999, 1, 0,  0.0, chigh, debug)

### interface functions

function check(::Type{Kraken}, env::Union{<:UnderwaterEnvironment,Missing})
  if env === missing
    mktempdir(prefix="kraken_") do dirname
      try
        kraken(dirname, false)
      catch e
        e isa KrakenError && e.details == ["Unable to execute kraken.exe"] && throw(e)
      end
    end
  else
    seabed(env) isa RayleighReflectionCoef || throw(ArgumentError("Seabed type not supported"))
    seasurface(env) === Vacuum || throw(ArgumentError("Only vacuum seasurface supported"))  #TODO: incoperate other surface types 
  end
  env
end


function transfercoef(model::Kraken, tx1::AcousticSource, rx::AcousticReceiverGrid2D; mode=:coherent)
  if mode === :coherent
    taskcode = "C"
  elseif mode === :incoherent
    throw(ArgumentError("Incoherent mode addition is not yet supported")) #TODO: incoherent mode addition
    taskcode = "I"
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  mktempdir(prefix="kraken_") do dirname
    xrev, zrev = writeenv(model, [tx1], rx, dirname)
    writeflp(model, [tx1], rx, dirname)
    kraken(dirname, model.debug)
    field(dirname, model.debug)
    readshd_kraken(joinpath(dirname, "model.shd"); xrev=xrev, zrev=zrev)
  end
end

function transfercoef(model::Kraken, tx1::AcousticSource, rx1::AcousticReceiver; mode=:coherent)
  if mode === :coherent
    taskcode = "C"
  elseif mode === :incoherent
    throw(ArgumentError("Incoherent mode addition is not yet supported")) #TODO: incoherent mode addition
    taskcode = "I"
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  mktempdir(prefix="kraken_") do dirname
    writeenv(model, [tx1], [rx1], dirname)
    writeflp(model, [tx1], [rx1], dirname)
    kraken(dirname, model.debug)
    field(dirname, model.debug)
    readshd_kraken(joinpath(dirname, "model.shd"))[1]
  end
end

### helper functions

struct KrakenError <: Exception
  details::Vector{String}
end

function Base.show(io::IO, e::KrakenError)
  if length(e.details) == 1
    println(io, e.details[1])
  else
    println(io, "Kraken said:")
    for s ∈ e.details
      println(io, "  ", s)
    end
  end
end

function kraken(dirname, debug)
  infilebase = joinpath(dirname, "model")
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(ignorestatus(`kraken.exe $infilebase`); stdout=outfilename, stderr=outfilename))
    if debug
      @info "Kraken run completed in $dirname, press ENTER to delete intermediate files..."
      readline()
    end
  catch
    throw(KrakenError(["Unable to execute kraken.exe"]))
  end
  err = String[]
  checkerr_kraken!(err, outfilename)
  checkerr_kraken!(err, joinpath(dirname, "model.prt"))
  if length(err) > 0
    throw(KrakenError(err))
  end
end

function field(dirname, debug)
  infilebase = joinpath(dirname, "model")
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(ignorestatus(Cmd(`field.exe $infilebase`; dir=dirname)); stdout=outfilename, stderr=outfilename))
    if debug
      @info "Field run completed in $dirname, press ENTER to delete intermediate files..."
      readline()
    end
  catch
    throw(KrakenError(["Unable to execute field.exe"]))
  end
  err = String[]
  checkerr_kraken!(err, outfilename)
  checkerr_kraken!(err, joinpath(dirname, "field.prt"))
  if length(err) > 0
    throw(KrakenError(err))
  end
end

function checkerr_kraken!(err, filename)
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

function print_distance(io, a::AbstractVector) #prints start and end values in an array
  println(io, length(a))
  @printf(io, "%0.6f ", a[1])
  (length(a) > 1) && @printf(io, "%0.6f ", a[end])
  println(io, "/")
end

# Implementation of env file is based on Kraken user manual section 4.2.2 on p79:
function writeenv(model::Kraken, tx::Vector{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver}, dirname; nmedia=model.nmedia, nmesh = model.nmesh, clow = model.clow, chigh = model.chigh)
  all(location(tx1)[1] == 0.0 for tx1 ∈ tx) || throw(ArgumentError("Kraken requires transmitters at (0, 0, z)")) 
  all(location(tx1)[2] == 0.0 for tx1 ∈ tx) || throw(ArgumentError("Kraken 2D requires transmitters in the x-z plane"))
  all(location(rx1)[1] >= 0.0 for rx1 ∈ rx) || throw(ArgumentError("Kraken requires receivers to be in the +x halfspace"))
  all(location(rx1)[2] == 0.0 for rx1 ∈ rx) || throw(ArgumentError("Kraken 2D requires receivers in the x-z plane"))
  xrev = false
  zrev = false
  env = model.env
  name = split(basename(dirname), "_")[end]
  filename = joinpath(dirname, "model.env")
  open(filename, "w") do io
    println(io, "'", name, "'")                   # 1) title                
    flist = [nominalfrequency(tx1) for tx1 ∈ tx]    
    f = sum(flist) / length(flist)
    maximum(abs.(flist .- f))/f > 0.2 && @warn("Source frequency varies by more than 20% from nominal frequency")
    @printf(io, "%0.6f\n", f) # 2) frequency (Hz)
    @printf(io, "%i\n", nmedia)    # 3) number of media
    if length(rx) == 1 
      maxr = location(rx[1])[1]
    elseif rx isa AcousticReceiverGrid2D
      maxr = maximum(rx.xrange)
    else
      throw(ArgumentError("Receivers must be on a 2D grid"))
    end
    ss = ssp(env) 
    sspi = "S" 
    ss isa SampledSSP1D && ss.interp === :linear && (sspi = "C")
    print(io, "'", sspi, "VWT'\n")   #4) options  #TODO: allow options for other halfspace properties
    bathy = bathymetry(env)
    waterdepth = maxdepth(bathy)
    @printf(io, "%i 0.0 %0.6f\n", nmesh, waterdepth) # 5) medium info
    if ss isa IsoSSP
      @printf(io, "0.0 %0.6f /\n", soundspeed(ss, 0.0, 0.0, 0.0), )   # 5a) SSP
      @printf(io, "%0.6f %0.6f /\n", waterdepth, soundspeed(ss, 0.0, 0.0, 0.0))
    elseif ss isa SampledSSP1D
      (abs(ss.z[1]) == 0.0) || (throw(ArgumentError("user defined ssp should start at depth = 0")))
      (abs(ss.z[end]) == waterdepth)  || (throw(ArgumentError("end poiont of user defined ssp should not be smaller than water depth")))
      for i ∈ 1:length(ss.z)
        @printf(io, "%0.6f %0.6f /\n", -ss.z[i], ss.c[i])
      end
    else
      for d ∈ range(0.0, waterdepth; length=recommendlength(waterdepth, f))
        @printf(io, "%0.6f %0.6f /\n", d, soundspeed(ss, 0.0, 0.0, -d))
      end
      floor(waterdepth) != waterdepth && @printf(io, "%0.6f %0.6f /\n", waterdepth, soundspeed(ss, 0.0, 0.0, -waterdepth))
    end
    print(io, "'A")  # 6) Bottom boundary condition #TODO: allow other options
    println(io, "' 0.0") # bottom roughness = 0
    bed = seabed(env)
    c2 = soundspeed(ss, 0.0, 0.0, -waterdepth) * bed.cᵣ
    α = bed.δ * 40π / log(10)      # based on APL-UW TR 9407 (1994), IV-9 equation (4)
    @printf(io, "%0.6f %0.6f 0.0 %0.6f %0.6f 0.0/\n", waterdepth, c2, bed.ρᵣ, α) #additional line because "A" in 6)
    if !(bathy isa ConstantDepth)  #TODO: Coupled mode
      throw(ArgumentError("only range-independent batheymetry is allowed"))
    end
    @printf(io, "%0.6f  %0.6f\n",clow, chigh)  #7) phase speed limit: CLOW and CHIGH (m/s) 
    @printf(io,  "%0.6f\n", maxr / 1000.0) # 8) maximum range (km)
    print_distance(io, [-location(tx1)[3] for tx1 ∈ tx])  # 9) number of source depth, source depths (m)
    if length(rx) == 1
      print_distance(io, [-location(rx[1])[3]]) # 9) number of receiver depth, receiver depths (m)
    elseif rx isa AcousticReceiverGrid2D
      d = reverse(-rx.zrange)
      if first(d) > last(d)
       d = reverse(d)
       zrev = true
      end
      r = rx.xrange ./ 1000.0
      if first(r) > last(r)
       r = reverse(r)
       xrev = true
      end
      print_distance(io, d) # 9) receiver depths (m)
    end
  end
  xrev, zrev
end


#Implementation of .flp file is based on Kraken manual section 4.3.1 on p93
function writeflp(model::Kraken, tx::Vector{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver}, dirname; nmodes=model.nmodes)
  filename = joinpath(dirname, "model.flp")
  name = split(basename(dirname), "_")[end]
  open(filename, "w") do io
    println(io, "'", name, "'")        # 1) title                
    (length(tx) > 1) && (src = "X")  # 2) line source
    (length(tx) == 1) && (src = "R")  # 2) point source
    # 2)  Adiabatic mode thoery TODO: add coupled model
    #TODO: check why coherent option cause error
    print(io, "'", src, "A", "'\n") #2) Coherent / incoherent mode addition
    @printf(io, "%i\n", nmodes) #3) Number of modes
    println(io, "1") #4) Number of profile  #TODO: hardcoded, check for general setting
    println(io, "0.0") #4) Profile ranges (km)  #TODO: hardcoded, check for general setting
    if length(rx) == 1
      print_distance(io, [location(rx[1])[1]./ 1000.0] )  # 6)receiver ranges(km)
      print_distance(io, [-location(tx1)[3] for tx1 ∈ tx])  # Source depth (m) #TODO: check multiple sources how
      print_distance(io, [-location(rx1)[3] for rx1 ∈ rx] )  # 6)receiver depth(m)
      nr = 1
    elseif rx isa AcousticReceiverGrid2D
      d = reverse(-rx.zrange)
      if first(d) > last(d)
       d = reverse(d)
       zrev = true
      end
      r = rx.xrange ./ 1000.0
      if first(r) > last(r)
       r = reverse(r)
       xrev = true
      end
      print_distance(io, r)  # Receiver ranges (km)
      print_distance(io, [-location(tx1)[3] for tx1 ∈ tx])  # Source depth (m)
      print_distance(io, d)  # 6)receiver depth(m)
      nr = length(d)
    end
    print_distance(io, zeros(nr))# Number of receiver range displacement (same as NRD)
  end
end



function readshd_kraken(filename; xrev=false, zrev=false)
  open(filename, "r") do io
    r = read(io, Int32)#read(io, UInt32)
    seek(io, 4r)
    b = Array{UInt8}(undef, 10)
    read!(io, b)
    seek(io, 8r)
    nfreq = read(io, Int32)
    nfreq == 1 || error("Bad shd file format: incorrect nfreq")
    nθ = read(io, Int32)
    nθ == 1 || error("Bad shd file format: incorrect nθ")
    nsx = read(io, Int32)
    nsy = read(io, Int32)
    nsd = read(io, Int32)
    nsd == 1 || error("Bad shd file format: incorrect nsd")
    nrd = read(io, Int32)
    nrr = read(io, Int32)
    pressure = Array{ComplexF32}(undef, nrr, nrd)
    for ird ∈ 0:nrd-1
      recnum = 10 + ird
      seek(io, recnum * 4r)
      temp = Array{Float32}(undef, nrr * 2)
      read!(io, temp)
      pressure[:,ird+1] .= complex.(temp[1 : 2 : nrr * 2], temp[2 : 2:  nrr * 2])
    end
    xrev && (pressure = reverse(pressure; dims=1))
    zrev || (pressure = reverse(pressure; dims=2))
    pressure
  end
end