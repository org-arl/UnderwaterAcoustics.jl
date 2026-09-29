using TestItems

@testitem "uacr names" begin
  dir = mktempdir()
  @test_throws "Unknown UACR file \"red\"" uacr_download("red"; dir)
  @test isempty(readdir(dir))
end

@testitem "uacr cache" begin
  # a file of the right size is used as is, without going to the network
  dir = mktempdir()
  path = joinpath(dir, "red_noise.mat")
  write(path, zeros(UInt8, 5176))
  @test uacr_download("red_noise"; dir) == path
  @test uacr_download("red_noise.mat"; dir) == path
  @test read(path) == zeros(UInt8, 5176)
  # a file of any other size may be the user's own, so it is not replaced
  write(path, zeros(UInt8, 10))
  @test_throws "delete it" uacr_download("red_noise"; dir)
  @test filesize(path) == 10
  withenv("UWA_CHANNELS_CACHE" => dir) do
    @test UnderwaterAcoustics._uacr_dir() == dir
  end
  withenv("UWA_CHANNELS_CACHE" => nothing) do
    @test UnderwaterAcoustics._uacr_dir() == joinpath(homedir(), ".cache", "uwa-channels")
  end
end

@testitem "uacr fetch" begin
  import UnderwaterAcoustics: _uacr_fetch, Downloads
  fileurl(path) = "file:///" * replace(lstrip(replace(path, '\\' => '/'), '/'), ' ' => "%20")
  data = UInt8.((1:1000) .% 256)
  src = joinpath(mktempdir(), "src.mat")
  write(src, data)
  dir = mktempdir()
  path = joinpath(dir, "sub", "x.mat")
  @test (@test_logs (:info, r"Downloading x\.mat \(1 kB\)") _uacr_fetch(fileurl(src), path, 1000)) == path
  @test read(path) == data
  # neither a truncated nor a failed download leaves a file behind
  y = joinpath(dir, "y.mat")
  @test_logs (:info, r"Downloading y\.mat") @test_throws "Incomplete download of y.mat (1000 of 1001 bytes)" _uacr_fetch(fileurl(src), y, 1001)
  @test_logs (:info, r"Downloading y\.mat") @test_throws Downloads.RequestError _uacr_fetch(fileurl(src * ".missing"), y, 1000)
  @test readdir(dir) == ["sub"]
  @test readdir(joinpath(dir, "sub")) == ["x.mat"]
end

@testitem "uacr progress" begin
  io = IOBuffer()
  progress = UnderwaterAcoustics._download_progress(io)
  progress(0, 0)  # before the size of the file is known
  foreach(now -> progress(1000, now), (0, 5, 10, 999, 1000))
  @test String(take!(io)) == "\r  0%\r  1%\r  99%\r  100%"
end

@testitem "uacr download" begin
  using MAT: matread
  dir = mktempdir()
  path = try
    uacr_download("red_noise"; dir)
  catch e
    # no HTTP response means no network, which says nothing about the package
    e isa UnderwaterAcoustics.Downloads.RequestError && e.response.status == 0 || rethrow()
    @warn "Skipping UACR download test, since Zenodo is unreachable" exception=e
    nothing
  end
  if path !== nothing
    @test path == joinpath(dir, "red_noise.mat")
    @test issubset(["version", "Fs", "alpha", "beta"], keys(matread(path)))
  end
end
