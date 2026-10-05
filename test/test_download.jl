using TestItems

@testsnippet DownloadSetup begin
  import UnderwaterAcoustics: _url, _filepath, _localfile, Downloads
  fileurl(path) = "file:///" * replace(lstrip(replace(path, '\\' => '/'), '/'), ' ' => "%20")
  withhome(f) = (home = mktempdir(); withenv(f, "HOME" => home, "USERPROFILE" => home))
end

@testitem "download-paths" setup=[DownloadSetup] begin
  @test _url("uacr://red_1") == "https://zenodo.org/records/21287414/files/red_1.mat"
  @test _url("uacr://red_1.mat") == _url("uacr://red_1")
  @test _filepath("red_1.mat") == "red_1.mat"
  cache = joinpath(homedir(), ".cache", "uwa-channels")
  @test _filepath("uacr://red_1") == joinpath(cache, "zenodo.org", "records", "21287414", "files", "red_1.mat")
  @test _filepath("https://example.com:8080/a/b.mat?x=1") == joinpath(cache, "example.com_8080", "a", "b.mat")
  @test _filepath("https://example.com/../../b.mat") == joinpath(cache, "example.com", "_", "_", "b.mat")
end

@testitem "download-cache" setup=[DownloadSetup] begin
  data = UInt8.((1:1000) .% 256)
  src = joinpath(mktempdir(), "src.mat")
  write(src, data)
  withhome() do
    path = @test_logs (:info, r"Downloading file://") _localfile(fileurl(src))
    @test read(path) == data
    rm(src)
    @test (@test_logs _localfile(fileurl(src))) == path
    @test_logs (:info, r"Downloading") @test_throws Downloads.RequestError _localfile(fileurl(src * ".x"))
    @test readdir(dirname(path)) == ["src.mat"]
  end
end

@testitem "download-replay-channel" setup=[DownloadSetup] begin
  using MAT: matwrite
  src = joinpath(mktempdir(), "channel.mat")
  matwrite(src, Dict(
    "version" => 1.0,
    "h_hat" => ones(ComplexF64, 4, 1, 10),
    "params" => Dict("fs_delay" => 1000.0, "fs_time" => 100.0, "fc" => 300.0),
  ))
  withhome() do
    ch = @test_logs (:info, r"Downloading") BasebandReplayChannel(fileurl(src))
    @test ch.h == BasebandReplayChannel(src).h
  end
end

@testitem "download-progress" begin
  io = IOBuffer()
  progress = UnderwaterAcoustics._download_progress(io)
  progress(0, 0)
  foreach(now -> progress(1000, now), (0, 5, 10, 999, 1000))
  @test String(take!(io)) == "\r  0%\r  1%\r  99%\r  100%"
end

@testitem "download-zenodo" setup=[DownloadSetup] begin
  using MAT: matread
  withhome() do
    path = try
      _localfile("uacr://red_noise")
    catch e
      # no network
      e isa Downloads.RequestError && e.response.status == 0 || rethrow()
      @warn "Zenodo unreachable; skipping test" exception=e
      nothing
    end
    path === nothing || @test issubset(["version", "Fs", "alpha", "beta"], keys(matread(path)))
  end
end
