import Downloads

export uacr_download

"""
    uacr_download(name; dir)

Download a channel or noise file from the underwater acoustic channel
repository (UACR), and return its local path. `name` is the name of the file,
with or without the `.mat` extension, e.g. `"red_1"` for a channel or
`"red_noise"` for a noise model.

The available files are:
- channels `blue_1` … `blue_20`, `red_1` … `red_4`, `yellow_1` … `yellow_6`,
  `purple_1` … `purple_15`, `green`, `black`, `pink_1` … `pink_3` and `brown`
- noise models `blue_noise`, `red_noise`, `yellow_noise_1` … `yellow_noise_4`,
  `purple_noise_1` … `purple_noise_5`, `green_noise`, `black_noise` and
  `pink_noise_1` … `pink_noise_3`

See https://uwa-channels.github.io/ for a description of each channel, and the
noise model that goes with it. Channel files are large (160 MB to 1.5 GB), while
noise files are under 1 MB.

Files are downloaded from version 1.0 of the repository on Zenodo
(https://doi.org/10.5281/zenodo.21287414) to directory `dir`. A file that is
already in `dir` is not downloaded again. `dir` defaults to the directory named
by the `UWA_CHANNELS_CACHE` environment variable if it is set, and to
`~/.cache/uwa-channels` otherwise.

# Examples
```julia-repl
julia> using UnderwaterAcoustics, MAT

julia> ch = BasebandReplayChannel(uacr_download("red_1"))
[ Info: Downloading red_1.mat (192 MB) from Zenodo to /home/user/.cache/uwa-channels
BasebandReplayChannel(3 × 47.8 s, 25000.0 Hz, 19200.0 Sa/s)
```
"""
function uacr_download(name::AbstractString; dir=_uacr_dir())
  name = chopsuffix(name, ".mat")
  i = findfirst(f -> f.first == name, _UACR_FILES)
  i === nothing && error("Unknown UACR file \"$name\"; see ?uacr_download for the available files")
  url = "https://zenodo.org/api/records/$_UACR_RECORD/files/$name.mat/content"
  _uacr_fetch(url, joinpath(abspath(expanduser(dir)), "$name.mat"), _UACR_FILES[i].second)
end

# the cache used by the tests of the UACR Python toolbox, so that the two share
# downloaded files
function _uacr_dir()
  dir = get(ENV, "UWA_CHANNELS_CACHE", "")
  isempty(dir) ? joinpath(homedir(), ".cache", "uwa-channels") : dir
end

# downloads to a uniquely named file next to path and then renames it, so that a
# failed download never leaves a partial file at path
function _uacr_fetch(url, path, nbytes)
  if isfile(path)
    filesize(path) == nbytes && return path
    error("$path has $(filesize(path)) bytes, but $(basename(path)) in the UACR has $nbytes; delete it to download a fresh copy")
  end
  dir = dirname(path)
  mkpath(dir)
  @info "Downloading $(basename(path)) ($(_bytes2str(nbytes))) from Zenodo to $dir"
  tmp = tempname(dir)
  # only on a terminal, since in logs and notebooks each update would be a new line
  progress = stderr isa Base.TTY ? _download_progress(stderr) : nothing
  try
    Downloads.download(url, tmp; progress)
    filesize(tmp) == nbytes || error("Incomplete download of $(basename(path)) ($(filesize(tmp)) of $nbytes bytes); please try again")
    mv(tmp, path; force=true)
  finally
    progress === nothing || print(stderr, "\r\e[K")
    rm(tmp; force=true)
  end
  path
end

function _download_progress(io)
  shown = Ref(-1)
  function (total, now)
    total > 0 || return
    pct = floor(Int, 100 * now / total)  # 100 * now overflows a 32-bit Int
    pct > shown[] || return
    shown[] = pct
    print(io, "\r  $pct%")
  end
end

function _bytes2str(n)
  n < 1e6 && return "$(round(Int, n / 1e3)) kB"
  n < 1e9 && return "$(round(Int, n / 1e6)) MB"
  "$(round(n / 1e9; digits=1)) GB"
end

# a version-specific record, whose files never change, so a downloaded file
# stays valid for as long as the package points to this record
const _UACR_RECORD = "21287414"

# file sizes in bytes, as listed at https://zenodo.org/api/records/21287414; the
# size is checked rather than Zenodo's MD5 checksum, since Julia has no MD5 stdlib
const _UACR_FILES = [
  "blue_1"         => 206132594,
  "blue_2"         => 206132594,
  "blue_3"         => 206132594,
  "blue_4"         => 206132594,
  "blue_5"         => 206132594,
  "blue_6"         => 206132594,
  "blue_7"         => 206132594,
  "blue_8"         => 206132594,
  "blue_9"         => 206132594,
  "blue_10"        => 206132602,
  "blue_11"        => 206132602,
  "blue_12"        => 206132602,
  "blue_13"        => 206132602,
  "blue_14"        => 206132602,
  "blue_15"        => 206132602,
  "blue_16"        => 206132602,
  "blue_17"        => 206132602,
  "blue_18"        => 206132602,
  "blue_19"        => 206132602,
  "blue_20"        => 206132602,
  "red_1"          => 191566442,
  "red_2"          => 191566442,
  "red_3"          => 191566442,
  "red_4"          => 191566442,
  "yellow_1"       => 821046493,
  "yellow_2"       => 821046493,
  "yellow_3"       => 984922035,
  "yellow_4"       => 360679714,
  "yellow_5"       => 834576462,
  "yellow_6"       => 394887977,
  "purple_1"       => 712812521,
  "purple_2"       => 712812521,
  "purple_3"       => 534541341,
  "purple_4"       => 534541341,
  "purple_5"       => 354079789,
  "purple_6"       => 712812521,
  "purple_7"       => 712812521,
  "purple_8"       => 534541341,
  "purple_9"       => 534541341,
  "purple_10"      => 354079789,
  "purple_11"      => 712812521,
  "purple_12"      => 712812521,
  "purple_13"      => 534541341,
  "purple_14"      => 534541341,
  "purple_15"      => 354079789,
  "green"          => 1094029080,
  "black"          => 161710585,
  "pink_1"         => 1496605964,
  "pink_2"         => 1496605964,
  "pink_3"         => 1440393938,
  "brown"          => 660486317,
  "blue_noise"     => 155097,
  "red_noise"      => 5176,
  "yellow_noise_1" => 614833,
  "yellow_noise_2" => 614833,
  "yellow_noise_3" => 313722,
  "yellow_noise_4" => 313722,
  "purple_noise_1" => 580006,
  "purple_noise_2" => 580006,
  "purple_noise_3" => 306031,
  "purple_noise_4" => 306031,
  "purple_noise_5" => 117081,
  "green_noise"    => 7416,
  "black_noise"    => 104793,
  "pink_noise_1"   => 614833,
  "pink_noise_2"   => 614833,
  "pink_noise_3"   => 514576,
]
