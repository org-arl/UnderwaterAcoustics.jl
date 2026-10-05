import Downloads

const _UACR_URL = "https://zenodo.org/records/21287414/files/"

_isurl(filename) = occursin(r"^[a-z][a-z0-9+.-]*://"i, filename)

function _url(filename)
  startswith(filename, "uacr://") || return filename
  name = chopprefix(filename, "uacr://")
  _UACR_URL * (endswith(name, ".mat") ? name : name * ".mat")
end

function _filepath(filename)
  url = _url(filename)
  _isurl(url) || return filename
  parts = split(replace(url, r"^[^:]*://" => "", r"[?#].*" => ""), '/'; keepempty=false)
  joinpath(homedir(), ".cache", "uwa-channels", (replace(p, r"^\.\.?$" => "_", r"[^A-Za-z0-9._-]" => "_") for p ∈ parts)...)
end

function _localfile(filename)
  url = _url(filename)
  path = _filepath(filename)
  if _isurl(url) && !isfile(path)
    mkpath(dirname(path))
    @info "Downloading $url to $path"
    tmp = tempname(dirname(path))
    progress = stderr isa Base.TTY ? _download_progress(stderr) : nothing
    try
      Downloads.download(url, tmp; progress)
      mv(tmp, path; force=true)
    finally
      progress === nothing || print(stderr, "\r\e[K")
      rm(tmp; force=true)
    end
  end
  path
end

function _download_progress(io)
  shown = Ref(-1)
  function (total, now)
    total > 0 || return
    pct = floor(Int, 100 * now / total)
    pct > shown[] || return
    shown[] = pct
    print(io, "\r  $pct%")
  end
end
