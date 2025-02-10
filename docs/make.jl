using Documenter
using UnderwaterAcoustics

makedocs(
  sitename = "UnderwaterAcoustics.jl",
  format = Documenter.HTML(prettyurls = false),
  modules = [UnderwaterAcoustics],
  warnonly = [:linkcheck, :missing_docs],
  pages = Any[
    "Home" => "index.md",
    "Manual" => Any[
      "basic.md",
      "quickstart.md",
      "replay.md",
      "api.md",
      "ext.md",
      "utils.md",
      "porting.md",
      "reference.md"
    ]
  ]
)

deploydocs(
  repo = "github.com/org-arl/UnderwaterAcoustics.jl.git",
  branch = "gh-pages",
  devbranch = "master",
  devurl = "dev",
  versions = ["stable" => "v^", "v#.#", "dev" => "dev"]
)
