using Documenter
using Plots

push!(LOAD_PATH,"../src/")
using UnderwaterAcoustics

makedocs(
  sitename = "UnderwaterAcoustics.jl",
  format = Documenter.HTML(prettyurls = false),
  linkcheck = !("skiplinks" in ARGS),
  pages = Any[
    "Home" => "index.md",
    "Manual" => Any[
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
