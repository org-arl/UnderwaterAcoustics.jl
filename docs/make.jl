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
      "uw_basic.md",
      "pm_basic.md",
      "pm_envref.md",
      "pm_api.md"
    ],
    "Propagation models" => Any[
      "pm_pekeris.md",
      "pm_rays.md",
      "pm_bellhop.md"
    ],
    "Tutorials" => Any[
      "tut_turing.md",
      "tut_autodiff.md"
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
