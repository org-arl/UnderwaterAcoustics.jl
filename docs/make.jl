using Documenter
# using Plots

# push!(LOAD_PATH,"../src/")
using UnderwaterAcoustics

makedocs(
  sitename = "UnderwaterAcoustics.jl",
  format = Documenter.HTML(prettyurls = false),
  modules = [UnderwaterAcoustics],
  warnonly = [:linkcheck, :missing_docs],
  pages = Any[
    "Home" => "index.md",
    "Manual" => Any[
      "uw_basic.md",
      "pm_api.md",
      "pm_ext.md",
      "utils.md"
    ],
    # "Propagation models" => Any[
    #   "pm_pekeris.md"
    # ],
    # "Tutorials" => Any[
    #   "tut_turing.md",
    #   "tut_autodiff.md"
    # ]
  ]
)

deploydocs(
  repo = "github.com/org-arl/UnderwaterAcoustics.jl.git",
  branch = "gh-pages",
  devbranch = "master",
  devurl = "dev",
  versions = ["stable" => "v^", "v#.#", "dev" => "dev"]
)
