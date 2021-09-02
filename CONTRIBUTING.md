# Contributing to UnderwaterAcoustics.jl

Contributions in terms of bug reports, feature requests, ideas/suggestions, documentation improvements, bug fixes, and new feature implementations are most welcome! Even contributions on improving this document are welcome!

If you intend to contribute, please read the guidelines below:

## Getting started

We assume familiarity with git, Github, markdown and Julia. If you need to brush up on these, here are some good starting points:

* [Github](https://docs.github.com/en/get-started/quickstart/set-up-git)
* [Markdown](https://guides.github.com/features/mastering-markdown/)
* [Julia](https://julialang.org/learning/)

While we haven't adopted a formal code of conduct, there is an implicit expectation that all of us shall be professional, respectful of differing opinions and viewpoints, empathetic and kind, and open to giving and gracefully accepting constructive feedback. By contributing, you accept to abide by this code.

Finally, read and understand the [ROADMAP](ROADMAP.md) document to understand what UnderwaterAcoustics.jl is meant to be, and meant not to be, so that all of us work from a common set of expectations. If you have suggestions on the scope, directions or roadmap for the project, we are open to [discussing](https://github.com/org-arl/UnderwaterAcoustics.jl/discussions) those too.

## Bug reports, features requests & discussions

### Bug reports

Bug reports are an important part of making UnderwaterAcoustics.jl more stable and reliable. Having a good bug report will allow others to reproduce it and provide insight into fixing. See [this stackoverflow article](https://stackoverflow.com/help/mcve) for tips on writing a good bug report.

Trying out the bug-producing code on the `master` branch is often a worthwhile exercise to confirm that the bug still exists. You'd also want to search for existing bug reports and pull requests to see if the issue has already been reported and/or fixed before [filing a new issue](https://github.com/org-arl/UnderwaterAcoustics.jl/issues/new/choose).

### Feature requests

If you have ideas for new features that fit within the [scope of the project](ROADMAP.md), feel free to [file a new issue](https://github.com/org-arl/UnderwaterAcoustics.jl/issues/new/choose) with your idea. Providing some background on the motivation behind the feature request, and some use cases of how the feature might be used, is very valuable for others to understand your idea.

### Discussions

If you have things you'd like to talk about that don't naturally fit in as a bug report or feature request, [open a new discussion](https://github.com/org-arl/UnderwaterAcoustics.jl/discussions) for it. Discussions may include Q&A, project scope/roadmap discussions, ideas that aren't concrete enough to be feature requests yet, cool things you've done with UnderwaterAcoustics.jl that you'd want to show others, etc.

## Bug fixes, new features & documentation enhancements

In order to contribute to the code and/or documentation for UnderwaterAcoustics.jl, you'll need to be familiar with git and have a Github account. We have some useful links in the [Getting Started](#getting-started) section above, if you're new to git or Github.

The general process to contribute code or documentation is as follows:
1. Fork the repository, and make a clone of the fork on your own machine.
2. Create a new branch based on `master` and checkout the new branch. The name of the branch should be descriptive but short, with only lowercase letters and hyphens in it. Each branch should only contain changes for a specific issue, or a related set of isses.
3. Make changes to the codebase (`src` folder) and/or documentation (`docs/src` folder) locally. You may find it useful to have [VSCode](https://code.visualstudio.com) with the [Julia extension](https://github.com/julia-vscode/julia-vscode) installed to work with the codebase, although this is not strictly necessary.
4. If you're fixing a bug, it may be worth adding test cases (`test` folder) that catch the bug before you fix it, and then running those test cases to ensure that your fix works.
5. If you're adding a new feature, you may want to develop test cases for the feature before or along with your development. Include your test cases in the test suite (in the `test` folder). You'd also want to add documentation in the form of [docstrings](https://docs.julialang.org/en/v1/manual/documentation/) within your code and/or markdown documentation (`docs/src` folder). We use [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) to generate the online documentation automatically.
6. Test the changes locally, making sure that your code works as advertised, and that you haven't broken anything else. In order to do that:
    - Run the regression test suite locally:
      ```sh
      $ julia --project                 # run Julia
      julia>                            # press ] for package mode
      (UnderwaterAcoustics) pkg> test   # start the test
      ```
      Make sure that all tests pass.
    - Build the documentation locally and check it:
      ```sh
      $ cd docs
      $ julia --project make.jl
      ```
      After the build process, the documentation is in the `docs/build` folder. Open it with a browser and check it.
7. Commit your changes (you may want to do this often during your development and testing phases) and push them to your forked repository. Ensure you use [good commit messages](#commit-messages).
8. Raise a [pull request (PR)](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) from your branch in your forked repository to the `master` branch in the main repository. Provide details of your changes, and links to one or more original issues that may have led to this PR. You may raise PRs before you are ready with the final changes, but please mark them as _draft_ (Github allows you to create a _draft PR_) until it's ready for others to review.
9. Once your PR is ready for review, one of the maintainers of the repository will assign a reviewer. The job of the reviewer is to ensure code and documentation quality before the PR is merged. The reviewer will provide constructive feedback to help address any concerns with the PR.
10. Once all concerns are addressed and marked _resolved_, the PR will be ready for merging. One of the maintainers will merge the PR, and the changes will become part of the next release of UnderwaterAcoustics.jl.

## Commit messages

We have guidelines on how our git commit messages should be formatted. This leads to more readable messages that are easy to follow when looking through the project history, and also allow automated generation of change logs.

Each commit message consists of a **header**, and optionally, a **body**. The header has a special format that includes a **type**, a **scope** and a **summary**:
```
type(scope): summary

body
```
The summary should be kept to about 50 characters, and each line in the body should not exceed 80 characters. If the commit fixes or addresses any issue, it should be referenced in the summary or the body. Subject lines should be succint and use imperative voice and present tense. For the subject line, don't capitalize the first letter and don't add a dot (.) at the end.

We should strive to maintain backward compatibility as much as possible. However, if a commit breaks backward compatibility, it MUST be flagged. To do so, start the body with `"BREAKING CHANGE:"`.

Allowable **type**s include:
- **feat**: new feature
- **fix**: bug fix
- **docs**: documentation-only changes
- **test**: new test cases or fixes to existing test cases
- **perf**: code changes that improve performance
- **refactor**: code change that improves codebase without adding features or fixing bugs
- **style**: changes that do not affect the code (e.g. whitespace, indentation, etc)
- **chore**: changes to build scripts, CI configuration, gitignore, etc.
- **revert**: commit reverts a previous commit (provide SHA of previous commit)

Allowable **scope** include:
- **uwa**: basic underwater acoustic functionality
- **pm**: core infrastructure for propagation modeling
- **pekeris**: Pekeris ray propagation model
- **bellhop**: Bellhop ray propagation model
- **raysolv** RaySolver ray propagation model
- **plot**: plotting support functionality

NOTE: As we evolve the UnderwaterAcoustics.jl library, the allowable scope list will evolve. If you find something you're doing doesn't fit in this scope list, please [open an issue](https://github.com/org-arl/UnderwaterAcoustics.jl/issues/new/choose) to propose addition of a new scope.

In some cases, the scope of a commit cuts across multiple scopes. In that case, a comma-separated scope list may be used. In cases, where all scopes are affected, a "`*`" may be used as scope, or the scope may be omitted.

Some examples of good commit messages:
```
fix(raysolv): fix type instability for Float32 (issue #17)
```
```
feat(uwa): add bubble resonance calculator

Compute resonance frequency of a freely oscillating has bubble in water
using equation from based on Medwin & Clay (1998).
```
```
fix(pm): return complex transmission loss from models

BREAKING CHANGE: Propagation models used to return just the magnitude of
transmission loss, and ignore phase. We have updated the API to now return
complex transmission loss instead, and let users take the `abs()`
of it, if they need only the magnitude.
```

## Coding standards

We have not evolved a formal coding standard yet, but good coding practices common in the Julia community are adopted in this project. At minimum, you should be familiar (and comply to the extent reasonable) with [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/), [Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/) and the [Julia documentation guidelines](https://docs.julialang.org/en/v1/manual/documentation/).

### Some guidelines:
- Properly format and indent your code. Use 2 spaces for intentation (DO NOT use tabs).
- Use blank lines and/or comments to separate code sections for readability. But do not use excess whitespace or blank lines, as this limits the amount of visible code on the screen and makes the code harder to navigate.
- Do not leave chunks of commented code in the codebase. Delete unused code, and rely on git history for access to old code when necessary.
- Type names start with a capital letter and use CamelCase (e.g. `RaySolver`).
- Function and variable names are generally all lowercase. Multiple words are concatenated together without a underscore. However, the use of underscore or uppercase is permitted when concatenated words are difficult to read (e.g. `bubbleresonance()`, `surfaceloss`, `transmissionloss_dB()`).
- Use descriptive function and variable names, and avoid abbreviations unless obvious or common (e.g. prefer `count` over `cnt`). Local index variables may use single letter variables (e.g. `i`, `j`, `k`, etc).
- When implementing mathematical algorithms, it is fine to use the notation used in the original papers, including single letter variables, Greek letters, unicode characters, and capital letters for matrices. When doing so, however, it is recommended that a citation to the original paper be made available in the code as a comment, or in the docstring.
- Code should be as self-documenting as possible, and comments that could go out of sync with code should be avoided. However, in cases where additional information would be useful for the reader, comments may be included.
- Comments may be prefixed with `NOTE:`, `FIXME:` or `TODO:`, if the developer wishes to leave some thoughts for oneself or a later developer to be reminded of, understand or address.
- When docstrings are included, they MUST either include the method signature, or the text "`$(SIGNATURES)`" (or equivalent). We use [DocStringExtensions.jl](https://juliadocs.github.io/DocStringExtensions.jl/stable/) to expand "`$(SIGNATURES)`" out into the method signatures, when generating the online documentation.
- Be mindful of performance (memory organization, type stability, etc) when coding, as much of the code in this project gets used in computationally sensitive applications. The [Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/) section provides a good guidelines on this.
- When in doubt, be guided by other parts of the existing codebase. If still unsure, [ask](https://github.com/org-arl/UnderwaterAcoustics.jl/discussions).

## Acknowledgments

In preparing this document, we used the following references:
- [Contributing to xarray](http://xarray.pydata.org/en/stable/contributing.html)
- [Contributing to Github](https://github.com/github/docs/blob/main/CONTRIBUTING.md)
- [Angular commit format reference sheet](https://gist.github.com/brianclements/841ea7bffdb01346392c)
