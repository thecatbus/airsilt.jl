# docs/make.jl

using AIRSilt, Documenter, Oscar

makedocs(sitename="AIRSilt.jl",
    modules=[AIRSilt],
    repo="https://github.com/thecatbus/airsilt.jl")

deploydocs(
    repo="://github.com",
    devbranch="main"
)
