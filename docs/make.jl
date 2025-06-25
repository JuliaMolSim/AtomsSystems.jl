using Documenter
using AtomsSystems
using AtomsBase
using Unitful

makedocs(
    sitename = "AtomsSystems",
    format = Documenter.HTML(),
    modules = [AtomsSystems],
    pages = [
        "Home" => "index.md",
    ],
    checkdocs = :public,
)


deploydocs(
    repo = "JuliaMolSim/AtomsSystems.jl",
)