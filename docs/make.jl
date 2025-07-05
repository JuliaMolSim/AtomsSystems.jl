using Documenter
using AtomsSystems
using AtomsBase
using Unitful

makedocs(
    sitename = "AtomsSystems",
    format = Documenter.HTML(),
    modules = [AtomsSystems],
    pages = [
        "Home" => "main_page.md",
        "Atoms" => "atoms.md",
        "Systems" => "systems.md",
        "Utilities" => "utilities.md",
        "Index" => "index.md",
    ],
    checkdocs = :public,
)


deploydocs(
    repo = "JuliaMolSim/AtomsSystems.jl",
)