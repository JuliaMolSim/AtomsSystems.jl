using Documenter
using AtomsSystems
using AtomsSystems.AtomsTrajectories
using AtomsBase
using Unitful

makedocs(
    sitename = "AtomsSystems",
    format = Documenter.HTML(),
    modules = [AtomsSystems, AtomsSystems.AtomsTrajectories ],
    pages = [
        "Home" => "index.md",
        "Atoms" => "atoms.md",
        "Systems" => "systems.md",
        "Utilities" => "utilities.md",
        "Trajectories" => "trajectories.md",
        "Chemfiles Extension" => "chemfiles.md",
        "Index" => "end_index.md",
    ],
    checkdocs = :public,
)


deploydocs(
    repo = "JuliaMolSim/AtomsSystems.jl",
)