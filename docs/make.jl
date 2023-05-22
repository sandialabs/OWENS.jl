using Documenter, Literate, OWENS

# Build documentation
makedocs(;
    modules = [OWENS],
    pages = [
        "Home" => "index.md",
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENS.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
)