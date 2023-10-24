using Documenter, Literate, OWENS

# Build documentation
makedocs(;
    modules = [OWENS],
    pages = [
        "Home" => "index.md",
        "Installation" => "setup.md",
        "Developer Guide" => "OWENS_Dev_Guide.md",
        "Frames of Reference" => "FramesOfReference.md",
        "API Reference" => joinpath("reference", "reference.md"),
        "Legacy User Guide" => "legacyUserGuide.md",
        "Legacy VAWTGen Guide" => "VAWTGenUserGuide.md",
    ],
    sitename = "OWENS.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
    remotes = nothing
)