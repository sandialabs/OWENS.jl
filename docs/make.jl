using Documenter, Literate, OWENS

# Generate examples
include("generate.jl")

# Build documentation
makedocs(;
    modules = [OWENS],
    pages = [
        "Home" => "index.md",
        "Installation" => "setup.md",
        "Examples" => [
            joinpath("examples", "A_simplyRunningOWENS.md"),
            joinpath("examples", "B_detailedInputs.md"),
            joinpath("examples", "C_customizablePreprocessing.md"),
        ],
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

deploydocs(
    repo = "github.com/sandialabs/OWENS.jl.git",
)

# ## Documentation
# Until public hosting of the documentation is set up, a readthedocs style webpage can be built via:

#     cd path2OWENS.jl/OWENS.jl/docs
#     julia --project make.jl

# and then a local server can be started via

#     cd ..
#     julia -e 'using LiveServer; serve(dir="docs/build")'

# then open your favorite browser and open the following (or what is indicated in the terminal output if different)

#     http://localhost:8000/