using Documenter, Literate, OWENS

# Generate examples
include("generate.jl")

# Build documentation
makedocs(;
    modules = [OWENS],
    pages = [
        "Home" => "index.md",
        "Installation" => "indepth_installation.md",
        "Examples" => [
            joinpath("examples", "A_simplyRunningOWENS.md"),
            joinpath("examples", "B_detailedInputs.md"),
            joinpath("examples", "C_customizablePreprocessing.md"),
            joinpath("examples", "D_simulatingFloatingPlatforms.md"),
            joinpath("examples", "E_RM2_Medium.md"),
        ],
        "Developer Guide" => "OWENS_Dev_Guide.md",
        "API Reference" => [
            joinpath("reference", "reference.md"),
            joinpath("reference", "referenceAero.md"),
            joinpath("reference", "referenceFEA.md"),
            joinpath("reference", "referenceOpenFASTWrappers.md"),
            joinpath("reference", "referencePreComp.md"),
        ],
        "Legacy User Guide" => "legacyUserGuide.md",
        "Legacy VAWTGen Guide" => "VAWTGenUserGuide.md",
    ],
    sitename = "OWENS.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
    remotes = nothing,
    format = Documenter.HTML(
        size_threshold = 5 * 1024 * 1024,  # 5 MiB
        size_threshold_warn = 1 * 1024 * 1024  # 1 MiB
    )
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