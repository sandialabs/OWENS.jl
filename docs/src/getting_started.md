# Getting Started

OWENS Studio currently provides a dependency-light Julia command surface and
static HTML workbench artifacts. It is not yet a packaged desktop application or
launched HTTP server.

## Studio Diagnostics

From the OWENS.jl repository, run:

```julia
julia --project=. app/OWENS_APP/src/OWENS_APP.jl studio-doctor .
```

The command returns an `owens-studio-doctor/v1` YAML payload with Julia version,
route/catalog availability, template/example availability, dependency visibility,
and whether the requested output directory can be written.

Inspect the current Studio quality gates:

```julia
julia --project=. app/OWENS_APP/src/OWENS_APP.jl project-quality-gates
```

The `owens-studio-quality-gates/v1` payload lists the evidence required before a
Studio feature should be called done, including service contracts, file
provenance, generated-script equivalence, physical validation, browser tests,
visual regression, accessibility, documentation, and performance coverage.
The same payload is available to the dependency-light route layer at
`/api/quality-gates`.

The same payload is also exposed by the dependency-light route contract as
`/api/doctor` for future local-server wrappers.

## Static Studio Start

Create a static project gallery:

```julia
julia --project=. app/OWENS_APP/src/OWENS_APP.jl studio-home owens_studio_home.html
```

Create an RM2 Studio project fixture:

```julia
julia --project=. app/OWENS_APP/src/OWENS_APP.jl project-template rm2 owens-rm2-studio
```

Inspect the project:

```julia
julia --project=. app/OWENS_APP/src/OWENS_APP.jl project-open owens-rm2-studio/owens_project.yml
```

Inspect editable inputs with an explicit preview/validation byte limit:

```julia
julia --project=. app/OWENS_APP/src/OWENS_APP.jl project-inputs owens-rm2-studio/owens_project.yml true 200000
```

Write a static workbench bundle:

```julia
julia --project=. app/OWENS_APP/src/OWENS_APP.jl project-bundle owens-rm2-studio/owens_project.yml owens-rm2-bundle
```

The generated bundle contains HTML plus YAML payloads that preserve project,
input, session, and run-manifest health. Run-manifest artifact drift is reported
as project `attention` until the stale input/output/generated artifact is
restored or the run is regenerated.
