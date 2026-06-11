using Test

const OWENS_EXAMPLE_ROOT = normpath(joinpath(@__DIR__, "..", "examples"))

function _owens_example_julia_files()
    files = String[]
    for (root, _, names) in walkdir(OWENS_EXAMPLE_ROOT)
        for name in names
            endswith(name, ".jl") || continue
            push!(files, joinpath(root, name))
        end
    end
    return sort(files)
end

function _heavyweight_example_smoke_files()
    files = filter(
        endswith(".jl"),
        readdir(joinpath(OWENS_EXAMPLE_ROOT, "literate"); join = true),
    )
    push!(files, joinpath(OWENS_EXAMPLE_ROOT, "Optimization", "optimizationExample.jl"))
    return sort(normpath.(files))
end

function _windio_reversed_call_locations(text::AbstractString)
    pattern = r"runOWENSWINDIO\(\s*(?:windio|WINDIO_filename|designparams)\s*,\s*(?:OWENS_Options|modelopt)\s*,"
    locations = Tuple{Int,String}[]
    lines = split(text, '\n'; keepempty = true)
    for (line_number, line) in enumerate(lines)
        occursin(pattern, line) || continue
        push!(locations, (line_number, strip(line)))
    end
    return locations
end

function _parse_example_failure(path::AbstractString)
    try
        Meta.parseall(read(path, String); filename = path)
        return nothing
    catch err
        return sprint(showerror, err)
    end
end

@testset "runOWENSWINDIO example argument order" begin
    offenders = String[]
    for file in _owens_example_julia_files()
        locations = _windio_reversed_call_locations(read(file, String))
        isempty(locations) && continue
        relative = relpath(file, OWENS_EXAMPLE_ROOT)
        for (line_number, line) in locations
            push!(offenders, "$relative:$line_number: $line")
        end
    end

    @test offenders == String[]
end

@testset "heavyweight examples parse" begin
    failures = String[]
    for file in _heavyweight_example_smoke_files()
        failure = _parse_example_failure(file)
        isnothing(failure) && continue
        push!(failures, "$(relpath(file, OWENS_EXAMPLE_ROOT)): $failure")
    end

    @test failures == String[]
end
