# generate examples
import Literate

is_ci = haskey(ENV, "GITHUB_ACTIONS")

EXAMPLEDIR = joinpath(@__DIR__, "src", "literate")
GENERATEDDIR = joinpath(@__DIR__, "src", "examples")
mkpath(GENERATEDDIR)

# Run Literate on all examples
for myfile in readdir(EXAMPLEDIR)
    if endswith(myfile, ".jl")
        input = abspath(joinpath(EXAMPLEDIR, myfile))
        script = Literate.script(input, GENERATEDDIR)
        code = strip(read(script, String))

        # remove "hidden" lines which are not shown in the markdown
        line_ending_symbol = occursin(code, "\r\n") ? "\r\n" : "\n"
        code_clean = join(filter(x->!endswith(x,"#hide"),split(code, r"\n|\r\n")), line_ending_symbol)
        code_clean = replace(code_clean, r"^# This file was generated .*$"m => "")
        code_clean = strip(code_clean)

        mdpost(str) = replace(str, "@__CODE__" => code_clean)

        Literate.markdown(input, GENERATEDDIR, postprocess = mdpost)
        Literate.notebook(input, GENERATEDDIR, execute = is_ci) # Don't execute locally
    elseif myfile != "vtk"
        @warn "ignoring literate conversion of $myfile, but copying to $GENERATEDDIR"
        cp("$EXAMPLEDIR/$myfile","$GENERATEDDIR/$myfile";force=true)
    end
end

# remove any .vtu files in the generated dir (should not be deployed)
cd(GENERATEDDIR) do
    foreach(file -> endswith(file, ".vtu") && rm(file), readdir())
    foreach(file -> endswith(file, ".pvd") && rm(file), readdir())
end