import OWENS

runpath = path = @__DIR__
run_full_example = get(ENV, "OWENS_RUN_DOC_EXAMPLES", "false") == "true"

modelopt = OWENS.ModelingOptions("$(path)/OWENS_Opt.yml")
designparams = OWENS.Design_Data("$path/WINDIO_example.yaml")

if run_full_example
    OWENS.runOWENSWINDIO(modelopt,designparams,runpath)
end

nothing

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
