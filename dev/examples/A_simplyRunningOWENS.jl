import OWENS

runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]
# runpath = path = splitdir(@__FILE__)[1]

modelopt = OWENS.ModelingOptions("$(path)/OWENS_Opt.yml")
designparams = OWENS.Design_Data("$path/WINDIO_example.yaml")

OWENS.runOWENSWINDIO(modelopt,designparams,runpath)

nothing

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
