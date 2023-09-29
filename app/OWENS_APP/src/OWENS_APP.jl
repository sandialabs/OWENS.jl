module OWENS_APP
import OWENS
import VAWTAero
import YAML


# pathtest,_ = splitdir(@__FILE__)

function julia_main()::Cint
    try
        real_main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function real_main()
    # println(path)
    path = ARGS[1]
    println(path)
    # println(pathtest)
    yamlInputfile = "$path/sampleOWENS.yml"#ARGS[1]

    Inp = OWENS.MasterInput(yamlInputfile)
    OWENS.runOWENS(Inp,path)
end

end #module