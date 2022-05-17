module OWENS

# Custom unregistered (julia) packages
try
    import VAWTAero
catch
    @info "Using one way CACTUS Coupling"
    #TODO: propogate logic
end

include("C:/code/GyricFEA.jl/src/GyricFEA.jl")
include("C:/code/VAWTHydro.jl/src/VAWTHydro.jl")  # import VAWTHydro

# Github packages
import Statistics
import DelimitedFiles
import SparseArrays
import LinearAlgebra
import FLOWMath
import HDF5

# using pyfloater
# pyFloater = pyfloater.pyFloater.pyFloater #simplify the call
# wave = pyimport("mhkit.wave")
# ------------ GLOBAL VARIABLES ------------------------------------------------
const module_path = splitdir(@__FILE__)[1]          # Path to this module

Modal = GyricFEA.modal
export Unsteady #, UnsteadyROM
export owens
export Modal#, Flutter
# export Steady

include("Steady.jl")
include("Unsteady.jl")
include("UnsteadyROM.jl")
include("utilities.jl")
include("structs.jl")

end
