module OWENS

# Custom unregistered (julia) packages
try
    import VAWTAero
catch
    @info "Using one way CACTUS Coupling"
    #TODO: propogate logic
end

import GyricFEA
import VAWTHydro

# Github packages
import Statistics
import DelimitedFiles
import LinearAlgebra
import SparseArrays
import FLOWMath
import HDF5
import GXBeam

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
include("utilities.jl")
include("structs.jl")

end
