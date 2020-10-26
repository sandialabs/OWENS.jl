module OWENS

import Statistics
import DelimitedFiles
import LinearAlgebra
import FLOWMath
# import VAWTAero
# using MATLAB #if this is used, must run from src location

export Unsteady #, UnsteadyROM

# export Modal, Flutter
#
# export Steady
#
# export GyricFEA

include("GyricFEA.jl")
include("Modal.jl")
include("Steady.jl")
include("Unsteady.jl")
include("UnsteadyROM.jl")
include("meshing_utilities.jl")
include("structural_utilities.jl")
include("file_reading.jl")
include("structs.jl")




end
