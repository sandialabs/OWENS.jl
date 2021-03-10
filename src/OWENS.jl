module OWENS

import Statistics
import DelimitedFiles
import LinearAlgebra
import FLOWMath
import VAWTAero
import OptimizationParameters
import HDF5
import SparseArrays
import ArnoldiMethod
import Printf
import PyPlot
import PreComp
import Composites
import MAT #for saving .mat H5 files
using pyfloater
pyFloater = pyfloater.pyFloater.pyFloater #simplify the call
using PyCall
wave = pyimport("mhkit.wave")
# ------------ GLOBAL VARIABLES ------------------------------------------------
const module_path = splitdir(@__FILE__)[1]          # Path to this module


# mat_path = string(module_path,"/Matlab_Cxx")
# using MATLAB #if this is used, must run from src location
# mat"addpath($mat_path)"

export Unsteady #, UnsteadyROM
export owens
export Modal#, Flutter
#
# export Steady
#
# export GyricFEA

# "include" is basically the same as copying and pasting the code here
include("GyricFEA.jl")
include("Modal.jl")
include("Steady.jl")
include("Unsteady.jl")
include("UnsteadyROM.jl")
include("meshing_utilities.jl")
include("aero_utilities.jl")
include("structural_utilities.jl")
include("file_io.jl")
include("structs.jl")
include("depreciated.jl")
include("visualization.jl")

end
