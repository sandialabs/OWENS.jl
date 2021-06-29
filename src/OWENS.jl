module OWENS

# Custom unregistered (julia) packages
import VAWTAero
import GyricFEA
using pyfloater
import PreComp
import Composites
import OptimizationParameters

# Github packages
import Statistics
import DelimitedFiles
import LinearAlgebra
import FLOWMath
import HDF5
import PyPlot
import MAT #for saving .mat H5 files
pyFloater = pyfloater.pyFloater.pyFloater #simplify the call
using PyCall
wave = pyimport("mhkit.wave")
# ------------ GLOBAL VARIABLES ------------------------------------------------
const module_path = splitdir(@__FILE__)[1]          # Path to this module

# mat_path = string(module_path,"/Matlab_Cxx")
# using MATLAB #if this is used, must run from src location
# mat"addpath($mat_path)"
Modal = GyricFEA.modal
export Unsteady #, UnsteadyROM
export owens
export Modal#, Flutter
# export Steady

include("Steady.jl")
include("Unsteady.jl")
include("UnsteadyROM.jl")
include("meshing_utilities.jl")
include("aero_utilities.jl")
include("file_io.jl")
include("structs.jl")
include("depreciated.jl")
include("visualization.jl")

end
