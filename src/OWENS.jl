module OWENS

using PyCall
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

# ------------ GLOBAL VARIABLES ------------------------------------------------
const module_path = splitdir(@__FILE__)[1]          # Path to this module


# mat_path = string(module_path,"/Matlab_Cxx")
# using MATLAB #if this is used, must run from src location
# mat"addpath($mat_path)"

export Unsteady #, UnsteadyROM
export owens #TODO: do this right

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
include("aero_utilities.jl")
include("structural_utilities.jl")
include("file_io.jl")
include("structs.jl")
include("depreciated.jl")
include("visualization.jl")


# ------------ LOAD airfoilprep.py ---------------------------------------------
path_hydro = module_path                    # Path to tlp_platform.py
hydro = PyNULL()                                    # tlp_platform module

function __init__()
    imp = pyimport("imp")
    (file, filename, data) = imp.find_module("tlp_platform", [path_hydro])
    copy!(hydro, imp.load_module("tlp_platform", file, filename, data))
end

end
