module OWENS

# Custom unregistered (julia) packages
import GyricFEA
import OpenFASTWrappers

# Github packages
import Statistics
import DelimitedFiles
import LinearAlgebra
import SparseArrays
using StaticArrays
import FLOWMath
import HDF5
import GXBeam
import QuadGK
import Composites
import YAML

import PreComp
import Dierckx
using WriteVTK

# using pyfloater
# pyFloater = pyfloater.pyFloater.pyFloater #simplify the call
# wave = pyimport("mhkit.wave")
# ------------ GLOBAL VARIABLES ------------------------------------------------
const module_path = splitdir(@__FILE__)[1]          # Path to this module

Modal = GyricFEA.modal
FEAModel = GyricFEA.FEAModel
Mesh = GyricFEA.Mesh
El = GyricFEA.El
initialElementCalculations = GyricFEA.initialElementCalculations
export Unsteady #, UnsteadyROM
export owens
export Modal#, Flutter
# export Steady
function __init__()

println("                                           NOTICE:
For five (5) years from 10/15/2020, the United States Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government. There is provision for the possible extension of the term of this license. Subsequent to that period or any extension granted, the United States Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so. The specific term of the license can be identified by inquiry made to National Technology and Engineering Solutions of Sandia, LLC or DOE.
NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
Any licensee of \"Offshore Wind Energy Simulator (OWENS)\" has the obligation and responsibility to abide by the applicable export control laws, regulations, and general prohibitions relating to the export of technical data. Failure to obtain an export control license or other authority from the Government may result in criminal liability under U.S. laws.
                                             (End of Notice)")
end

include("Steady.jl")
include("Unsteady.jl")
include("Unsteady_Land.jl")
include("utilities.jl")
include("Unsteady_utilities.jl")
include("structs.jl")
include("./PostProcessing.jl")
include("./visualization.jl")
include("./AeroMapping.jl")
include("./fileio.jl")
include("./meshing_utilities.jl")
include("./gxbeam_conversion.jl")
include("./SetupTurbine.jl")
include("./topRunDLC.jl")

end
