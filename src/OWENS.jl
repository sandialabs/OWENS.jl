__precompile__(false)
module OWENS

# Custom unregistered (julia) packages
import OWENSFEA
import OWENSOpenFASTWrappers
import OWENSAero

# Github packages
import Statistics
import Statistics:mean
import DelimitedFiles
import LinearAlgebra
import SparseArrays
using StaticArrays
import FLOWMath
import HDF5
import JLD2
import GXBeam
import QuadGK
import Composites
import YAML
import OrderedCollections
import Base.show
import ProgressBars

import OWENSPreComp
import Dierckx
using WriteVTK

# using pyfloater
# pyFloater = pyfloater.pyFloater.pyFloater #simplify the call
# wave = pyimport("mhkit.wave")
# ------------ GLOBAL VARIABLES ------------------------------------------------
const module_path = splitdir(@__FILE__)[1]          # Path to this module

Modal = OWENSFEA.modal
AutoCampbellDiagram = OWENSFEA.autoCampbellDiagram
FEAModel = OWENSFEA.FEAModel
Mesh = OWENSFEA.Mesh
El = OWENSFEA.El
initialElementCalculations = OWENSFEA.initialElementCalculations
export Unsteady #, UnsteadyROM
export owens
export Modal#, Flutter
# export Steady
function __init__()

println("\nThis program is running OWENS.jl, the Offshore/Onshore Wind/Water Energy Simulator for turbine type devices

Copyright 2013-2025
National Technology & Engineering Solutions of Sandia, LLC (NTESS). 
Under the terms of Contract DE-NA0003525 with NTESS, 
the U.S. Government retains certain rights in this software.

Licensed under the LGPL GNU LESSER GENERAL PUBLIC LICENSE V3.0 License;
you may not use this file except in compliance with the License.

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an \"AS IS\" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.\n")
end

# Core functionality - fundamental data structures and configuration
include("core/ModelingOptions.jl")
include("core/structs.jl")
include("core/utilities.jl")
include("core/TurbineModel.jl")
include("core/InputTranslation.jl")
include("core/UnifiedDriver.jl")


# Input/Output functionality
include("io/fileio.jl")
include("io/windio.jl")

# Preprocessing functionality
include("preprocessing/gxbeam_conversion.jl")
include("preprocessing/meshing_utilities.jl")
include("preprocessing/setup_utilities.jl")
include("preprocessing/SetupTurbine.jl")
# include("preprocessing/SetupTurbineHAWT.jl")  # Temporarily disabled

# Runtime functionality
include("runtime/Steady.jl")
include("runtime/Unsteady.jl")
include("runtime/Unsteady_Land.jl")
include("runtime/Unsteady_utilities.jl")
include("runtime/DLCAnalysis.jl")
include("postprocess/AeroMapping.jl")

``
# Postprocessing functionality
include("postprocess/PostProcessing.jl")
include("postprocess/visualization.jl")

end
