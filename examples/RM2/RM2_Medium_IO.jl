import OWENS
import YAML
import OrderedCollections
using Test
import HDF5

runpath = splitdir(@__FILE__)[1]

OWENS_Options = "$runpath/modeling_options_OWENS_RM2.yml"

WINDIO_filename = "$runpath/WINDIO_RM2.yaml"

OWENS.runOWENSWINDIO(OWENS_Options,WINDIO_filename,runpath)