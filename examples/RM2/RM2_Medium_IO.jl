import OWENS
import YAML
import OrderedCollections

runpath = splitdir(@__FILE__)[1]

OWENS_Options = "$runpath/modeling_options_OWENS_RM2.yml"

WINDIO_filename = "$runpath/WINDIO_RM2.yaml"
windio = YAML.load_file(WINDIO_filename; dicttype=OrderedCollections.OrderedDict{Symbol,Any})

OWENS.runOWENSWINDIO(WINDIO_filename,OWENS_Options,runpath)
