# Relocatable OWENS Studio RM2 fixture driver.
using OWENS

rm2_input_root = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "RM2"))
modeling_options_file = joinpath(rm2_input_root, "modeling_options_OWENS_RM2.yml")
windio_file = joinpath(rm2_input_root, "WINDIO_RM2.yaml")
run_path = @__DIR__

OWENS.runOWENSWINDIO(modeling_options_file, windio_file, run_path)
