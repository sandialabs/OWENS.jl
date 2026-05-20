# OWENS
```@meta
CurrentModule = OWENS
```

## Index

```@index
```

## Types and functions

```@autodocs
Modules = [OWENS]
```

## Output Data Channels

`OWENS.outputData` writes a root-level HDF5 file and annotates each dataset
with channel metadata. Use `output_data_channels()` when postprocessing tools,
validation dashboards, or GUI code need stable units and axis labels without
hard-coding dataset positions.

The public helpers are:

- `ResultChannel`
- `output_data_channels()`
- `output_data_channel(name)`
- `output_data_channel_names()`
- `annotate_output_data_channels!(file)`
- `WindIORunSpec`
- `windio_run_spec(modeling_options_file, windio_file, run_path)`
- `render_windio_run_script(spec)`
- `write_windio_run_script(path, spec)`

The generated API listing above includes their full docstrings.
