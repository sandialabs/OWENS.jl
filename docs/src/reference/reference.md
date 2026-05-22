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

Use `output_data_summary(path)` to inspect which registered channels are
present in an HDF5 result file, including shape, element type, units, and
metadata-attribute mismatches, without reading the full result arrays. This is
the preferred entry point for GUI channel pickers and validation dashboards.

Use `output_data_channel_metrics(reference_path, candidate_path; channels, atol,
rtol)` to compare selected numeric channels from two result files. It reports
bias, RMSE, maximum absolute error, mean absolute error, reference RMS, relative
RMSE, and a pass/fail result using `atol + rtol * max_abs_reference`.

The public helpers are:

- `ResultChannel`
- `output_data_channels()`
- `output_data_channel(name)`
- `output_data_channel_names()`
- `annotate_output_data_channels!(file)`
- `output_data_summary(path)`
- `output_data_channel_metrics(reference_path, candidate_path)`
- `WindIORunSpec`
- `windio_run_spec(modeling_options_file, windio_file, run_path)`
- `render_windio_run_script(spec)`
- `write_windio_run_script(path, spec)`
- `build_windio_run_manifest(spec)`
- `write_windio_run_manifest(path, spec)`

The generated API listing above includes their full docstrings.
