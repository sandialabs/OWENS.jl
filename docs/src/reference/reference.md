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

Use `run_manifest_issues(manifest_or_path)` and
`validate_run_manifest(manifest_or_path)` to check that provenance manifests
retain the expected schema, file-record shape, byte counts, and SHA-256 digest
format before a GUI or validation dashboard trusts them.

Use `verify_file_provenance(record; root)` and
`run_manifest_health(manifest_or_path)` to detect missing, modified, or
malformed run artifacts from a manifest. `run_manifest_health` can also attach
`output_data_summary` rows for healthy HDF5 outputs so a GUI can surface stale
files and missing channel metadata in one pass.

Use `build_output_data_validation_report(reference_path, candidate_path; ...)`
and `write_output_data_validation_report(path, reference_path, candidate_path;
...)` to package selected-channel metric rows into a durable YAML artifact for
regression dashboards. The report intentionally performs elementwise
comparisons only; validation cases that need time alignment, phase correction,
or whole-revolution windows should preprocess the selected channels before
report generation.

Use `build_studio_project(root; ...)`, `studio_project_health(project_or_path)`,
and `render_studio_workbench_html(project_or_health)` as the first OWENS Studio
service boundary. These helpers collect project-level input files, WindIO files,
run manifests, file-provenance checks, and output-channel health into a stable
dictionary schema that the GUI, CLI, or documentation examples can share.
`create_studio_project_template(target; template="rm2")` creates the first
file-backed GUI template: an RM2 WindIO project with copied inputs, a generated
Julia run script, a run manifest, and a Studio project manifest.
Use `studio_project_template_catalog()` to drive GUI "new project" controls
without duplicating template metadata.
The `OWENS_APP` package exposes dependency-light route handlers around the same
services so a future Genie shell can serve health YAML and workbench HTML
without duplicating project logic.
Use `open_studio_project(project_path)` in `OWENS_APP` as the one-call
workbench bootstrap payload after a user selects a project.
Use `studio_route_catalog()` in `OWENS_APP` to keep the future Genie route table
aligned with the tested service handlers.
Use `dispatch_studio_route(route; method, params)` to resolve those catalog
entries by name or path without duplicating route dispatch logic in the web
shell.
Use `studio_project_generated_script_path(project_or_path)` and
`read_studio_project_generated_script(project_or_path)` when the GUI needs to
show or export the exact Julia driver attached to a project.
Use `write_studio_workbench_bundle(output_dir, project_path)` to create a
server-free workbench directory with HTML, health YAML, and the generated Julia
driver. `OWENS_APP.write_studio_project_bundle` also writes the Studio
open-project bootstrap YAML that the app shell can hydrate from.

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
- `run_manifest_issues(manifest_or_path)`
- `validate_run_manifest(manifest_or_path)`
- `verify_file_provenance(record)`
- `run_manifest_health(manifest_or_path)`
- `build_output_data_validation_report(reference_path, candidate_path)`
- `write_output_data_validation_report(path, reference_path, candidate_path)`
- `read_output_data_validation_report(path)`
- `build_studio_project(root)`
- `write_studio_project(path, project_or_root)`
- `read_studio_project(path)`
- `studio_project_issues(project_or_path)`
- `validate_studio_project(project_or_path)`
- `studio_project_health(project_or_path)`
- `studio_project_generated_script_path(project_or_path)`
- `read_studio_project_generated_script(project_or_path)`
- `render_studio_workbench_html(project_or_health)`
- `write_studio_workbench_html(path, project_or_health)`
- `write_studio_workbench_bundle(output_dir, project_path)`
- `studio_project_template_catalog()`
- `studio_project_template_names()`
- `create_studio_project_template(target)`

The generated API listing above includes their full docstrings.
