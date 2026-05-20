# OWENS GUI Plan

Updated: 2026-05-20

This plan describes a QBlade-class OWENS GUI that is intentionally better in workflow: every GUI action must be reproducible as Julia code, WindIO/YAML, or a checked project manifest; every simulation must carry provenance and validation status; and every advanced model option must be traceable to a solver API, test, or documented limitation.

## Product Goal

Build `OWENS Studio`: a Julia-native engineering workbench for VAWT, HAWT, floating, aeroelastic, and validation workflows. It should compete directly with QBlade on setup, visualization, controller/wind/wave workflows, simulation replay, and postprocessing, while beating QBlade on reproducibility, scripting, validation dashboards, package interoperability, and research extensibility.

The GUI should not become a separate simulation engine. It should orchestrate OWENS, OWENSAero, OWENSFEA, OWENSOpenFASTWrappers, OWENSPreComp, CCBlade, GXBeam, WindIO, OpenFAST modules, and postprocessing APIs through explicit, tested service boundaries.

## Recommended Julia Stack

Primary path:

- Julia backend: a new package/app layer, tentatively `OWENSStudio.jl`, in `OWENS.jl/app/OWENSStudio` or a sibling package once stable.
- UI/server: Genie.jl for the app shell, routing, upload/download endpoints, WebSockets, job state, and project APIs. Genie is a full-stack Julia web framework where controllers are Julia modules and apps are normal Julia projects.
- Visualization: Makie recipes as the core plotting/3D abstraction, with GLMakie for desktop interactive views, WGLMakie for browser/WebGL views, and CairoMakie for publication/vector export. Makie provides interactive visualization, 2D/3D plots, UI blocks, scenes, events, layouts, mesh, volume, surface, streamplot, and animations.
- Native desktop packaging: initially run as a local web app launched by Julia. Later package with a sysimage and lightweight desktop shell if needed.
- Native fallback: Gtk4.jl is viable for a pure native Julia desktop shell, but it should not be the first implementation unless browser/WebGL embedding proves blocking. Gtk4 has Julia bindings and normal desktop widgets, but the OWENS workflow needs rich plots, replay, docs links, uploads, validation dashboards, and web-style layouts.

Why this stack:

- Browser UI gives modern workflow design, tables, trees, drag/drop, docs links, report export, and easy remote/HPC access.
- Julia owns model state, units, validation, and solver calls.
- Makie gives one visualization recipe layer for GUI, docs, tests, and publication outputs.
- A web app can later be deployed for remote runs without rewriting the model layer.

Useful ecosystem anchors:

- Makie docs: https://docs.makie.org/stable/
- Genie docs: https://genieframework.github.io/Genie.jl/dev/tutorials/1--Overview.html
- Gtk4.jl docs: https://juliagtk.github.io/Gtk4.jl/dev/

## Core Design Principle

Every user-facing operation has four artifacts:

1. A validated in-memory Julia model.
2. A serialized project representation.
3. A generated script or command path that can be rerun headlessly.
4. A result/provenance record that includes package versions, git SHAs, solver settings, input hashes, output hashes, warnings, and validation metrics.

This is the biggest opportunity to beat QBlade: QBlade is highly capable, but OWENS can make the GUI a transparent front end for reproducible research rather than a black-box project file workflow.

## User Workflow

The first screen should be the workbench, not a landing page.

Primary navigation:

- Project
- Geometry
- Airfoils
- Structure
- Environment
- Controls
- Simulation
- Validation
- Results
- Reports

### 1. Project Hub

Purpose: Create, open, duplicate, import, validate, and compare complete turbine projects.

Expected controls:

- New project from templates: RM2 VAWT, SNL 5MW HAWT, 34m VAWT, CCT2 floating, custom blank.
- Import from WindIO, OWENS YAML, OpenFAST/AeroDyn/HydroDyn/MoorDyn files, and QBlade text exports where feasible.
- Project health panel: missing files, invalid units, unsupported solver paths, stale generated files, failing validation metrics.
- One-click script export: generate a Julia driver that exactly reproduces the GUI setup.
- Dependency/status panel: package versions, OpenFAST artifact status, optional XFOIL availability, GPU/WebGL availability.

Better than QBlade:

- The GUI always shows the equivalent Julia call path.
- Project state can be reviewed in git.
- Validation failures are visible before expensive runs.

### 2. Geometry And Mesh Builder

Purpose: Build HAWT/VAWT geometry with explicit frames, coordinate conventions, and mesh ownership.

Features:

- HAWT and VAWT mode selector.
- Blade, tower, strut, shaft, platform, cable, and point-mass object tree.
- Interactive 3D model with frame axes, node IDs, element IDs, master/slave joints, span/height/chord/twist labels.
- Segment and node-generation tools with deterministic IDs.
- Explicit global/local frame inspector for every component.
- Helical blade preview and AeroDyn export preview.
- Mesh quality checks: duplicate points, overconstrained joints, unexpected span/height conventions, non-monotone stations, invalid element lengths.
- Export to OWENS, GXBeam, AeroDyn, and VTK.

Non-GUI prerequisite:

- A stable project/geometry schema and transform API that separates stationary tower/platform frames from rotating rotor frames.

### 3. Airfoil And Polar Lab

Purpose: Manage airfoils, polars, extrapolation, dynamic-stall metadata, and structural moment transfer.

Features:

- Import airfoil coordinates and validate orientation, trailing-edge closure, flat-back handling, chord normalization, and thickness.
- Generate XFOIL polars when Xfoil.jl is available, but keep XFOIL optional.
- Read/write OWENSAero and AeroDyn polar formats.
- Multi-Reynolds polar library with alpha/Re extrapolation guards.
- 360-degree extrapolation workflow.
- Static and dynamic-stall model comparison panels.
- `Cm25`/moment coefficient visualization and structural moment propagation checks.
- Damaged-blade polar regions and active flow-control polar-state regions as future advanced features.

Better than QBlade:

- Polar format conversion is logged and testable.
- Every extrapolation/clamping event can be surfaced as a validation warning.
- Dynamic-stall fixtures and reference metrics are attached to the polar/model.

### 4. Structure And Materials Lab

Purpose: Configure structural properties, composite sections, mass/stiffness matrices, constraints, and modal checks.

Features:

- NuMAD/PreComp-style laminate workflow.
- GXBeam/OWENSFEA section-property import/export.
- Timoshenko, GXBeam, and ROM solver selectors with documented limitations.
- Concentrated mass and platform mass-offset matrix builder.
- Constraint/joint visual inspection.
- Gravity direction and frame checks.
- Modal/Campbell analysis dashboard.
- Stress, fatigue, buckling, Goodman, rainflow, and safety-factor output setup.
- Material/source metadata and unit validation.

Non-GUI prerequisite:

- Public APIs for structural model construction, modal runs, and stress/safety-factor output that do not require examples to mutate globals or parse side-effect files.

### 5. Environment Lab

Purpose: Build wind, wave, current, and offshore inputs with validation and visualization.

Features:

- Uniform, prescribed, WindIO, TurbSim, and InflowWind sources.
- TurbSim template editor with exact-key rendering, `.bts` import, preview, probes, and cut planes.
- Wind shear, turbulence, gust, directional shear, and inflow-angle controls.
- Linear wave generator and HydroDyn input staging.
- Current profile editor.
- MoorDyn line/cable layout preview.
- Seabed/contact and Morison member authoring as long-term native OWENS features.
- PotFile root resolution status and relative-path diagnostics.

Better than QBlade:

- Environment inputs can be round-tripped through WindIO/OpenFAST paths.
- GUI previews share the same sampling code as tests/docs.

### 6. Controls And Actuators Lab

Purpose: Configure generator, drivetrain, controller, pitch/yaw/RPM schedules, and actuator studies.

Features:

- Prescribed RPM/Vinf controls with monotonic time checks.
- Simple generator and drivetrain parameter validation.
- Controller plugin interface for ROSCO/OpenFAST-style controllers and future user DLL/shared-library adapters.
- Command/measurement channel table with units, sign conventions, and validation.
- Sinusoidal step/ramp actuator studies for pitch, RPM, torque, yaw, platform motion, and active flap states.
- Controller logs and output channel plots.

Non-GUI prerequisite:

- A stable controller API with pinned channel names, units, frame/sign conventions, and scriptable run hooks.

### 7. Simulation Lab

Purpose: Configure, launch, monitor, cancel, restart, and compare simulations.

Features:

- Solver path selector: DMS, AC, CCBlade HAWT, AeroDyn/OpenFAST wrapper, HydroDyn/MoorDyn coupled, structures-only, modal, validation script.
- Time-step and convergence controls with estimated runtime and stability warnings.
- Local run queue with logs, progress, warnings, and artifact output paths.
- Optional remote/HPC backend later.
- Restart files and failure-capture output, including VTK dump on failure.
- Deterministic run manifest: inputs, hashes, environment variables, package versions, git SHAs.
- Side-by-side scenario comparison.

Better than QBlade:

- The run queue can generate CI-ready commands.
- Each warning links to the exact model field or source file that caused it.

### 8. Validation Lab

Purpose: Make research credibility visible in the product.

Features:

- Built-in RM2, SNL 5m/34m, SNL 5MW, CCT2 floating, Campbell, dynamic-stall, and HAWT AeroDyn comparison cases.
- Metrics: RMSE, bias, max error, peak error, phase error, nondimensional power/thrust/torque, modal MAC, frequency error, load-channel sign checks.
- Acceptance thresholds by validation case.
- "Validated for this solver path" badge on models/runs.
- Regression fixture browser and diff viewer.
- OpenFAST/AeroDyn station-level comparison for HAWT CCBlade path.
- 2022 validation dashboard with simulated vs experimental overlays.

Better than QBlade:

- The GUI should show exactly where OWENS is validated, where it is only regression-pinned, and where it is unsupported.

### 9. Results, Replay, And Reports

Purpose: Analyze results without leaving the OWENS workflow.

Features:

- Time-series explorer with derived channels.
- Whole-revolution averaging tools.
- PSD, spectra, Campbell, modal, and ensemble plots.
- 3D replay of structure, deformation, loads, wake diagnostics, windfield slices, wave surface, and mooring/platform state.
- VTK, CSV, HDF5/JLD2, MAT, and report export.
- Stress/safety-factor surface maps.
- Comparison reports against validation baselines.
- Publication-quality figure presets.

Non-GUI prerequisite:

- Result stores should be queryable by channel metadata rather than hard-coded array positions.

## Feature Parity Matrix

| Capability | QBlade has | OWENS status | OWENS Studio target |
| --- | --- | --- | --- |
| HAWT BEM/UBEM | Mature | CCBlade path emerging | CCBlade native, AeroDyn validation |
| VAWT DMS/AC | Present | Strong OWENS focus | Preserve, validate, expose clearly |
| Free vortex/LLFVW | Strong | Gap/FLOWUnsteady candidate | Add validated wake/replay workflow |
| Airfoil/polar generation | Strong | Improving | Better metadata, validation, script export |
| Active flow control | Present | Gap | Add after polar data model stabilizes |
| Damaged-blade polars | Present | Gap | Add as alternate-polar regions |
| Structural dynamics | Chrono workflow | OWENSFEA/GXBeam/ROM | Solver-transparent, validation-first UI |
| Offshore hydro | Native authoring | HydroDyn/MoorDyn wrappers | Wrapper-first, native authoring later |
| Controllers | Bundled external controllers | Basic generator/drivetrain | Stable controller API and examples |
| Windfield authoring | Strong GUI | TurbSim/OpenFAST helpers | TurbSim/WindIO diagnostics and previews |
| Noise/ice throw | Present | Gap | Lower-priority modules |
| Replay/cut planes | Strong | Partial VTK/postprocessing | First-class result explorer |
| Validation transparency | Limited | Tests/scripts | User-visible validation dashboards |
| Script reproducibility | Limited | Strong potential | Mandatory script export for every workflow |

## Architecture

Recommended package/app split:

- `OWENSProject`: project schema, units, file references, input validation, migrations.
- `OWENSStudio`: Genie app, route handlers, UI state, auth/security if remote.
- `OWENSRunner`: local/remote job runner, cancellation, restart, logs, progress, artifact manifest.
- `OWENSViz`: Makie recipes for geometry, meshes, sections, frames, polars, time histories, validation comparisons, stress maps, and replay.
- `OWENSResults`: channel registry, HDF5/JLD2 store, metadata, provenance, derived quantities, export.
- `OWENSValidation`: validation case registry, metrics, thresholds, baseline management.
- `OWENSController`: controller API, schedule/actuator commands, ROSCO/OpenFAST-style adapter path.
- `OWENSPlugins`: plugin/extension contract for new solvers, templates, validation cases, and report generators.

These can start as modules under `app/OWENSStudio` and graduate to packages only when APIs stabilize.

## Data Model

Use a project manifest that can be reviewed in git:

- `Project.toml`: Julia environment for the project.
- `owens_project.yml`: canonical model, file references, units, solver choices.
- `inputs/`: WindIO, OWENS, OpenFAST, airfoil, structure, wave/current, controller inputs.
- `generated/`: generated OpenFAST/AeroDyn/HydroDyn/MoorDyn/TurbSim files.
- `runs/<run_id>/run_manifest.yml`: hashes, versions, commands, warnings, outputs.
- `runs/<run_id>/results.h5`: channelized results.
- `runs/<run_id>/figures/`: generated plots.
- `runs/<run_id>/vtk/`: VTK outputs.
- `validation/<case_id>/metrics.yml`: validation metrics and thresholds.

The GUI should never silently overwrite user input files. Generated files go under `generated/` or the run directory with hashes.

## UI Design Rules

- Dense engineering UI, not a marketing page.
- Left navigation for workflow stages; center canvas/table editor; right inspector/validation panel.
- Every numerical field shows units and validation state.
- Invalid entries fail early with specific messages.
- Use frame/axis icons, not prose, for orientation controls.
- Geometry and results views should be interactive but not decorative.
- Tooltips explain conventions, not basic UI.
- Provide "show generated Julia" and "show generated files" from every stage.
- The app should work on a laptop for setup and postprocessing, while heavier runs can be queued locally or remotely later.

## Implementation Phases

### Phase 0 - Non-GUI Foundations

Goal: make existing OWENS APIs GUI-safe.

Tasks:

- Define project schema and result manifest.
- Create pure constructors/validators for geometry, airfoils, environment, controls, and solver options.
- Add channel metadata for common outputs.
- Add Makie recipe skeletons for geometry, polars, and time histories.
- Add validation registry format.
- Add run manifest/provenance writer.

Exit criteria:

- A simple RM2 project can be loaded, validated, run headlessly, and produce a result manifest without GUI-specific code.

### Phase 1 - Workbench Shell

Goal: first usable local GUI.

Tasks:

- Genie app shell with project browser.
- Upload/import/open existing OWENS/WindIO project.
- Project health panel.
- Geometry view using Makie.
- Airfoil/polar table and plot.
- Run button for one existing example.
- Results plot for CP/torque/power.
- Script export for the run.

Exit criteria:

- A user can open RM2, inspect geometry/polars, run a short case, and view results.

### Phase 2 - Model Builders

Goal: QBlade-like setup capability for VAWT and HAWT.

Tasks:

- Component tree editor for blade/tower/strut/platform/cable.
- Mesh/joint visual diagnostics.
- Polar generation/conversion workflow.
- Structure/material section editor.
- HAWT CCBlade setup path.
- WindIO/OpenFAST import/export checks.

Exit criteria:

- New VAWT and HAWT projects can be created without manually editing files.

### Phase 3 - Run Management And Validation

Goal: make the GUI credible for research.

Tasks:

- Run queue, cancellation, logs, warnings, restart artifacts.
- Validation dashboard.
- RM2/SNL/CCT2/HAWT validation cases.
- Whole-revolution averaging and metric overlays.
- Report export.

Exit criteria:

- A model can be marked validated/regression-pinned/unsupported by solver path.

### Phase 4 - Advanced Physics Parity

Goal: close QBlade capability gaps.

Tasks:

- Controller API and actuator studies.
- Windfield probes/cut planes.
- Wake replay/cut planes through FLOWUnsteady/OpenFAST/OWENSAero hooks.
- Offshore environment authoring around HydroDyn/MoorDyn.
- Active flow-control and damaged-blade polar regions.
- Native Morison/environment preprocessing if wrapper workflows are not enough.

Exit criteria:

- OWENS Studio can set up and inspect the major workflows users currently choose QBlade for.

### Phase 5 - Distribution

Goal: make it easy to install and run.

Tasks:

- PackageCompiler sysimage for lower startup latency.
- One-command local app launcher.
- Binary/artifact checks for OpenFAST wrappers.
- Optional desktop wrapper.
- Example project gallery.
- GUI smoke tests and screenshot tests in CI where possible.

Exit criteria:

- A new user can install, launch, open examples, and run a short validation case without hand-editing Julia environments.

## Testing Strategy

- Unit tests for every validator and project-schema migration.
- Golden project fixtures for RM2, SNL 5MW HAWT, CCT2 floating, and one blank/custom case.
- Headless GUI API tests for route handlers.
- Makie recipe tests that check generated objects/data, not fragile pixels, plus a few screenshot tests for layout regressions.
- Round-trip tests: GUI model -> serialized project -> generated Julia -> same solver inputs.
- Result manifest tests with exact hashes for small fixtures.
- Validation metric tests with pinned values and thresholds.

## Immediate Next Tasks

Current implementation branch stack:

1. Complete on PR #141-#148: output channel registry, channel summaries, metrics, validation reports, run manifests, manifest schema checks, and run-manifest health checks.
2. In progress on `codex/gui-project-health-shell`: dependency-light OWENS Studio project manifest plus static health/workbench HTML and `app/OWENS_APP` service commands.

Updated task breakdown:

| Step | Task | State | Acceptance |
| --- | --- | --- | --- |
| G1 | Add `OWENSProject`/Studio project schema prototype for project files and run manifests. | Complete on `codex/gui-project-health-shell` | `build_studio_project`, YAML round trip, schema diagnostics, and project health tests. |
| G2 | Expose first workbench health panel without Genie/Makie. | Complete on `codex/gui-project-health-shell` | Static HTML renders project status, file/run records, schema issues, and can be generated from a project path. |
| G3 | Add app service commands for manifest health, output summaries, WindIO script prep, project health, and static HTML export. | Complete on `codex/gui-project-health-shell` | `app/OWENS_APP` no longer calls stale `runOWENS`; commands return structured YAML-compatible data. |
| G4 | Add a local Genie shell that calls the service layer. | In progress on `codex/gui-project-health-shell` | Open-project bootstrap payloads, static workbench bundles that link `open.yml`, route handlers, and a route catalog are now tested service contracts for Genie to wrap next. |
| G5 | Add `OWENSViz` Makie recipes for geometry/frames, airfoil/polar curves, and time-history/validation overlays. | Pending | Recipe unit tests on generated plot data; screenshot checks only after layout stabilizes. |
| G6 | Wire RM2 as the first end-to-end GUI workflow. | In progress on `codex/gui-project-health-shell` | Template catalog, RM2 creation, health, generated Julia script, and static bundle workflows are pinned in tests. |
| G7 | Add "export generated Julia script" before adding more UI controls. | In progress on `codex/gui-project-health-shell` | WindIO run specs render scripts; RM2 Studio templates write a script artifact; app services and routes now expose that script directly. |
| G8 | Convert the static workbench shell into the first Genie route group. | In progress on `codex/gui-project-health-shell` | Dependency-light route handlers, route parameter metadata, a dispatcher, and a pinned route catalog now define project health, workbench, script, bundle, and template endpoints for Genie to wrap. |
| G9 | Add first GUI fixture project under `examples/gui` or `docs/src/assets/gui`. | In progress on `codex/gui-project-health-shell` | RM2 template creation is deterministic in tests; a committed fixture can be added once the schema stabilizes. |

Near-term rule: keep Genie and Makie out until G1-G3 are merged. The service boundary must be stable and testable first; the UI framework should call it rather than own solver or file semantics.

## Open Risks

- Julia package startup latency can make a GUI feel slow; sysimages and persistent sessions are likely required.
- WGLMakie/browser 3D may lag GLMakie for heavy mesh/wake replay; the visualization layer should support both.
- Existing OWENS APIs still mix file side effects, globals, and example-specific conventions; Phase 0 is mandatory.
- OpenFAST native library setup remains a user-support risk until artifact checks and install docs are fully robust.
- A GUI without validation dashboards would only copy QBlade's surface value; the validation layer is not optional.

## Definition Of Done For A Competitive Release

- VAWT and HAWT project setup can be completed from the GUI.
- RM2, SNL 5MW HAWT, and one floating case run from templates.
- Airfoil/polar management supports import, validation, conversion, plotting, and script export.
- Structure/material setup supports section properties, modal checks, and stress output.
- Wind/wave/current/controller setup is visible, validated, and reproducible.
- Runs have manifests, logs, warnings, result stores, and replay.
- Validation dashboards show metrics against baselines.
- Every GUI action is reproducible as Julia code or WindIO/YAML.
- Documentation includes tutorial videos/screenshots plus generated scripts.
