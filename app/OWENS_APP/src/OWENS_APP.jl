module OWENS_APP
import OWENS
import YAML
import OrderedCollections

include("studio_services.jl")
include("studio_routes.jl")

const USAGE = """
OWENS_APP Studio service commands:

  project-routes
  studio-doctor [output_dir]
  studio-home <output.html>
  project-capabilities
  project-quality-gates
  project-templates
  project-examples
  project-open <owens_project.yml>
  project-inputs <owens_project.yml> [include_text] [max_text_bytes]
  project-session <owens_project.yml> [include_text] [max_text_bytes]
  project-editor-html <owens_project.yml> <output.html>
  project-save-input <owens_project.yml> <role> <text_file>
  manifest-health <run_manifest.yml>
  output-summary <results.h5>
  windio-script <modeling_options.yml> <windio.yml> <run_path>
  project-template <template> <target_dir>
  project-health <owens_project.yml>
  project-script <owens_project.yml>
  project-html <owens_project.yml> <output.html>
  project-bundle <owens_project.yml> <output_dir>
"""

function julia_main()::Cint
    try
        real_main(ARGS)
        return 0
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
end

function real_main(args = ARGS; io = stdout)
    result = _dispatch_command(args)
    YAML.write(io, result)
    write(io, "\n")
    return result
end

function _dispatch_command(args)
    isempty(args) && throw(ArgumentError(USAGE))
    command = args[1]

    if command == "project-routes" && length(args) == 1
        return studio_route_catalog()
    elseif command == "studio-doctor" && length(args) in (1, 2)
        output_dir = length(args) == 2 ? args[2] : pwd()
        return studio_doctor(; output_dir)
    elseif command == "studio-home" && length(args) == 2
        html = write_studio_home(args[2])
        return OrderedCollections.OrderedDict{String,Any}(
            "output_html" => abspath(args[2]),
            "bytes" => sizeof(html),
        )
    elseif command == "project-capabilities" && length(args) == 1
        return list_studio_gui_capabilities()
    elseif command == "project-quality-gates" && length(args) == 1
        return list_studio_quality_gates()
    elseif command == "project-templates" && length(args) == 1
        return list_studio_project_templates()
    elseif command == "project-examples" && length(args) == 1
        return list_studio_example_projects()
    elseif command == "project-open" && length(args) == 2
        return open_studio_project(args[2])
    elseif command == "project-inputs" && length(args) in 2:4
        include_text = length(args) >= 3 ? _cli_bool(args[3]) : false
        max_text_bytes = length(args) == 4 ? parse(Int, args[4]) : 200_000
        return inspect_studio_project_inputs(args[2]; include_text, max_text_bytes)
    elseif command == "project-session" && length(args) in 2:4
        include_text = length(args) >= 3 ? _cli_bool(args[3]) : false
        max_text_bytes = length(args) == 4 ? parse(Int, args[4]) : 200_000
        return inspect_studio_project_session(args[2]; include_text, max_text_bytes)
    elseif command == "project-editor-html" && length(args) == 3
        html = write_studio_project_editor(args[3], args[2])
        return OrderedCollections.OrderedDict{String,Any}(
            "output_html" => abspath(args[3]),
            "bytes" => sizeof(html),
        )
    elseif command == "project-save-input" && length(args) == 4
        return save_studio_project_input(args[2], args[3], read(args[4], String))
    elseif command == "manifest-health" && length(args) == 2
        return inspect_run_manifest(args[2])
    elseif command == "output-summary" && length(args) == 2
        return inspect_output_data(args[2])
    elseif command == "windio-script" && length(args) == 4
        return prepare_windio_run(args[2], args[3], args[4])
    elseif command == "project-template" && length(args) == 3
        return create_studio_template_project(args[3]; template = args[2])
    elseif command == "project-health" && length(args) == 2
        return inspect_studio_project(args[2])
    elseif command == "project-script" && length(args) == 2
        return inspect_studio_project_script(args[2])
    elseif command == "project-html" && length(args) == 3
        project_health = inspect_studio_project(args[2])
        html = write_studio_project_workbench(args[3], project_health)
        return OrderedCollections.OrderedDict{String,Any}(
            "output_html" => abspath(args[3]),
            "bytes" => sizeof(html),
            "project_status" => project_health["status"],
        )
    elseif command == "project-bundle" && length(args) == 3
        return write_studio_project_bundle(args[3], args[2])
    end

    throw(ArgumentError(USAGE))
end

function _cli_bool(value)
    normalized = lowercase(strip(string(value)))
    normalized in ("true", "1", "yes", "y") && return true
    normalized in ("false", "0", "no", "n") && return false
    throw(ArgumentError("Expected boolean CLI value, got: $value"))
end

end #module
