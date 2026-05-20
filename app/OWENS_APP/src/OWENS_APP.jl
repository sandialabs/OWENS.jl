module OWENS_APP
import OWENS
import YAML
import OrderedCollections

include("studio_services.jl")
include("studio_routes.jl")

const USAGE = """
OWENS_APP Studio service commands:

  project-routes
  project-templates
  project-open <owens_project.yml>
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
    elseif command == "project-templates" && length(args) == 1
        return list_studio_project_templates()
    elseif command == "project-open" && length(args) == 2
        return open_studio_project(args[2])
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

end #module
