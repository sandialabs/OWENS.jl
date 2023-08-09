using OWENS_APP

push!(ARGS, "arg")
OWENS_APP.julia_main()

using Example
Example.hello("PackageCompiler")

# It is ok to use stdlibs that are not in the project dependencies
using Test
@test 1==1