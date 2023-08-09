julia -e 'using PackageCompiler; create_app("OWENS_APP","OWENS_APP_COMPILED";force=true)'
OWENS_APP_COMPILED/bin/OWENS_APP foo bar --julia-args -t4