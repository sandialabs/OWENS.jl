julia -e 'using PackageCompiler; create_app("OWENS_APP","OWENS_APP_COMPILED";force=true)'
OWENS_APP_COMPILED/bin/OWENS_APP "OWENS_APP/src/data/SNL34mGeom.csv" bar --julia-args -t4