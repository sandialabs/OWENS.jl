# OWENS (Offshore Wind ENergy Simulator)

This repository is based on the original structural dynamics solver by Brian Owens (see dissertation: http://hdl.handle.net/1969.1/151813).
The original code has been translated to Julia and revised for simplicity and performance while maintaining accuracy.  GXBeam.jl has also been coupled for geometrically exact beam solutions
The aerodynamics are provided by the VAWTAero.jl module (https://gitlab.sandia.gov/8821-vawt-tools/VAWTAero.jl) in addition to a coupling to the OpenFAST AeroDyn module. All codes that can be standalone (like the aerodynamics and structures) should be separate and handled through the dependency manager, other functions specific to the OWENS ontology should remain in this repository.

## Documentation
Until public hosting of the documentation is set up, a readthedocs style webpage can be built via:

    cd path2OWENS.jl/OWENS.jl/docs
    julia --project make.jl

and then a local server can be started via

    cd ..
    julia -e 'using LiveServer; serve(dir="docs/build")'

then open your favorite browser and open the following (or what is indicated in the terminal output if different)

    http://localhost:8000/

Note that you may need to install julia packages as directed in the terminal output if you are building the docs for the first time.

## Contributing
Please make all feature changes and bug fixes as branches and then create pull requests against the dev branch.  The dev branch will be periodically pulled into master for significant version changes.

## Software License

Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
Export of this data may require a license from the United States Government.

See Copyright.txt file for more information
