# OWENS (Offshore Wind ENergy Simulator)

This repository is based on the original structural dynamics solver by Brian Owens (see dissertation: http://hdl.handle.net/1969.1/151813).
The original code has been translated to Julia and revised for simplicity and performance while maintaining accuracy.  GXBeam.jl has also been integrated in beta form for geometrically exact beam solutions
The aerodynamics are provided by the VAWTAero.jl module (https://gitlab.sandia.gov/8821-vawt-tools/VAWTAero.jl) in addition to a coupling to the OpenFAST AeroDyn module.

## Documentation

In Work: Documentation can be found in the docs folder along with the validation paper(s).

-	All of the functions have docstrings describing the i/o and function purpose, which can be accessed by:
    * import module
    * ? module.function() 				
-	A note about julia debuggers – if you don’t want it to step through everything, you need to tell it what packages to compile vs while packages to step through. This will make the debugger comparable (if not faster) than Matlab in speed. In VSCode, this can be done in the debug pane.

## Software License

Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
Export of this data may require a license from the United States Government.

See Copyright.txt file for more information
