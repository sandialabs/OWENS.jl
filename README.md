# OWENS (Offshore Wind ENergy Simulator)

[![pipeline status](https://gitlab.sandia.gov/8821-vawt-tools/OWENS.jl/badges/master/pipeline.svg)](https://gitlab.sandia.gov/8821-vawt-tools/OWENS.jl/-/commits/master) [![coverage report](https://gitlab.sandia.gov/8821-vawt-tools/OWENS.jl/badges/master/coverage.svg)](https://gitlab.sandia.gov/8821-vawt-tools/OWENS.jl/-/commits/master)

This repository is based on the original structural dynamics solver by Brian Owens (see dissertation: http://hdl.handle.net/1969.1/151813).
The original code has been translated to Julia and revised for simplicity and performance while maintaining accuracy.
This repository consists of the core of the aero-hydro-servo-elastic solver and structural dynamics solver (elastic).
The aerodynamics are provided by the VAWTAero.jl module (https://gitlab.sandia.gov/8821-vawt-tools/VAWTAero.jl), the hydro by TBD, and the servo by TBD.

## Documentation

The theory manual can be found in the docs folder along with the validation paper(s).
Additionally, there is a lessons learned document regarding Matlab to C++ "automatic" translation.
The test cases are the base material for the validation paper(s).


## Software License

Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
Export of this data may require a license from the United States Government.