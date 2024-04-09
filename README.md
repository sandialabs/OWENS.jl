# OWENS (Onshore/Offshore Wind/Water ENergy Simulator)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sandialabs.github.io/OWENS.jl)
![](https://github.com/sandialabs/OWENS.jl/workflows/CI/badge.svg)

This repository is based on the original structural dynamics solver by Brian Owens (see dissertation: http://hdl.handle.net/1969.1/151813).
The original code has been translated to Julia and revised for simplicity and performance while maintaining accuracy.  GXBeam.jl has also been coupled for geometrically exact beam solutions
The aerodynamics are provided by the OWENSAero.jl module (https://gitlab.sandia.gov/8821-vawt-tools/OWENSAero.jl) in addition to a coupling to the OpenFAST AeroDyn module. All codes that can be standalone (like the aerodynamics and structures) should be separate and handled through the dependency manager, other functions specific to the OWENS ontology should remain in this repository.

## Examples
Please see the documentation under examples.

## Installation 
Please see the documentation under setup. 

## Contributing
Please make all feature changes and bug fixes as branches and then create pull requests against the dev branch.  The dev branch will be periodically pulled into master for significant version changes.
