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

## Getting started

1. **Clone repository**: TODO

2. **Install Julia and Julia packages**: TODO

3. **Install Python and Python packages**: The `environment.yml` file lists the working Python version and packages.
To set up your system, you can use the following steps. If you are an adept Python user, you can certainly figure out other ways to do this.

	1. **Install Anaconda**: Anaconda is a tool for managing Python packages and installations (it can also handle non-Python packages). You can download an install from https://www.anaconda.com/products/individual.

	2. **Create a dedicated conda environment & install depedencies**: From your terminal, call the following to create a dedicated [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) with the necessary dependecies.

	```bash
	conda env create -f environment.yml
	```

	3. **Activate conda environment**: Whenever you want to do work on this project, activate the `OWENS` conda environment.

	```bash
	conda activate OWENS
	```

4. **TODO**: TODO

## Software License

Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
Export of this data may require a license from the United States Government.