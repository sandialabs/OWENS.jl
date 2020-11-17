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

## Using the Package (without developing)
1. Get access to repository and set up ssh keys
	Create a ssh key via:

```bash
ssh-keygen -t rsa -m PEM -C username@sandia.gov
```
	Copy the id_rsa.pub key found in ~./ssh into your GitLab profile public key

2. Install Julia such that it is on your system path so it can be called via "julia" in the terminal (see https://julialang.org/downloads/platform/ for platform specific instructions, for mac this is as simple as: brew cask install julia)

3. Install the custom dependencies and OWENS, then build and test using the following command in the terminal:

```bash
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/byuflowlab/OptimizationParameters.jl")); Pkg.add(PackageSpec(url="git@gitlab.sandia.gov:8821-vawt-tools/VAWTAero.jl.git")); Pkg.add(PackageSpec(url="git@gitlab.sandia.gov:8821-vawt-tools/PreComp.jl.git")); Pkg.add(PackageSpec(url="git@gitlab.sandia.gov:8821-vawt-tools/OWENS.jl.git")); Pkg.build("OWENS"); Pkg.test("OWENS";coverage=true)'
```

4. Run your desired cases by: TODO:

Also see example: TODO:

## Developing (not just using the final package)

1. Follow steps 1 and 2 from above to get access to the repo and julia running

2. Clone this repository (into a convenient location of your choice) and enter it

```bash
git clone git@gitlab.sandia.gov:8821-vawt-tools/OWENS.jl.git
cd ./OWENS.jl
```

3. Start Julia and Enter Development Mode

```bash
julia
```
```julia
] activate .
```

4. Build the package (uses the Manifest.toml so you have the exact same build as everyone else)
```julia
(OWENS) pkg> build
```
Note: If you're installing julia for the first time, there are quite a few dependencies of dependencies and it will take some time to install it all

5. Run tests to verify it built correctly
```julia
(OWENS) pkg> test
```

6. Make changes/additions to you heart's desire (following general best practices of course!)
To exit the package manager (but retain the OWENS environment), press backspace.  From this julia console, you'll be able to run any file e.g. include("./test/runtests.jl") with the OWENS environment active.  You can run OWENS with the OWENS environment deactivated via either 1) Installing it as a backend package, or 2) adding all of the dependencies to your base julia environment (the environment that we were in before activating the current directory in step 3.).  

If you'd like a more interactive IDE and Matlab-like debugger, consider Juno.  Here are setup instructions: http://docs.junolab.org/stable/man/installation/

There are a TON of questions already answered via a simple google search about Julia code development etc, but also feel free to reach out to the developers of this repository for julia specific help.  

## Software License

Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
Export of this data may require a license from the United States Government.

See Copyright.txt file for more information
