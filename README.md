# OWENS (Offshore Wind ENergy Simulator)

This repository is based on the original structural dynamics solver by Brian Owens (see dissertation: http://hdl.handle.net/1969.1/151813).
The original code has been translated to Julia and revised for simplicity and performance while maintaining accuracy.  GXBeam.jl has also been integrated in beta form for geometrically exact beam solutions
The aerodynamics are provided by the VAWTAero.jl module (https://gitlab.sandia.gov/8821-vawt-tools/VAWTAero.jl) in addition to a coupling to the OpenFAST AeroDyn module.

## Documentation

In Work: Documentation can be found in the docs folder along with the validation paper(s).

## Using/Setting Up the Package
1. Get access to repository and set up ssh keys
	Create a ssh key via:

```bash
ssh-keygen -t rsa -m PEM -C username@sandia.gov
```
	Copy the resulting id_rsa.pub key found in ~./ssh into your Github profile public key

# Mac, to just install OWENS as a regular user, find the shetup.sh script in the test folder.  It will run homebrew to install the required software, build and compile the openfast inflowwind, hydrodyn, and moordyn dynamic libraries, and install the OWENS packages in the correct order.  There is also a VScode profile in the docs/setup folder that can be loaded into VSCode to set up the IDE for julia.
```bash
chmod +x setup.sh
./setup.sh
```

# Linux, same as mac except change the brew statements in the .sh script to apt-get.
```bash
chmod +x setup.sh
./setup.sh
```

# Windows, install the software manually by downloading the windows executables, be sure julia is on your path, and follow the windows compilation instructions for the openfast inflowwind, libraries. At this stage in the software's maturity, it is recommended to use mac or linux environments unless the user is experienced with compiled software development in a windows environment.
    - https://julialang.org/downloads/
    - https://www.paraview.org/download/
    - https://visualstudio.microsoft.com/downloads/
    - https://github.com/OpenFAST/openfast

    Then start julia, either from the command line, or from vscode
    ```bash
    julia
    ```

    And run the following:
    ```julia
    # Install OWENS and non-registered dependencies as a regular user

    using Pkg
    Pkg.add(PackageSpec(url="https://github.com/byuflowlab/Composites.jl.git"))
    Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/PreComp.jl.git"))
    Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/OpenFASTWrappers.jl.git"))
    Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/VAWTAero.jl.git"))
    Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/GyricFEA.jl.git"))
    Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/OWENS.jl.git"))

    # Add other registered packages for running the example scripts
    Pkg.add("PyPlot")
    Pkg.add("Statistics")
    Pkg.add("DelimitedFiles")
    Pkg.add("Dierckx")
    Pkg.add("QuadGK")
    Pkg.add("FLOWMath")
    Pkg.add("HDF5")
    # Build and test OWENS
    Pkg.test("OWENS")
    ```


2. To develop OWENS, the recommended setup is as follows:

    -	Use linux, mac, or unix based system – the code works on windows, but you’ll find yourself taking 10x longer to setup and troubleshoot any custom compiled fortran libraries.

    -	Install Julia 1.7+ (brew/apt-get install, or download the binary from julialang.org, more detail here https://julialang.org/downloads/platform/) such that it is on the system path

    -	Clone the packages you will be developing into an easy to access working folder (e.g. OWENS, GyricFEA, VAWTAero, etc)

    -	If you are using a proxy, be sure that the proxy variables are declared/exported in your .bash_profile or the equivalent
        * http_proxy, https_proxy, HTTP_PROXY, HTTPS_PROXY, no_proxy, NO_PROXY
        * git config --global http.proxy http://user:nopass@proxy.yourorg:number
        * git config --global https.proxy http://user:nopass@proxy.yourorg:number
        * export JULIA_SSL_CA_ROOTS_PATH=""
        * export JULIA_SSL_NO_VERIFY_HOSTS="*.yourorgurl"
        * export JULIA_PKG_USE_CLI_GIT=true 				


    -	Install the custom packages in the following order to prevent precompilation dependency errors.
        * OpenFASTWrappers.jl
        * VAWTAero.jl
        * GyricFEA.jl
        * PreComp.jl
        * Composites.jl (from https://github.com/byuflowlab/Composites.jl.git)
        * OWENS.jl			


    -	Install custom repositories you want to develop by starting Julia from the cloned directory and using the command:
        * ] dev .
        * This type of installation will cause the module to reload each time Julia starts without needing to tell Julia to update 	


    -	Install custom repositories you don’t want to develop directly within Julia using:
        * ] add url2yourrepo
        * You may need to set up ssh keys and use the ssh url depending on the repo setup. To get updates to these packages, you need to tell julia to update, ie ] update thepackage


    -	Install other nonregistered packages using:
        * ] add url2yourpackage 		


    -	Install other registered packages as prompted using:
        * ] add packagename 			


    -	Be sure to run update before building/testing to make sure the other ~120 open source dependencies are up to date.
        * ] update 				


    -	Run your code either from a terminal via julia scriptname.jl or from within an IDE like VSCode where you can see variables, use the debugger with breakpoints similar to matlab and run sections of code iteratively as desired.
        * OWENS can be run from files, like in the test scripts, or without files (except the composite input and airfoils)
        * There are several plotting options in Julia, note that if you use PyPlot it does download anaconda and set up its own conda behind Julia. If you are using a proxy and have your proxy settings right (as described above) it will install smoothly (though it takes a while to download everything initially and when you use it the first time). 	


    -	All of the functions have docstrings describing the i/o and function purpose, which can be accessed by:
        * import module
        * ? module.function() 				
    -	A note about julia debuggers – if you don’t want it to step through everything, you need to tell it what packages to compile vs while packages to step through. This will make the debugger comparable (if not faster) than Matlab in speed. In VSCode, this can be done in the debug pane.
```Julia
## Software License

Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
Export of this data may require a license from the United States Government.

See Copyright.txt file for more information
