
#!/usr/bin/env bash

# This script assumes you already have the required ssh keys and proxy (if needed) all set up.
# Also assumes you have the necessary fortran compilers
# xcode-select -install
# brew install cmake
# brew install gfortran

# Update brew packages
brew upgrade

# Install julia
brew install --cask julia

# Install openfast coupled libraries !NOTE!: if you change the location of the compiled libraries, you may need to update the rpath variable, or recompile.
git clone --depth 1 https://github.com/OpenFAST/openfast.git
mkdir openfast/build
cd openfast/build
cmake -DBUILD_SHARED_LIBS=ON ..
make ifw_c_binding
make moordyn_c_binding
make hydrodyn_c_binding
cd ../../

# Install OWENS and non-registered dependencies as a regular user
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/byuflowlab/Composites.jl.git")); Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/PreComp.jl.git")); Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/OpenFASTWrappers.jl.git")); Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/VAWTAero.jl.git")); Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/GyricFEA.jl.git")); Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/OWENS.jl.git")); Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/ModelGen.jl.git"))'

# Add other registered packages for running the example scripts
julia -e 'using Pkg; Pkg.add("PyPlot");Pkg.add("Statistics");Pkg.add("DelimitedFiles");Pkg.add("Dierckx");Pkg.add("QuadGK");Pkg.add("FLOWMath");Pkg.add("HDF5")'
# Build and test OWENS
julia -e 'using Pkg; Pkg.test("OWENS")'

# Run the example script
julia ExampleSNL5MW_turbulent.jl

# Install Paraview
brew install --cask paraview

# Install visual studio
brew install --cask visual-studio-code