
#!/usr/bin/env bash

# This script assumes you already have the required ssh keys and proxy (if needed) all set up.
# Also assumes you have the necessary fortran compilers
# xcode-select -install
# brew install cmake
# brew install gfortran

# Update brew packages

# Install julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.3-linux-x86_64.tar.gz
tar zxvf julia-1.9.3-linux-x86_64.tar.gz

# Load required modules to compile openfast libraries
module load sparc-cmake/3.23.2
module load sparc-dev/gcc  
module load struct-mech-tools/python3.9
module load paraview

# Install openfast coupled libraries !NOTE!: if you change the location of the compiled libraries, you may need to update the rpath variable, or recompile.
cd ../../
git clone --depth 1 https://github.com/OpenFAST/openfast.git
mkdir openfast/build
cd openfast/build
cmake -DBUILD_SHARED_LIBS=ON ..
make ifw_c_binding
make moordyn_c_binding
make hydrodyn_c_binding
make aerodyn_inflow_c_binding
cd ../../

# Install OWENS and non-registered dependencies as a regular user
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/byuflowlab/Composites.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/PreComp.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/OpenFASTWrappers.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/VAWTAero.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/GyricFEA.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/ModelGen.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/OWENS.jl.git"))'

# Add other registered packages for running the example scripts
julia -e 'using Pkg;Pkg.add("PyCall");ENV["PYTHON"] = "/projects/struct_mech_tools/bin/python3.9/bin/python";
 Pkg.add("PyPlot");Pkg.add("Statistics");Pkg.add("DelimitedFiles");Pkg.add("Dierckx");Pkg.add("QuadGK");Pkg.add("FLOWMath");
 Pkg.add("HDF5");Pkg.add("GXBeam")'

using PyCall
pygui(:qt5)
using PyPlot


# Run the example script
julia -e 'using PyCall; pygui(:qt5); include("ExampleSNL5MW_turbulent.jl")'

#for VPM 
#module load sparc-dev/intel-2021.3.0_mpich2-3.2
# follow ed's instructions, but must remove compat for CXXWrap, run update in julia, and then dev . for FLOWExaFMM
# clone GXBeam, remove implicitAD compat, add main ImplicitAD, dev GXBeam



# Install Paraview
# brew install --cask paraview

# Install visual studio
# brew install --cask visual-studio-code