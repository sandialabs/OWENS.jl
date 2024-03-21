
#!/usr/bin/env bash

# Install julia
curl -fsSL https://install.julialang.org | sh

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
make aerodyn_driver
make turbsim
cd ../../

julia -e 'using Pkg;Pkg.develop(url="../../GXBeam.jl")'

# Add registered packages for running the example scripts
julia -e 'using Pkg;Pkg.add("PyCall");ENV["PYTHON"] = "/projects/struct_mech_tools/bin/python3.9/bin/python";
 Pkg.add("PyPlot");Pkg.add("Statistics");Pkg.add("DelimitedFiles");Pkg.add("Dierckx");Pkg.add("QuadGK");Pkg.add("FLOWMath");
 Pkg.add("HDF5")'

# Install OWENS and non-registered dependencies as a regular user
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/byuflowlab/Composites.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/OWENSPreComp.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/OWENSOpenFASTWrappers.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/OWENSAero.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/OWENSFEA.jl.git")); Pkg.add(PackageSpec(url="git@cee-gitlab.sandia.gov:8921-VAWT-TOOLS/OWENS.jl.git"))'

using PyCall
pygui(:qt5)
using PyPlot


# Run the example script
julia -e 'using PyCall; pygui(:qt5); include("ExampleSNL5MW_turbulent.jl")'

#for VPM 
#module load sparc-dev/intel-2021.3.0_mpich2-3.2
# follow ed's instructions, but must remove compat for CXXWrap, run update in julia, and then dev . for FLOWExaFMM
# clone GXBeam, remove implicitAD compat, add main ImplicitAD, dev GXBeam