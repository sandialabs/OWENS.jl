println("Starting")
println("Setting Up Conda Environment for Hydrodynamics")
import Conda
# (NOTE: If you are installing Python packages for use with PyCall, you must use the root environment.)
# And we are, so we use the root environment to build in
Conda.runconda(`config --set ssl_verify False`) #TODO: FIX THIS!!!!
# Conda.runconda(`env create -f environment.yml`)
Conda.add_channel("conda-forge")
# Conda.add("python==3.8")
Conda.add("capytaine==1.2")
Conda.add("numpy==1.19.2")
Conda.pip_interop(true)
Conda.pip("install", "pygmsh==6.1.0")
Conda.pip("install", "meshmagick @ https://github.com/LHEEA/meshmagick/archive/2.0.tar.gz")
# Conda.pip("install","pymap")
# # Resolve mkl error
Conda.add("nomkl")
try
    Conda.rm("mkl")
catch
end
Conda.update()
#
# println("Linking PyCall to the Conda Environment")
# using Pkg
# Pkg.add("PyCall")
# ENV["PYTHON"] = "$(Conda.ROOTENV)"
# Pkg.build("PyCall")
