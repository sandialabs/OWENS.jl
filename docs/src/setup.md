
## Windows

At this stage in the software's maturity, it is recommended to use mac or linux environments unless the user is experienced with compiled software development in a windows environment. We are also currently testing the Linux Shell for Windows, which should be the same process as the Linux instructions below, but has not been tested yet.

Install julia, paraview, and visual studio manually by downloading the windows executables for

- https://julialang.org/downloads/
- https://www.paraview.org/download/
- https://visualstudio.microsoft.com/downloads/

Also download OpenFAST:
- https://github.com/OpenFAST/openfast

Be sure julia is on your path, and follow the windows compilation instructions for the openfast Inflowwind, AeroDyn, MoorDyn and HydroDyn libraries. Installation is otherwise the same as the Linux instructions below

## Mac

Same installation as Linux except we recommend using the homebrew package manager, so exchange all "apt-get" with "brew" 

## Linux

This was tested using a M1 mac using docker so I had a fresh build, this will take 1-2 hours assuming no issues. If you do decide to use Docker, be sure to specify the correct platform so you don't loose performance via emulation.  This is what I used on my M1 mac:

docker run -i -t --platform linux/amd64 --name myUbuntuLinux2 ubuntu:latest

you should use whatever linux environment/workstation is most comfortable.  Note that the Docker instance can save plots, but I don't believe it can show them, at least not my setup.

# Required Compilers
    apt-get update -y
    apt-get install git -y
    apt-get install wget -y
    apt-get install vim -y
    apt-get install cmake -y
    apt-get install gfortran -y
    apt-get install build-essential -y
    apt-get install libblas-dev liblapack-dev -y

# Install julia
    cd ~
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.3-linux-x86_64.tar.gz
    tar zxvf julia-1.9.3-linux-x86_64.tar.gz
    rm julia-1.9.3-linux-x86_64.tar.gz
    export PATH="$PATH:/root/julia-1.9.3/bin"

in your ~/.bashrc file, tell julia to use the command line git by inserting the following:

export JULIA_PKG_USE_CLI_GIT=true

# Proxy Setup
If you are using a proxy, be sure that the proxy variables are declared/exported in your .bash_profile or .bashrc or the equivalent

    http_proxy, https_proxy, HTTP_PROXY, HTTPS_PROXY, no_proxy, NO_PROXY
    git config --global http.proxy http://user:nopass@proxy.yourorg:number
    git config --global https.proxy http://user:nopass@proxy.yourorg:number
    export JULIA_SSL_CA_ROOTS_PATH=""
    export JULIA_SSL_NO_VERIFY_HOSTS="*.yourorgurl"
    export JULIA_PKG_USE_CLI_GIT=true 	

# Test That Julia Runs
the following should get you in and out of the julia interactive repl

    julia 
    exit()

# Set up SSH Keys

    # Make ssh keys and put in the correct places
    # Go to your gihub account settings
    # left side, SSH and GPG keys
    # new ssh key
    # name: owensrepos # or whatever you'd like
    # back in the linux terminal
    ssh-keygen -t rsa -m PEM -C username@youremail.gov
    # enter, enter, enter (i.e. use defaults)
    cd ~
    ls -a
    cd .ssh
    vim id_rsa.pub
    #copy the contents and paste them back in your browser to the ssh key box and create the key
    # esc : q enter # to get out of vim
    cd ~

# Install Optional OpenFAST Dependices
    mkdir coderepos
    cd coderepos
    # Install openfast coupled libraries !NOTE!: if you change the location of the compiled libraries, you may need to update the rpath variable, or recompile.
    git clone --depth 1 git@github.com:andrew-platt/openfast.git
    # if this errors, you can clone git@github.com:OpenFAST/openfast.git it just doesn't have the latest updates from Andy, but the interface should be the same and should run.
    cd openfast
    git remote set-branches origin '*'
    git fetch --depth 1 origin f/ADI_c_binding_multiRotor
    git checkout f/ADI_c_binding_multiRotor
    mkdir build
    cd build
    # can also add -DOPENMP=ON if desired for acceleration of OLAF
    # you can rebuild later by removing the build folder and following these instructions again.
    cmake -DBUILD_SHARED_LIBS=ON ..
    make ifw_c_binding
    # make moordyn_c_binding
    # make hydrodyn_c_binding
    make aerodyn_inflow_c_binding
    make aerodyn_driver
    make turbsim
    cd ../../

# Brief Julia Primer and OWENS Installation
Now open the julia interactive repl and run the following blocks, obviously a multi-line block should be entered as one.

julia

```julia

println("This is intended just to get you rolling, comprehensive documentation can be found at https://docs.julialang.org/en/v1/")
println("More detail on major differences between codes can be found at https://docs.julialang.org/en/v1/manual/noteworthy-differences/")

###############################################################
###############################################################
###############################################################

println("#####################")
println("Basic Data Handling")
println("#####################")
# Create a multidimensional array
newMatrix = zeros(2,3)
# Mutate the contents
newMatrix[1,2] = 1.0

#Create New Scalar
newScalar = 5.0

println("Print the first matrix")
# Print the contents and observe the mutation
for irow = 1:length(newMatrix[:,1])
    println(newMatrix[irow,:])
end

# Copy the newMatrix and scalar
newMatrix2 = newMatrix
newScalar2 = newScalar
# Mutate the second matrix and scalar
newMatrix[1,2] = 2.0
newScalar2 = 2.0

#Make a printing function inline
function printme(matrix, scalar)
    for irow = 1:length(matrix[:,1])
        println(matrix[irow,:])
    end
    println("Scalar: $scalar")
end

println("Printing the second matrix and scalar")
printme(newMatrix2,newScalar2)
println("Printing the first matrix again")
printme(newMatrix,newScalar)
if newMatrix[1,2] == newMatrix2[1,2]
    println("B=A references the arrays")
    println("B=copy(A) breaks the reference and does a true copy")
    println("B=deepcopy(A) is needed if it is a multi-level type, like a struct or dictionary")
    println("However, scalars are not linked.  This is because a scalar is directly looking at a memory element, while arrays are pointing to the memory elements")
end

###############################################################
###############################################################
###############################################################

println("\n#####################")
println("Scope of Functions")
println("#####################")


function coolfunction(input1,input2; mykeyname="default",mykeyname2=5.0)
    if mykeyname=="default"
        return input1.+input2[1,1,1], 1.0
    else
        input2[1,1,1] = 1.0 # Since arrays are always passed by reference, we can mutate it here and it will be mutated above
        return input1.+mykeyname2, 0.0
    end
end

outputs = coolfunction(ones(3).*5,zeros(4,5,2)) 
println("use the defaults for the optional args and dump the output into a tuple")
println("First output $(outputs[1]), Second output $(outputs[2])")

println("supply the optional args and dump the output into newly allocated items")
myinput = zeros(4,5,2)
mykeyname = "notdefault"
testinput = 5.0
outputs1, output2 = coolfunction(ones(3).*5,myinput;mykeyname,mykeyname2=testinput)
println("First output $(outputs1), Second output $(output2)")

println("now show that myinput was mutated within the function since it was passed by reference")
println(myinput)


###############################################################
###############################################################
###############################################################

println("\n#####################")
println("Types")
println("#####################")

first = 1.0
println(typeof(first))

second = 2
println(typeof(second))

###############################################################
###############################################################
###############################################################

println("\n#####################")
println("Structs")
println("#####################")

mutable struct mystruct
    coolterm1
    othercoolterm
end

newStruct = mystruct(1.0,2.0)

println(newStruct.coolterm1)
println(newStruct.othercoolterm)

###############################################################
###############################################################
###############################################################

println("\n#####################")
println("Package Manager with Standard Packages")
println("#####################")

using Pkg
Pkg.add("PyPlot") #Note, this will take a while (maybe 10 min depending on your connection) since it is pulling conda and installing it behind the ~/.julia folder 
# if you want to use your already installed python, you can instead run
# ENV["PYTHON"] = "path to your desired python install"
# Pkg.add("PyCall")
# Pkg.add("PyPlot")

Pkg.add("DelimitedFiles")

import PyPlot
import DelimitedFiles

x = [1; 2; 3; 4];

y = [5; 6; 7; 8];

open("delim_file.txt", "w") do io
        DelimitedFiles.writedlm(io, [x y])
       end

data = DelimitedFiles.readdlm("delim_file.txt", '\t', Int, '\n')

PyPlot.figure()
PyPlot.plot(data[:,1],data[:,2],label="data")
PyPlot.xlabel("x")
PyPlot.ylabel("y")
PyPlot.legend()

thisFilesPath = splitdir(@__FILE__)[1]

PyPlot.savefig("$(thisFilesPath)/saveme.pdf",transparent = true)

run(`rm $(thisFilesPath)/saveme.pdf`) #system run
rm("delim_file.txt") # julia's function that does the same thing

###############################################################
###############################################################
###############################################################

println("\n#####################")
println("Package Manager with Custom Packages")
println("#####################")

Pkg.add("Statistics");Pkg.add("Dierckx");Pkg.add("QuadGK");Pkg.add("FLOWMath");Pkg.add("HDF5");Pkg.add("ImplicitAD");
Pkg.add(PackageSpec(url="git@github.com:kevmoor/GXBeam.jl.git"))
Pkg.add(PackageSpec(url="https://github.com/byuflowlab/Composites.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/PreComp.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/OpenFASTWrappers.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/VAWTAero.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/GyricFEA.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:SNL-WaterPower/OWENS.jl.git"))

```

# Testing Your Build of OWENS

clone the owens repository which contains example run scripts, the turbine mesh generator

    git clone git@github.com:SNL-WaterPower/OWENS.jl

change the path to openfast in the scripts, most of these are: 
```julia
adi_lib = "path/to/openfast/build/modules/libraryfolder/libraryname"
```

    cd OWENS.jl/test
    julia Fig5_4_Fig5_3_normaloperation_torque_flatwise_generatorControl_Land.jl

You can visualize the output vtk/vtu/pvd paraview files with paraview, install paraview via
    apt-get -y install paraview

or you can run it more interactively by starting the repl

    julia

```julia
    include("path/to/file.jl")
```

or you can install VScode and get a debugger etc.  In VScode, there is a setting file which sets up VS code for julia and sets some quick keys that can be changed if desired (OWENS.jl/docs/setup/OWENS_Julia_VS.code-profile). 

With the sample profile loaded into VSCode, If you want to clear out all the variables and restart do 

    cmd-j cmd-k, 
    
if you want to clear out the console 

    cmd-j cmd-c

open the workspace 

    cmd-j cmd-w
    
run highlighted code 
    
    shift-enter

run the currently selected file
    
    cmd-shift-enter
    
You can also use the gui buttons.

**For debugging**, it is a lot like matlab, but you are working with a compiled code, so in the debugger portion of the vscode gui, you need to check which code is compiled and which is interpereted, and turn off the compilation for the packages you are trying to debug.  The debugger will not step through compiled code. Also, some lines of code have many instructions, particularly function calls with conversions on the inputs, etc, so if you are trying to use the step in gui button, you may need to click it multiple times.

If you are working on a module and want it to reload with your most recent local changes without committing to master, pushing, and telling julia to update the package which is pointing to the git repo:

# Install custom repositories you want to develop

start Julia from the cloned directory and use the command:

    ] dev .
    
This type of installation will cause the module to reload each time Julia starts without needing to tell Julia to update. You are developing the current directory

# alternatively, instead of using or import to get access to the module, within julia

```julia
include("path/to/module.jl/source/module.jl")
```
then you don't even have to restart julia when you make changes, but be careful to only do this for a limited number of modules, and if you are changing constants, like c library interfaces, or the libraries themselves, then you need to restart julia to get it to pick up the most recent changes.