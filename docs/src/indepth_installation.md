
# In-Depth Installation Instructions for New Users, Developer Install, and Julia 1.10 and Older

The OWENS software has been developed and designed to highly flexible and to operate in the paradigm similar to modern open source software, leveraging tools such as the terminal, git, public software repositories, and automated package management both for the operating system and the programming language. Before attempting these instructions, if you are not familiar with these types of tools, please consider becoming familiar with them prior to proceeding.  Here are some of the first google hits for guides:

- https://www.redhat.com/sysadmin/beginners-guide-vim
- https://www.freecodecamp.org/news/the-beginners-guide-to-git-github/
- https://www.howtogeek.com/63997/how-to-install-programs-in-ubuntu-in-the-command-line/

The use of precompiled binaries for the optional OpenFAST fortran libraries significantly simplifies use and installation.  However, if you need to make modifications to the OpenFAST coupling, it requires an entirely different level of software development skills, and would require the user to become proficient at developing and compiling the standalone OpenFAST code first.

## Main New Programs' Installation (Julia, Paraview, Optionally VSCode)
Install julia, paraview, and visual studio by downloading/installing the OS specific versions for each, either manually (links below for your convenience) or through your favorite package manager:

- https://julialang.org/downloads/
- https://www.paraview.org/download/
- https://visualstudio.microsoft.com/downloads/

Be sure julia is on your path (should happen automatically as part of the installation process)

It is recommended, but not required to set up VS code as also described below.

Windows notes: 
- details are key, hash mismatch may mean that you are behind a proxy and the proxy isn't setup correctly.  For visual studio, the default installation and recommendations should work.  If importing the provided vs-code profile, the quick keys may need to be remapped to control from cmd.  
- Set up the environment variables below and be sure Julia is also on the path (https://www.wikihow.com/Create-an-Environment-Variable-in-Windows-10)


# Environment Variables
If you are using a proxy, such as those at the national labs be sure that the proxy variables are also declared/exported in your .bash_profile or .bashrc or the equivalent list of environment variables in Windows.  Additionally, you may need to tell your system where the SSL cert file is, and CA bundle. For those at Sandia, more information can be found at: https://wiki.sandia.gov/pages/viewpage.action?pageId=227381234#SandiaProxyConfiguration,Troubleshooting&HTTPS/SSLinterception-SSLCertificate.1.  

You cay need to then tell julia where the SSL_CA_ROOTS are.  Then in some instances, the built in julia git system has trouble with proxies, so you can tell it to use the command line interface git.

Here is an example section of a mac .zshrc script with dummy variables filled in assuming the machine is behind a proxy. If you aren't behind a proxy you should be able to just ignore all of the first block, and possibly the SSL CA, and CLI_GIT variables depending on your operating system and setup.

    export HTTP_PROXY="http://user:nopass@proxy.yourproxy.com:80/"
    export HTTPS_PROXY="http://user:nopass@proxy.yourproxy.com:80/"
    export no_proxy=*.local,169.254/16,*.yourproxy.com,*.srn.yourproxy.com,localhost,127.0.0.1,::1
    export NO_PROXY=*.local,169.254/16,*.yourproxy.com,*.srn.yourproxy.com,localhost,127.0.0.1,::1
    git config --global http.proxy http://user:nopass@proxy.yourproxy.com:80
    git config --global https.proxy http://user:nopass@proxy.yourproxy.com:80

    export SSL_CERT_FILE=/Users/kevmoor/.ssh/your_root_ca.pem
    export REQUESTS_CA_BUNDLE=${SSL_CERT_FILE}

    export JULIA_SSL_CA_ROOTS_PATH=""
    export JULIA_PKG_USE_CLI_GIT=true
    export JULIA_NUM_THREADS=1
    export OMP_NUM_THREADS=10

# Test That Julia Runs
the following should get you in and out of the julia interactive repl, open up a terminal and type

    julia 
    exit()

# Set up SSH Keys
Note that for installation behind the Sandia network, you will need to be on the network and follow additional instructions at https://wiki.sandia.gov/pages/viewpage.action?pageId=227381234#SandiaProxyConfiguration,Troubleshooting&HTTPS/SSLinterception-SSLCertificate.1. For windows, note that you may have to use id_ecdsa keys.

    # Go to your gihub account settings
    # left side, SSH and GPG keys
    # new ssh key
    # name: owensrepos # or whatever you'd like
    # back in the terminal, use the following, or the equivalent for Windows found in the github docs (https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
    ssh-keygen -t rsa -m PEM -C username@youremail.gov
    # enter, enter, enter (i.e. use defaults with !NO PASSWORD! since it is not necessary)
    cd ~
    ls -a
    cd .ssh
    vim id_rsa.pub
    #copy the contents to github.com (User icon > Settings > SSH and GPG > New SSH Key) and paste them back in your browser to the ssh key box and create the key
    # esc : q enter # to get out of vim
    cd ~

Additionally, if you find that your ssh is erroring when you try to install packages, try editing your ~/.ssh/config and add:

    Host *
    PubkeyAcceptedAlgorithms +ssh-rsa
    PubkeyAcceptedAlgorithms +ssh-ed25519

# Brief Julia Primer
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
```

## OWENS Installation for Julia 1.11+

```julia
using Pkg
Pkg.add(PackageSpec(url="git@github.com:sandialabs/OWENS.jl.git"))
```

Note that there are many packages used in the examples.  While they are installed within the OWENS.jl environment, if you want to additionally install them in your 1.11+ environment where you will likely be running from:
```julia
using Pkg
Pkg.add("Statistics");Pkg.add("Dierckx");Pkg.add("QuadGK");Pkg.add("FLOWMath");Pkg.add("HDF5");Pkg.add("ImplicitAD");Pkg.add("GXBeam");
Pkg.add("DelimitedFiles");Pkg.add("Statistics");Pkg.add("FFTW")
```

If you want to show any of the plots in the examples, they currently use the PyPlot interface, which means that julia has to install its own conda in the back end, which can take some time.  Alternatively, you can point to your own python if desired.
```julia
using Pkg
# Install PyPlot if not already installed
Pkg.add("PyPlot") #Note, this will take a while (maybe 10 min depending on your connection) since it is pulling conda and installing it behind the ~/.julia folder 
# if you want to use your already installed python, you can instead run
# ENV["PYTHON"] = "path to your desired python install"
# Pkg.add("PyCall")
# Pkg.add("PyPlot")
```

## OWENS Installation for Julia 1.10 and older (down to 1.6)
These steps require a secure download, such as through the SSH keys detailed above, to avoid man-in-the-middle attacks. 
```julia

using Pkg

println("\n#####################")
println("Install OWENS")
println("#####################")

Pkg.add("Statistics");Pkg.add("Dierckx");Pkg.add("QuadGK");Pkg.add("FLOWMath");Pkg.add("HDF5");Pkg.add("ImplicitAD");Pkg.add("GXBeam");
Pkg.add(PackageSpec(url="https://github.com/byuflowlab/Composites.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:sandialabs/OWENSPreComp.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:sandialabs/OWENSOpenFAST_jll.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:sandialabs/OWENSOpenFASTWrappers.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:sandialabs/OWENSAero.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:sandialabs/OWENSFEA.jl.git"))
Pkg.add(PackageSpec(url="git@github.com:sandialabs/OWENS.jl.git"))

```

# Testing Your Build of OWENS

clone the owens repository which contains example scripts that will setup and run example turbines from end to end

    git clone git@github.com:sandialabs/OWENS.jl

If you get an error about attempting to access library xyz, but it doesn't exist, check the path to openfast in the scripts at the top level of the error to make sure the path and library file matches, most of these are: 
```julia
adi_lib = "path/to/openfast/build/modules/libraryfolder/libraryname"
```

    cd OWENS.jl/examples/literate/
    julia B_detailedInputs.jl

You can visualize the output vtk/vtu/pvd paraview files with paraview, install paraview via
    apt-get -y install paraview # or on mac, brew install paraview

You can also run julia more interactively and maintain variables in scope for inspections etc if you don't have an IDE set up (but be careful of assuming a variable was cleared when it wasn't!) by starting the repl and essentially copying and pasing the run script via

    julia

```julia
    include("path/to/file.jl")
```

# Visual Studio Code IDE

You can install VScode and get a debugger etc.  In VScode, there is a setting file which sets up VS code for julia and sets some quick keys that can be changed if desired (OWENS.jl/docs/OWENS.code-profile). 

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

## VSCode Julia Debugger
It is a lot like matlab, but you are working with a compiled code, so in the debugger portion of the vscode gui, you need to check which code is compiled and which is interpereted, and turn off the compilation for the packages you are trying to debug.  The debugger will not step through compiled code. Also, some lines of code have many instructions, particularly function calls with conversions on the inputs, etc, so if you are trying to use the step in gui button, you may need to click it multiple times.

If you are working on a module and want it to reload with your most recent local changes without committing to master, pushing, and telling julia to update the package which is pointing to the git repo:

## Install custom repositories you want to develop

start Julia from the cloned directory and use the command:

    ] dev .
    
This type of installation will cause the module to reload each time Julia starts without needing to tell Julia to update. NOTE: even though julia 1.11 can to a regular install with a single line and resolve the custom dependencies, for the dev install it doesn't appear to automatically do so; so you may need to be careful about the installation order of dev packages to ensure you have downstream dependencies satisfied.

Alternatively, you can effectively copy and paste the module in, like so:
```julia
include("path/to/module.jl/source/module.jl")
```
then you don't even have to restart julia when you make changes, but be careful to only do this for a limited number of modules, and if you are changing constants, like c library interfaces, or the libraries themselves, then you need to restart julia to get it to pick up the most recent changes.

You can also install a specific branch of a remote repository package without having to clone the repo and checkout the branch:

```julia
using Pkg
Pkg.add(url = "git@github.com:sandialabs/OWENS.jl.git", rev = "dev")
```

# Custom Install of Optional OpenFAST Dependices
If your system is already set up such that it is capable of compiling OpenFAST, and you are on mac or linux, then you may skip this and rely on the automatically compiled version that are created when the OWENSOpenFAST libraries are installed by Julia.
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
    make moordyn_c_binding
    make hydrodyn_c_binding
    make aerodyn_inflow_c_binding
    make aerodyn_driver
    make turbsim
    cd ../../

There is also a OWENSOpenFASTWrappers.jl/deps/legacy_build.jl script which can be useful.

Then, within OWENS, you would provide the path to the resulting binaries you'd like to use (e.g. adi_lib="buildpath/openfast/build/modules/AeroDyn/libaerodyn_inflow_c_binding.so") instead of the default, which points to the precompiled OWENSOpenFAST_jll.jl binaries.