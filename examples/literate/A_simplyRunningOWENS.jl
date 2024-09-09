# # [Simply Running OWENS](@id simple1)
# 
# In this example, we show the first level of what is going on behind the precompiled binary
# Running julia directly with this as a starting point could make things like automating many runs in 
# a way that is not compatible with the current interface, but your design design fits.
#
# OWENS is comprised of many building blocks.  These series of examples progressively shows the internals
# of several of the key building blocks a new user might employ for their projects.  Fundamentally, OWENS has been 
# built to be as generalizable as possible. The lowest level of building blocks enable this, however, there are many
# common use cases for which helper functions have been developed, such as for meshing certain standard architectures
# and calculating and applying sectional properties to these architectures. The figure below summarizes this at a high 
# level.  
# TODO: yml file definition and inputs expanded 
#
# ![](../assets/OWENS_Example_Figure_Building_Blocks.png)
#
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook todo: get link working:
#-

import OWENS

runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]
#md runpath = path = splitdir(@__FILE__)[1]

Inp = OWENS.MasterInput("$runpath/sampleOWENS.yml")

OWENS.runOWENS(Inp,runpath)

# Here is an example of using the same model against the automated DLC run script.
# Note that for a setup cutom to a specific design, you'll want to go to the B level to get all of the detailed inputs correct
# One of these is the controller where a discon controller library can be coupled instead of the specified RPM control.

simulated_time = 2.0 #seconds
 DLCs = ["1_1"] #"normal" 
#md DLCs = ["1_3"] #"normal" 
#md DLCs = ["1_4"] #"normal" 
#md DLCs = ["1_5"] #"normal" 
#md DLCs = ["2_1"] #"freewheelatNormalOperatingRPM" 
#md DLCs = ["2_3"] #"freewheelatNormalOperatingRPM" 
#md DLCs = ["3_1"] #"startup" 
#md DLCs = ["3_2"] #"startup" 
#md DLCs = ["3_3"] #"startup" 
#md DLCs = ["4_1"] #"shutdown" 
#md DLCs = ["4_2"] #"shutdown" 
#md DLCs = ["5_1"] #"emergencyshutdown" 
#md DLCs = ["6_1"] #"parked" 
#md DLCs = ["6_2"] #"parked_idle" 
#md DLCs = ["6_4"] #"parked" 
#md DLCs = ["7_1"] #"parked" 
#md DLCs = ["2_3","3_1","3_2","3_3","4_1","4_2","5_1"]


OWENS.runDLC(DLCs,Inp,runpath;
    IEC_std="\"1-ED3\"",
    WindChar="\"A\"",
    WindClass=1,
    NumGrid_Z=38,
    NumGrid_Y=26,
    Vdesign=11.0,
    grid_oversize=1.25,
    Vinf_range=[10.0],#LinRange(4,24,21),
    regenWindFiles=true,
    delta_t_turbsim=0.05,
    simtime_turbsim=30.0,
    pathtoturbsim="$runpath/../../../OWENSOpenFASTWrappers.jl/deps/openfast/build/modules/turbsim/turbsim",
    runScript=OWENS.runOWENS)

nothing