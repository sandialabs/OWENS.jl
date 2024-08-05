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
using HDF5
using Test
using YAML
using OrderedCollections

runpath = splitdir(@__FILE__)[1]

OWENS_Options = OWENS.MasterInput("$runpath/modeling_options_OWENS_windioExample.yml")

WINDIO_filename = "$runpath/WINDIO_example.yaml"

windio = YAML.load_file(WINDIO_filename; dicttype=OrderedCollections.OrderedDict{Symbol,Any})

numadIn_bld_old = OWENS.readNuMadGeomCSV("$(runpath)$(OWENS_Options.NuMad_geom_xlscsv_file_bld)";section=:blade)
numadIn_bld_new = OWENS.readNuMadGeomCSV(windio;section=:blade)

for key in fieldnames(typeof(numadIn_bld_old))
    old_data = getfield(numadIn_bld_old,key)
    new_data = getfield(numadIn_bld_new,key)

    if old_data!=new_data
        println(key)
        println(old_data)
        println(new_data)
    end
    # println(isapprox(old_data,new_data))
end

OWENS.runOWENSWINDIO(windio,OWENS_Options,runpath)

# Alternatively OWENS.runOWENSWINDIO(WINDIO_filename,OWENS_Options,runpath)

file = "$runpath/InitialDataOutputs_UNIT.h5"
t_UNIT = HDF5.h5read(file,"t")
aziHist_UNIT = HDF5.h5read(file,"aziHist")
OmegaHist_UNIT = HDF5.h5read(file,"OmegaHist")
OmegaDotHist_UNIT = HDF5.h5read(file,"OmegaDotHist")
gbHist_UNIT = HDF5.h5read(file,"gbHist")
gbDotHist_UNIT = HDF5.h5read(file,"gbDotHist")
gbDotDotHist_UNIT = HDF5.h5read(file,"gbDotDotHist")
FReactionHist_UNIT = HDF5.h5read(file,"FReactionHist")
FTwrBsHist_UNIT = HDF5.h5read(file,"FTwrBsHist")
genTorque_UNIT = HDF5.h5read(file,"genTorque")
genPower_UNIT = HDF5.h5read(file,"genPower")
torqueDriveShaft_UNIT = HDF5.h5read(file,"torqueDriveShaft")
uHist_UNIT = HDF5.h5read(file,"uHist")
uHist_prp_UNIT = HDF5.h5read(file,"uHist_prp")
epsilon_x_hist_UNIT = HDF5.h5read(file,"epsilon_x_hist")
epsilon_y_hist_UNIT = HDF5.h5read(file,"epsilon_y_hist") 
epsilon_z_hist_UNIT = HDF5.h5read(file,"epsilon_z_hist")
kappa_x_hist_UNIT = HDF5.h5read(file,"kappa_x_hist")
kappa_y_hist_UNIT = HDF5.h5read(file,"kappa_y_hist")
kappa_z_hist_UNIT = HDF5.h5read(file,"kappa_z_hist") 
massOwens_UNIT = HDF5.h5read(file,"massOwens")
stress_U_UNIT = HDF5.h5read(file,"stress_U")
SF_ult_U_UNIT = HDF5.h5read(file,"SF_ult_U")
SF_buck_U_UNIT = HDF5.h5read(file,"SF_buck_U")
stress_L_UNIT = HDF5.h5read(file,"stress_L")
SF_ult_L_UNIT = HDF5.h5read(file,"SF_ult_L")
SF_buck_L_UNIT = HDF5.h5read(file,"SF_buck_L")
stress_TU_UNIT = HDF5.h5read(file,"stress_TU")
SF_ult_TU_UNIT = HDF5.h5read(file,"SF_ult_TU")
SF_buck_TU_UNIT = HDF5.h5read(file,"SF_buck_TU")
stress_TL_UNIT = HDF5.h5read(file,"stress_TL")
SF_ult_TL_UNIT = HDF5.h5read(file,"SF_ult_TL")
SF_buck_TL_UNIT = HDF5.h5read(file,"SF_buck_TL")
topstrainout_blade_U_UNIT = HDF5.h5read(file,"topstrainout_blade_U")
topstrainout_blade_L_UNIT = HDF5.h5read(file,"topstrainout_blade_L")
topstrainout_tower_U_UNIT = HDF5.h5read(file,"topstrainout_tower_U")
topstrainout_tower_L_UNIT = HDF5.h5read(file,"topstrainout_tower_L")
topDamage_blade_U_UNIT = HDF5.h5read(file,"topDamage_blade_U")
topDamage_blade_L_UNIT = HDF5.h5read(file,"topDamage_blade_L")
topDamage_tower_U_UNIT = HDF5.h5read(file,"topDamage_tower_U")
topDamage_tower_L_UNIT = HDF5.h5read(file,"topDamage_tower_L")


file = "$runpath/InitialDataOutputs.h5"
t = HDF5.h5read(file,"t")
aziHist = HDF5.h5read(file,"aziHist")
OmegaHist = HDF5.h5read(file,"OmegaHist")
OmegaDotHist = HDF5.h5read(file,"OmegaDotHist")
gbHist = HDF5.h5read(file,"gbHist")
gbDotHist = HDF5.h5read(file,"gbDotHist")
gbDotDotHist = HDF5.h5read(file,"gbDotDotHist")
FReactionHist = HDF5.h5read(file,"FReactionHist")
FTwrBsHist = HDF5.h5read(file,"FTwrBsHist")
genTorque = HDF5.h5read(file,"genTorque")
genPower = HDF5.h5read(file,"genPower")
torqueDriveShaft = HDF5.h5read(file,"torqueDriveShaft")
uHist = HDF5.h5read(file,"uHist")
uHist_prp = HDF5.h5read(file,"uHist_prp")
epsilon_x_hist = HDF5.h5read(file,"epsilon_x_hist")
epsilon_y_hist = HDF5.h5read(file,"epsilon_y_hist")  
epsilon_z_hist = HDF5.h5read(file,"epsilon_z_hist")
kappa_x_hist = HDF5.h5read(file,"kappa_x_hist")
kappa_y_hist = HDF5.h5read(file,"kappa_y_hist")
kappa_z_hist = HDF5.h5read(file,"kappa_z_hist") 
massOwens = HDF5.h5read(file,"massOwens")
stress_U = HDF5.h5read(file,"stress_U")
SF_ult_U = HDF5.h5read(file,"SF_ult_U")
SF_buck_U = HDF5.h5read(file,"SF_buck_U")
stress_L = HDF5.h5read(file,"stress_L")
SF_ult_L = HDF5.h5read(file,"SF_ult_L")
SF_buck_L = HDF5.h5read(file,"SF_buck_L")
stress_TU = HDF5.h5read(file,"stress_TU")
SF_ult_TU = HDF5.h5read(file,"SF_ult_TU")
SF_buck_TU = HDF5.h5read(file,"SF_buck_TU")
stress_TL = HDF5.h5read(file,"stress_TL")
SF_ult_TL = HDF5.h5read(file,"SF_ult_TL")
SF_buck_TL = HDF5.h5read(file,"SF_buck_TL")
topstrainout_blade_U = HDF5.h5read(file,"topstrainout_blade_U")
topstrainout_blade_L = HDF5.h5read(file,"topstrainout_blade_L")
topstrainout_tower_U = HDF5.h5read(file,"topstrainout_tower_U")
topstrainout_tower_L = HDF5.h5read(file,"topstrainout_tower_L")
topDamage_blade_U = HDF5.h5read(file,"topDamage_blade_U")
topDamage_blade_L = HDF5.h5read(file,"topDamage_blade_L")
topDamage_tower_U = HDF5.h5read(file,"topDamage_tower_U")
topDamage_tower_L = HDF5.h5read(file,"topDamage_tower_L")

frac = 1e-5
@test isapprox(t_UNIT,t;atol=maximum(abs.(t_UNIT))*frac)
@test isapprox(aziHist_UNIT,aziHist;atol=maximum(abs.(aziHist_UNIT))*frac)
@test isapprox(OmegaHist_UNIT,OmegaHist;atol=maximum(abs.(OmegaHist_UNIT))*frac)
@test isapprox(OmegaDotHist_UNIT,OmegaDotHist;atol=maximum(abs.(OmegaDotHist_UNIT))*frac)
@test isapprox(gbHist_UNIT,gbHist;atol=maximum(abs.(gbHist_UNIT))*frac)
@test isapprox(gbDotHist_UNIT,gbDotHist;atol=maximum(abs.(gbDotHist_UNIT))*frac)
@test isapprox(gbDotDotHist_UNIT,gbDotDotHist;atol=maximum(abs.(gbDotDotHist_UNIT))*frac)
@test isapprox(FReactionHist_UNIT,FReactionHist;atol=maximum(abs.(FReactionHist_UNIT))*frac)
@test isapprox(FTwrBsHist_UNIT,FTwrBsHist;atol=maximum(abs.(FTwrBsHist_UNIT))*frac)
@test isapprox(genTorque_UNIT,genTorque;atol=maximum(abs.(genTorque_UNIT))*frac)
@test isapprox(genPower_UNIT,genPower;atol=maximum(abs.(genPower_UNIT))*frac)
@test isapprox(torqueDriveShaft_UNIT,torqueDriveShaft;atol=maximum(abs.(torqueDriveShaft_UNIT))*frac)
@test isapprox(uHist_UNIT,uHist;atol=maximum(abs.(uHist_UNIT))*frac)
@test isapprox(uHist_prp_UNIT,uHist_prp;atol=maximum(abs.(uHist_prp_UNIT))*frac)
@test isapprox(epsilon_x_hist_UNIT,epsilon_x_hist;atol=maximum(abs.(epsilon_x_hist_UNIT))*frac)
@test isapprox(epsilon_y_hist_UNIT,epsilon_y_hist;atol=maximum(abs.(epsilon_y_hist_UNIT))*frac)
@test isapprox(epsilon_z_hist_UNIT,epsilon_z_hist;atol=maximum(abs.(epsilon_z_hist_UNIT))*frac)
@test isapprox(kappa_x_hist_UNIT,kappa_x_hist;atol=maximum(abs.(kappa_x_hist_UNIT))*frac)
@test isapprox(kappa_y_hist_UNIT,kappa_y_hist;atol=maximum(abs.(kappa_y_hist_UNIT))*frac)
@test isapprox(kappa_z_hist_UNIT,kappa_z_hist;atol=maximum(abs.(kappa_z_hist_UNIT))*frac)
@test isapprox(massOwens_UNIT,massOwens;atol=maximum(abs.(massOwens_UNIT))*frac)
ipass = 0
for i = 1:length(stress_U_UNIT)
    # println("$i of $(length(stress_U_UNIT))")
    if isapprox(stress_U_UNIT[i],stress_U[i];atol=maximum(abs.(stress_U_UNIT[i]))*frac)
        ipass += 1
    end
end
println("Percent Stress Pass: $(ipass/length(stress_U_UNIT)*100)%")
# @test isapprox(SF_ult_U_UNIT,SF_ult_U;atol=maximum(abs.(SF_ult_U_UNIT))*frac)
# @test isapprox(SF_buck_U_UNIT,SF_buck_U;atol=maximum(abs.(SF_buck_U_UNIT))*frac)
# @test isapprox(stress_L_UNIT,stress_L;atol=maximum(abs.(stress_L_UNIT))*frac)
# @test isapprox(SF_ult_L_UNIT,SF_ult_L;atol=maximum(abs.(SF_ult_L_UNIT))*frac)
# @test isapprox(SF_buck_L_UNIT,SF_buck_L;atol=maximum(abs.(SF_buck_L_UNIT))*frac)
# @test isapprox(stress_TU_UNIT,stress_TU;atol=maximum(abs.(stress_TU_UNIT))*frac)
# @test isapprox(SF_ult_TU_UNIT,SF_ult_TU;atol=maximum(abs.(SF_ult_TU_UNIT))*frac)
# @test isapprox(SF_buck_TU_UNIT,SF_buck_TU;atol=maximum(abs.(SF_buck_TU_UNIT))*frac)
# @test isapprox(stress_TL_UNIT,stress_TL;atol=maximum(abs.(stress_TL_UNIT))*frac)
# @test isapprox(SF_ult_TL_UNIT,SF_ult_TL;atol=maximum(abs.(SF_ult_TL_UNIT))*frac)
# @test isapprox(SF_buck_TL_UNIT,SF_buck_TL;atol=maximum(abs.(SF_buck_TL_UNIT))*frac)
# @test isapprox(topstrainout_blade_U_UNIT,topstrainout_blade_U;atol=maximum(abs.(topstrainout_blade_U_UNIT))*frac)
# @test isapprox(topstrainout_blade_L_UNIT,topstrainout_blade_L;atol=maximum(abs.(topstrainout_blade_L_UNIT))*frac)
# @test isapprox(topstrainout_tower_U_UNIT,topstrainout_tower_U;atol=maximum(abs.(topstrainout_tower_U_UNIT))*frac)
# @test isapprox(topstrainout_tower_L_UNIT,topstrainout_tower_L;atol=maximum(abs.(topstrainout_tower_L_UNIT))*frac)
# @test isapprox(topDamage_blade_U_UNIT,topDamage_blade_U;atol=maximum(abs.(topDamage_blade_U_UNIT))*frac)
# @test isapprox(topDamage_blade_L_UNIT,topDamage_blade_L;atol=maximum(abs.(topDamage_blade_L_UNIT))*frac)
# @test isapprox(topDamage_tower_U_UNIT,topDamage_tower_U;atol=maximum(abs.(topDamage_tower_U_UNIT))*frac)
# @test isapprox(topDamage_tower_L_UNIT,topDamage_tower_L;atol=maximum(abs.(topDamage_tower_L_UNIT))*frac)