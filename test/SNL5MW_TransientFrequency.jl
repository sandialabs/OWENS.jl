path = splitdir(@__FILE__)[1]

# import OWENS
include("$path/../src/OWENS.jl")

# import PyPlot

using Test
import HDF5

test_transient = true
test_modal = true
test_flutter = false
println("Starting")

## ****************** COORDINATE SYSTEM DEFINITION *********************** #
# 1 - x -  surge (pointed downstream)
# 2 - y -  sway  (right hand rule)
# 3 - z -  heave (pointed upwards)
# 4 - Ox - roll
# 5 - Oy - pitch
# 6 - Oz - yaw
# *********************************************************************** #

# use this benchmark file
bmOwens = "$path/data/_15mTower_transient_dvawt_c_2_lcdt"
# append this name to the end of the saved files
outFileExt = "_15mTowerExt_NOcentStiff"

# convert rotational stiffness to N-m/rad
StiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6].*0
convRotStiff = [1 1 1 180/pi 180/pi 180/pi]
platStiffDiag = StiffDiag_Nm_deg .* convRotStiff

# filename root to save the created nodal file
platfileRoot = "$path/data/1_FourColumnSemi_2ndPass"
platMassDiag = [9.8088e6 9.7811e6 1.8914e7 3.6351e9 3.6509e9 2.4362e9].*0
platStiffDiag_Nm_deg = [1.329e5 1.329e5 1.985e6 3.993e6 3.995e6 1.076e6].*0

# define the filename saving convention
fname = string(platfileRoot,outFileExt)

# *********************************************************************
# perform operations for the nodal file generation
# *********************************************************************
MassVal = platMassDiag
StiffVal = platStiffDiag

nodes = [1 1]
cmkType = ["M6" "K6"]
cmkValues = zeros(6,6,2)
for dd = 1:6
    # set up mass matrix
    cmkValues[dd,dd,1] = MassVal[dd]
    # set up stiffness matrix
    cmkValues[dd,dd,2] = StiffVal[dd]
end
OWENS.writeOwensNDL(fname, nodes, cmkType, cmkValues)

# *********************************************************************
# perform operations for the aerodynamic forces file generation
# *********************************************************************
CACTUSfileRoot = "$path/data/DVAWT_2B_LCDT"
OWENSfileRoot = bmOwens

operatingRPM = 7.2 # rpm
Nrpm = 10    # number of rpm stations
maxRPM = 10
Nmodes = 4  # number of modes to calculate/extract:
timeStep = 0.01
timeSim = 10.0       # [sec]
n_t = timeSim/timeStep # length of time vector
timeArray = [0, timeSim+1]
rpmArray  = [operatingRPM, operatingRPM]
omegaArrayHz = rpmArray ./ 60
omegaArrayHz2 = [7.1]/60

## ************************************************************************
# perform the transient simulations using OWENS
# *************************************************************************

println("Running Transient")

OWENS.owens(string(fname, ".owens"),"TNB";iterationType="DI", delta_t=timeStep, numTS=floor(timeSim/timeStep), nlOn=false, turbineStartup=0, usingRotorSpeedFunction=false, tocp=timeArray, Omegaocp=omegaArrayHz)

# Perform Tests
tol = 1e-5
n_t = 50

new_filename = "$path/data/new5MWnofileio2.h5"



t = HDF5.h5read(new_filename,"t")
aziHist = HDF5.h5read(new_filename,"aziHist")
OmegaHist = HDF5.h5read(new_filename,"OmegaHist")
OmegaDotHist = HDF5.h5read(new_filename,"OmegaDotHist")
gbHist = HDF5.h5read(new_filename,"gbHist")
gbDotHist = HDF5.h5read(new_filename,"gbDotHist")
gbDotDotHist = HDF5.h5read(new_filename,"gbDotDotHist")
FReactionHist = HDF5.h5read(new_filename,"FReactionHist")
rigidDof = HDF5.h5read(new_filename,"rigidDof")
genTorque = HDF5.h5read(new_filename,"genTorque")
genPower = HDF5.h5read(new_filename,"genPower")
torqueDriveShaft = HDF5.h5read(new_filename,"torqueDriveShaft")
uHist = HDF5.h5read(new_filename,"uHist")
eps_xx_0_hist = HDF5.h5read(new_filename,"eps_xx_0_hist")
eps_xx_z_hist = HDF5.h5read(new_filename,"eps_xx_z_hist")
eps_xx_y_hist = HDF5.h5read(new_filename,"eps_xx_y_hist")
gam_xz_0_hist = HDF5.h5read(new_filename,"gam_xz_0_hist")
gam_xz_y_hist = HDF5.h5read(new_filename,"gam_xz_y_hist")
gam_xy_0_hist = HDF5.h5read(new_filename,"gam_xy_0_hist")
gam_xy_z_hist = HDF5.h5read(new_filename,"gam_xy_z_hist")

import MAT
file = MAT.matopen("$path/../../archiveowens.jl/src/input_files/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.mat")
tmat = MAT.read(file,"t")
aziHistmat = MAT.read(file,"aziHist")
OmegaHistmat = MAT.read(file,"OmegaHist")
OmegaDotHistmat = MAT.read(file,"OmegaDotHist")
gbHistmat = MAT.read(file,"gbHist")
gbDotHistmat = MAT.read(file,"gbDotHist")
gbDotDotHistmat = MAT.read(file,"gbDotDotHist")
FReactionHistmat = MAT.read(file,"FReactionHist")
rigidDofmat = MAT.read(file,"rigidDof")
genTorquemat = MAT.read(file,"genTorque")
genPowermat = MAT.read(file,"genPower")
torqueDriveShaftmat = MAT.read(file,"torqueDriveShaft")
uHistmat = MAT.read(file,"uHist")
# eps_xx_0_histmat = MAT.read(file,"eps_xx_0_hist")
# eps_xx_z_histmat = MAT.read(file,"eps_xx_z_hist")
# eps_xx_y_histmat = MAT.read(file,"eps_xx_y_hist")
# gam_xz_0_histmat = MAT.read(file,"gam_xz_0_hist")
# gam_xz_y_histmat = MAT.read(file,"gam_xz_y_hist")
# gam_xy_0_histmat = MAT.read(file,"gam_xy_0_hist")
# gam_xy_z_histmat = MAT.read(file,"gam_xy_z_hist")
close(file)
ipt = 24*6
import PyPlot
PyPlot.figure()
PyPlot.plot(t,FReactionHist[:,6])
PyPlot.plot(tmat[:],FReactionHistmat[:,6])

PyPlot.figure()
PyPlot.plot(t,uHist[ipt,:],label="New")
PyPlot.plot(tmat[:],uHistmat[ipt,:],label="old")
PyPlot.legend()


import FFTW


mylabel = ["origMat","new"]
myt = [tmat,t]
for (i,signal) in enumerate([-uHistmat[ipt,:],-uHist[ipt,:]])
    L = length(signal)
    if L%2 != 0
        signal = signal[1:end-1]
        L = length(signal)
    end
    # signal .-= mean(signal)
    Y = FFTW.fft(signal)

    t = myt[i]
    Fs = 1/(t[2]-t[1])
    P2 = abs.(Y./L)
    P1 = P2[1:Int(L/2)+1]
    P1[2:end-1] = 2*P1[2:end-1]

    f = Fs.*(0:Int(L/2))./L
    PyPlot.figure(ipt)
    PyPlot.plot(f,P1,label=mylabel[i])
    PyPlot.title("Single-Sided Amplitude Spectrum of X(t)")
    PyPlot.xlabel("f (Hz)")
    PyPlot.ylabel("|P1(f)|")
    PyPlot.xlim([0.0,10.0])
    PyPlot.legend()
    # PyPlot.savefig("$(path)/../figs/34m_fig5_6_upperRootFlatwiseStressSpectrumSteadyOperation.pdf",transparent = true)
end
