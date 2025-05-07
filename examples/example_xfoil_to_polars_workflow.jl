using Xfoil, Printf
import PyPlot
import DelimitedFiles
import Statistics:mean
import CCBlade
import FLOWMath
PyPlot.pygui(true)
PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=3.0)
PyPlot.rc("legend", frameon=false)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

path,_ = splitdir(@__FILE__)

# function airfoil(
yscale=[1.0]
i_station=1
filename = "$path/airfoils/NACA_0021_DMS.csv"
af_xy = DelimitedFiles.readdlm(filename,',',Float64,skipstart = 0)
# af_xy = DelimitedFiles.readdlm("$path/airfoils/naca2412.dat",Float64,skipstart = 0)
x_xf_1 = af_xy[:,1]
y_xf_1 = af_xy[:,2]

# Unpack
xnode = af_xy[:,1]
ynode = af_xy[:,2]

# Filter out repeated max and min indices
max_x = maximum(xnode)
min_x = minimum(xnode)
count_max = 0
count_min = 0
unique_idx = []
for i_x = 1:size(xnode)[1]
    # only push the max value once
    if xnode[i_x]==max_x && count_max==0
        global count_max += 1
        push!(unique_idx,i_x)
    end
    # only push the min value once
    if xnode[i_x]==min_x && count_min==0
        global count_min += 1
        push!(unique_idx,i_x)
    end
    # otherwise, if not max or min, push
    if xnode[i_x]!=max_x && xnode[i_x]!=min_x
        push!(unique_idx,i_x)
    end
end

xnode = xnode[unique_idx]
ynode = ynode[unique_idx]

# find leading edge
le_idx = argmin(xnode)
# shift so leading edge is first, in case the leading edge was defined slightly differently or mid-array
xnode = circshift(xnode,-(le_idx-1))
ynode = circshift(ynode,-(le_idx-1))

# Move the leading edge to 0 if not already there
xnode = xnode .- xnode[1]

# Check that the airfoil leading edge is at zero
le_idx = argmin(xnode)
if !isapprox(ynode[le_idx],0.0;atol=1e-3)
    @error "Airfoil leading edge not at y=0 for (numadIn.airfoil[i_station]).csv, the chord line must be at 0,0 (leading edge) to 1,0 (trailing edge) per standard airfoil definition."
end

# Check that the airfoil trailing edge is closed
te_idx = argmax(xnode)
if !isapprox(ynode[te_idx],0.0;atol=1e-3)
    @warn "Airfoil trailing edge not closed for (numadIn.airfoil[i_station]).csv, moving trailing edge airfoil point to the normalized 1,0 point, trailing edges should be closed, which can be accomplished with a flatback by inserting a point numerically beyond the trailing edge with a value of 0"
    # insert!(xnode,te_idx,maximum(xnode)+1e-5)
    # insert!(ynode,te_idx,0.0)
    ynode[te_idx] = 0.0
end

# Add a redundant leading edge at the end, we'll remove it if it doesn't end up as the leading edge below
push!(xnode,xnode[1])
push!(ynode,ynode[1])

# The data should go around the top/suction side first for precomp, so reverse the array if that isn't the case
if mean(ynode[1:te_idx])<mean(ynode[te_idx+1:end])
    reverse!(xnode)
    reverse!(ynode)
end


# now, as mentioned above, we need to remove the redundant leading edge
xnode = xnode[1:end-1]
ynode = ynode[1:end-1]

# normalize by the chord length
max_x = maximum(xnode)
x = xnode./max_x
y = ynode./max_x.*yscale[i_station] #optionally multiply by a thickness factor for design

# plot the airfoil geometry
PyPlot.figure("PreComp")
PyPlot.scatter(x, y,label="orig")#, label="", framestyle=:none, aspect_ratio=1.0, show=true)

# For XFoil, it goes from TE to LE to TE positive to negative, so we should just have to circleshift backwards, and then check the airfoil is closed with repeated trailing edge positions
# find trailing edge
te_idx = argmax(x)
# shift so leading edge is first, in case the leading edge was defined slightly differently or mid-array
x_xf = circshift(x,(te_idx-1))
y_xf = circshift(y,(te_idx-1))
# Add in the redundant trailing edge if it isn't there already
if x_xf[end] != 1.0
    push!(x_xf,x_xf[1])
    push!(y_xf,y_xf[1])
end

# The data should go around the top/suction side first for xfoil, so reverse the array if that isn't the case
if mean(y_xf[1:te_idx])<mean(y_xf[te_idx+1:end])
    reverse!(x_xf)
    reverse!(y_xf)
end

# Check if the airfoil trailing edge is a flat back, if it is, then remove those points and rescale
if isapprox(x_xf[1],x_xf[2];atol=1e-3)
    @warn "It looks like you have a flat back airfoil. To be compatible with xfoil, we're remving the trailing edge zero so xfoil's pane operation works correctly"
    x_xf = x_xf[2:end-1]
    x_xf = x_xf./maximum(x_xf)
    y_xf = y_xf[2:end-1]
end

# plot the airfoil geometry
PyPlot.figure("XFoil")
PyPlot.plot(x_xf, y_xf,".-",label="orig")#, label="", framestyle=:none, aspect_ratio=1.0, show=true)
PyPlot.plot(x_xf_1, y_xf_1,".-",label="other")#, label="", framestyle=:none, aspect_ratio=1.0, show=true)
PyPlot.legend()
# load airfoil coordinates into XFOIL
Xfoil.set_coordinates(x_xf,y_xf)



# # repanel using XFOIL's `PANE` command
xr, yr = Xfoil.pane()

# xr = x_xf
# yr = y_xf

# plot the refined airfoil geometry
PyPlot.plot(xr, yr,".",label="pane")#, label="", framestyle=:none, aspect_ratio=1.0, show=true)
PyPlot.legend()

Rerange = round.(unique(sort([logrange(1e4, 1e7, length=4);LinRange(1e5,3e5,10)]));digits=-3)

# Write the new file for OWENS
new_filename = "$(filename[1:end-3])dat"
open(new_filename, "w") do file
    # Write new data to file
        write(file, "Title: $new_filename
Thickness to Chord Ratio: $(maximum(yr)-minimum(yr))
Zero Lift AOA (deg): 0
Reverse Camber Direction: 0\n")

    for Re in Rerange
        # set operating conditions
        alpha = -20:1:20
        # re = 1e5

        c_l_raw, c_d_raw, c_dp_raw, c_m_raw, converged = Xfoil.alpha_sweep(xr, yr, alpha, Re, iter=100, zeroinit=true, printdata=true, reinit=false)



        c_l = FLOWMath.linear(alpha[converged.==1],c_l_raw[converged.==1],alpha)
        c_d = FLOWMath.linear(alpha[converged.==1],c_d_raw[converged.==1],alpha)
        c_dp = FLOWMath.linear(alpha[converged.==1],c_dp_raw[converged.==1],alpha)
        c_m = FLOWMath.linear(alpha[converged.==1],c_m_raw[converged.==1],alpha)
        
        maxcl,idxmax = findmax(c_l)
        aoamaxcl = alpha[idxmax]
        
        # Viterna extrapolation
        cr75 = 0.1
        alpha_ext, cl_ext, cd_ext = CCBlade.viterna(alpha ./ 360.0 .* 2π, c_l, c_d, cr75)
        alpha_ext = alpha_ext ./ (2π) .* 360
        
        # Sort by alpha
        sort_inds = sortperm(alpha_ext)
        alpha_ext = alpha_ext[sort_inds]
        cl_ext = cl_ext[sort_inds]
        cd_ext = cd_ext[sort_inds]
        
        PyPlot.figure("Re: $Re")
        PyPlot.plot(alpha,c_l)
        PyPlot.plot(alpha,c_d)
        PyPlot.plot(alpha,c_l_raw,label="raw")
        PyPlot.plot(alpha,c_d_raw,label="raw")
        PyPlot.plot(alpha_ext,cl_ext,".",label="cl_vit")
        PyPlot.plot(alpha_ext,cd_ext,".",label="cd_vit")
        PyPlot.plot(aoamaxcl,maxcl,".",label="max")
        PyPlot.ylim([-2.0,2.0])
        PyPlot.legend()


        write(file, "\nReynolds Number: $Re
BV Dyn. Stall Model - Positive Stall AOA (deg): $aoamaxcl
BV Dyn. Stall Model - Negative Stall AOA (deg): -$aoamaxcl
LB Dyn. Stall Model - Lift Coeff. Slope at Zero Lift AOA (per radian): 0.0
LB Dyn. Stall Model - Positive Critical Lift Coeff.: 0.0
LB Dyn. Stall Model - Negative Critical Lift Coeff.: 0.0
AOA (deg) CL CD Cm25\n")
        for (ialpha,alpha) in enumerate(alpha_ext)
            write(file, "$(alpha_ext[ialpha]) $(cl_ext[ialpha]) $(cd_ext[ialpha]) 0.0\n")
        end
    end
end



# # Generate the new file for AeroDyn
# new_filename = "$(filename[1:end-3])dat"
# open(new_filename, "w") do file
#     # Write the AirfoilInfo header
#     write(file, """
# ! ------------ AirfoilInfo v1.01.x Input File ----------------------------------
# ! Your airfoil description here (replace as needed)
# ! line
# ! line
# ! ------------------------------------------------------------------------------
# DEFAULT       InterpOrd     - ! Interpolation order to use for quasi-steady table lookup {1=linear; 3=cubic spline; "default"} [default=3]
# $(round(maximum(yr)-minimum(yr), digits=4))             RelThickness  - ! The non-dimensional thickness of the airfoil (thickness/chord)
# 1             NonDimArea    - ! The non-dimensional area of the airfoil (area/chord^2)
# @"NACA_0018.txt"         NumCoords   - ! Number of coordinates in the airfoil shape file.
# "unused"      BL_file       - ! The file name including the boundary layer characteristics of the profile.
# $(length(Rerange))            NumTabs       - ! Number of airfoil tables in this file.
# """)

#     for Re in Rerange
#         alpha = -30:1:30
#         c_l_raw, c_d_raw, c_dp_raw, c_m_raw, converged = Xfoil.alpha_sweep(xr, yr, alpha, Re, iter=100, zeroinit=true, printdata=true, reinit=false)

#         c_l = FLOWMath.linear(alpha[converged.==1],c_l_raw[converged.==1],alpha)
#         c_d = FLOWMath.linear(alpha[converged.==1],c_d_raw[converged.==1],alpha)
#         c_dp = FLOWMath.linear(alpha[converged.==1],c_dp_raw[converged.==1],alpha)
#         c_m = FLOWMath.linear(alpha[converged.==1],c_m_raw[converged.==1],alpha)

#         # Viterna extrapolation
#         cr75 = 0.1
#         alpha_ext, cl_ext, cd_ext = CCBlade.viterna(alpha ./ 360.0 .* 2π, c_l, c_d, cr75)
#         alpha_ext = alpha_ext ./ (2π) .* 360

#         # Sort by alpha
#         sort_inds = sortperm(alpha_ext)
#         alpha_ext = alpha_ext[sort_inds]
#         cl_ext = cl_ext[sort_inds]
#         cd_ext = cd_ext[sort_inds]

#         PyPlot.figure("Re: $Re")
#         PyPlot.plot(alpha,c_l)
#         PyPlot.plot(alpha,c_d)
#         PyPlot.plot(alpha,c_l_raw,label="raw")
#         PyPlot.plot(alpha,c_d_raw,label="raw")
#         PyPlot.plot(alpha_ext,cl_ext,".",label="cl_vit")
#         PyPlot.plot(alpha_ext,cd_ext,".",label="cd_vit")
#         PyPlot.legend()


#         # Write header for this Reynolds number table
#         write(file, """
# ! ------------------------------------------------------------------------------
# $(round(Re / 1e6, digits=3))         Re            - ! Reynolds number in millions
# 0             Ctrl          - ! Control setting
# False         InclUAdata    - ! Is unsteady aerodynamics data included?
# !........................................
# ! Table of aerodynamics coefficients
# $(length(alpha_ext))           NumAlf        - ! Number of data lines in the following table
# ! Alpha       Cl            Cd              Cm
# ! (deg)       (-)           (-)             (-)
# """)

#         for i in 1:length(alpha_ext)
#             write(file, @sprintf("%-8.1f\t%.4f\t%.4f\t%.4f\n", alpha_ext[i], cl_ext[i], cd_ext[i], 0.0))
#         end
#     end
# end
