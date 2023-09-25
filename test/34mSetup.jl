import QuadGK
##############################################
# Setup Structures
#############################################
starttime = time()
# SNL34_unit_xz = DelimitedFiles.readdlm("$(path)/data/SNL34m_unit_blade_shape.txt",'\t',skipstart = 0)
# include("$(path)/34mBladeshapeAnalytical2.jl")
# Optimized
# z_shape = collect(LinRange(0,41.9,12))
# x_shape =[0.0, 6.366614283571894, 10.845274468506549, 14.262299623370412, 16.295669420203165, 17.329613328115716, 17.158774783226303, 15.813537149178769, 13.479124754351849, 10.04333990055769, 5.606817958279066,0.0]
# Old 3.5%
# controlpts =[3.256322421104705, 5.774184885976885, 8.647885597459139, 11.186664047211988, 13.145770457622106, 14.674641597073201, 15.804683728982331, 16.61865064561603, 17.108011791043822, 17.296713936456687, 17.19399803467357, 16.803560229287708, 16.108240799302397, 15.118450192780204, 13.803814107334938, 12.21986771504199, 10.359459126160356, 8.205302489986666, 5.765600261682509, 2.975874178673999]
# New 2.15#
controlpts = [3.6479257474344826, 6.226656883619295, 9.082267631309085, 11.449336766507562, 13.310226748873827, 14.781369210504563, 15.8101544043681, 16.566733104331984, 17.011239869982738, 17.167841319391137, 17.04306679619916, 16.631562597633675, 15.923729603782338, 14.932185789551408, 13.62712239754136, 12.075292152969496, 10.252043906945818, 8.124505683235517, 5.678738418596312, 2.8959968657512207]

# z_shape = collect(LinRange(0,41.9,length(x_shape)))
z_shape1 = collect(LinRange(0,41.9,length(controlpts)+2))
x_shape1 = [0.0;controlpts;0.0]
z_shape = collect(LinRange(0,41.9,60))
x_shape = FLOWMath.akima(z_shape1,x_shape1,z_shape)#[0.0,1.7760245854312287, 5.597183088188207, 8.807794161662574, 11.329376903432605, 13.359580331518579, 14.833606099357858, 15.945156349709, 16.679839160110422, 17.06449826588358, 17.10416552269884, 16.760632435904647, 16.05982913536134, 15.02659565585254, 13.660910465851046, 11.913532434360155, 9.832615229216344, 7.421713825584581, 4.447602800040282, 0.0]
toweroffset = 4.3953443986241725
# Analytical
# z_shape = [0.0, 0.4027099927689326, 0.8054199855378652, 1.2081299783067978, 1.6108399710757304, 2.013549963844663, 2.4162599566135956, 2.818969949382528, 3.221679942151461, 3.624389934920394, 4.027099927689326, 4.429809920458259, 4.832519913227191, 5.235229905996124, 5.637939898765056, 6.221683770581164, 6.821719178571302, 7.437583162221221, 8.068800548394838, 8.714884317956601, 9.375335981533686, 10.049645964128134, 10.737293998282153, 11.437749525493246, 12.305529745531164, 13.201631116890074, 14.122956691386433, 15.066322345375763, 16.028467784161606, 17.00606780965343, 17.995743812332336, 18.994075447807884, 19.997612457611687, 21.002886593374797, 22.006423603178604, 23.004755238654155, 23.99443124133306, 24.97203126682489, 25.934176705610724, 26.877542359600053, 27.79886793409642, 28.694969305455324, 29.562749525493246, 30.263205052704336, 30.950853086858352, 31.62516306945281, 32.285614733029895, 32.93169850259165, 33.56291588876527, 34.1787798724152, 34.77881528040533, 35.362559152221436, 35.821422444858186, 36.280285737494935, 36.739149030131685, 37.198012322768435, 37.656875615405184, 38.11573890804193, 38.57460220067868, 39.03346549331543, 39.49232878595218, 39.951192078588925, 40.410055371225674, 40.868918663862424, 41.32778195649917, 41.78664524913592]
# x_shape = [0.0, 0.6201190084429031, 1.2402380168858063, 1.8603570253287094, 2.4804760337716125, 3.1005950422145157, 3.720714050657419, 4.340833059100322, 4.960952067543225, 5.581071075986128, 6.201190084429031, 6.821309092871934, 7.441428101314838, 8.061547109757742, 8.681666118200644, 9.276344926168852, 9.854581297981857, 10.41592909228786, 10.959955198206961, 11.48623986950018, 11.994377048426914, 12.4839746790409, 12.954655009683117, 13.406054884438076, 13.913532546749641, 14.369140295018864, 14.771303537761444, 15.118632389006988, 15.409926471778816, 15.644179066626025, 15.820580590870494, 15.938521396544312, 15.997593877348054, 15.997593877348054, 15.938521396544312, 15.820580590870494, 15.644179066626029, 15.409926471778816, 15.118632389006992, 14.771303537761447, 14.369140295018864, 13.913532546749645, 13.406054884438078, 12.954655009683123, 12.483974679040909, 11.994377048426916, 11.486239869500185, 10.959955198206966, 10.415929092287865, 9.85458129798186, 9.276344926168857, 8.681666118200653, 8.061547968581657, 7.441429818962661, 6.821311669343665, 6.201193519724669, 5.581075370105673, 4.960957220486675, 4.340839070867679, 3.7207209212486827, 3.1006027716296867, 2.480484622010689, 1.8603664723916928, 1.2402483227726968, 0.6201301731537008, 1.2023534704752592e-5]
# y_shape = zero(x_shape)

SNL34_unit_xz = [x_shape;;z_shape]
#Ensure the data is fully unitized since the hand picking process is only good to one or two significant digits.
SNL34x = SNL34_unit_xz[:,1]./maximum(SNL34_unit_xz[:,1])
SNL34z = SNL34_unit_xz[:,2]./maximum(SNL34_unit_xz[:,2])

#Scale the turbine to the full dimensions
height = 41.9 #m
radius = 17.1 #m
SNL34Z = SNL34z.*height
SNL34X = SNL34x.*radius

Nbld = 2

mymesh,myort,myjoint = OWENS.create_mesh_struts(;Ht=0.5,
    Hb = height, #blade height
    R = radius, # m bade radius
    nblade = Nbld,
    ntelem = 20, #tower elements
    nbelem = 60, #blade elements
    nselem = 3,
    strut_mountpointbot = 0.03,
    strut_mountpointtop = 0.03,
    bshapex = SNL34X,#cos.(LinRange(0,0,12)).*SNL34X, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = SNL34Z, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    # bshapey = sin.(LinRange(0,0,12)).*SNL34X,
    angularOffset = pi/2)

# PyPlot.figure()
# PyPlot.plot(mymesh.x,mymesh.z,"b-")
#  for myi = 1:length(mymesh.y)
#      PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.z[myi].+rand()/30,"$myi",ha="center",va="center")
#      PyPlot.draw()
#      sleep(0.1)
#  end
# PyPlot.xlabel("x")
# PyPlot.ylabel("y")
# # PyPlot.axis("equal")
# visualize_orts = true
# if visualize_orts
#     import PyPlot
#     PyPlot.pygui(true)
#     PyPlot.rc("figure", figsize=(15, 15))
#     PyPlot.rc("font", size=10.0)
#     PyPlot.rc("lines", linewidth=1.5)
#     PyPlot.rc("lines", markersize=4.0)
#     PyPlot.rc("legend", frameon=true)
#     PyPlot.rc("axes.spines", right=false, top=false)
#     PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
#     PyPlot.rc("figure",max_open_warning=500)
#     # PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
#     plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

#     ####################################################################################
#     ################## Plot for the elements #########################
#     ####################################################################################


#     PyPlot.figure()
#     PyPlot.scatter3D(mymesh.x,mymesh.y,mymesh.z,color=plot_cycle[1])
#     # PyPlot.zlim(([90,120]))

#     # Allocate the node angle arrays
#     Psi_d_Nodes = zeros(mymesh.numNodes)
#     Theta_d_Nodes = zeros(mymesh.numNodes)
#     Twist_d_Nodes = zeros(mymesh.numNodes)

#     function rotate_normal(i_el, mymesh, myort;vec=[0,1,0.0],normal_len=1)
#         # * `Psi_d::Array{<:float}`: length NumEl, element rotation about 3 in global FOR (deg) These angles are used to transform from the global coordinate frame to the local element/joint frame via a 3-2 Euler rotation sequence.
#         # * `Theta_d::Array{<:float}`: length NumEl, element rotation about 2 (deg)
#         # * `Twist_d::Array{<:float}`: length NumEl, element twist (deg)
        
#         # Map from element to node
#         nodenum1 = Int(mymesh.conn[i_el,1])
#         nodenum2 = Int(mymesh.conn[i_el,2])

#         # Extract the element angles
#         Psi_d_el = myort.Psi_d[i_el]
#         Theta_d_el = myort.Theta_d[i_el]
#         Twist_d_el = myort.Twist_d[i_el]

#         # Map the element angles to the node angles
#         Psi_d_Nodes[[nodenum1,nodenum2]] .= Psi_d_el
#         Theta_d_Nodes[[nodenum1,nodenum2]] .= Theta_d_el
#         Twist_d_Nodes[[nodenum1,nodenum2]] .= Twist_d_el

#         # Use a line and rotate it about the angles, different starting vectors show different angles.
#         myvec = vec.*normal_len

#         # apply the twist rotation, which is about the x (1) axis
#         DCM_roll = [1.0 0.0 0.0
#             0.0 cosd(Twist_d_el) -sind(Twist_d_el)
#             0.0 sind(Twist_d_el) cosd(Twist_d_el)]

#         # apply theta rotation, which is the tilt angle, or about the y (2) axis in global
#         DCM_pitch = [cosd(Theta_d_el) 0.0 sind(Theta_d_el)
#             0.0 1.0 0.0
#             -sind(Theta_d_el) 0.0 cosd(Theta_d_el)]

#         # apply Psi rotation, which is about Z (3) axis in global
#         DCM_yaw = [cosd(Psi_d_el) -sind(Psi_d_el) 0.0
#             sind(Psi_d_el) cosd(Psi_d_el) 0.0
#             0.0 0.0 1.0]

#         # Get the location of the element
#         x_el = (mymesh.x[nodenum1]+mymesh.x[nodenum2])/2
#         y_el = (mymesh.y[nodenum1]+mymesh.y[nodenum2])/2
#         z_el = (mymesh.z[nodenum1]+mymesh.z[nodenum2])/2

#         # Offset the myvector by the location of the element
#         myvec = DCM_yaw*DCM_pitch*DCM_roll*myvec + [x_el,y_el,z_el]
#         x_el_plot = [x_el,myvec[1]]
#         y_el_plot = [y_el,myvec[2]]
#         z_el_plot = [z_el,myvec[3]]
#         return x_el_plot, y_el_plot, z_el_plot
#     end

#     # Add the orientation vectors, ort is on elements
#     for i_el = 1:mymesh.numEl
#         x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, mymesh, myort;vec=[10,0,0.0])
#         PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[2])
#     end

#     for i_el = 1:mymesh.numEl
#         x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, mymesh, myort;vec=[0,5,0.0])
#         PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[3])
#     end

#     for i_el = 1:mymesh.numEl
#         x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, mymesh, myort;vec=[0,0,5.0])
#         PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[4])
#     end

#     PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[2],label="X-norm")
#     PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[3],label="Y-norm")
#     PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[4],label="Z-norm")
#     PyPlot.legend()
#     # for i_joint = 1:length(myjoint[:,1])
#     #     i_el = findall(x->x==myjoint[i_joint,2],mymesh.conn[:,1])
#     #     if length(i_el)==0 #Use the other element associated with the joint
#     #         i_el = findall(x->x==myjoint[i_joint,2],mymesh.conn[:,2])
#     #     end
#     #     if length(i_el)==0 #Use the other element associated with the joint
#     #         i_el = findall(x->x==myjoint[i_joint,3],mymesh.conn[:,2])
#     #     end
#     #     x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el[1], mymesh, myort;normal_len=3)
#     #     PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[5])
#     # end

#     PyPlot.xlabel("x")
#     PyPlot.ylabel("y")
#     PyPlot.zlabel("z")
#     PyPlot.axis("equal")
# end

nTwrElem = Int(mymesh.meshSeg[1])+2
nBldElem = Int(mymesh.meshSeg[2])+1
#Blades
NuMad_geom_xlscsv_file = "$path/data/SNL34mGeom.csv"
numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

for (i,airfoil) in enumerate(numadIn_bld.airfoil)
    numadIn_bld.airfoil[i] = "$path/airfoils/$airfoil"
end

NuMad_mat_xlscsv_file = "$path/data/SNL34mMaterials.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

bld1start = Int(mymesh.structuralNodeNumbers[1,1])
bld1end = Int(mymesh.structuralNodeNumbers[1,end])
spanpos = [0.0;cumsum(sqrt.(diff(mymesh.x[bld1start:bld1end]).^2 .+ diff(mymesh.z[bld1start:bld1end]).^2))]

bld_precompoutput,bld_precompinput = OWENS.getPreCompOutput(numadIn_bld;plyprops)
sectionPropsArray_bld = OWENS.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)

stiff_bld, mass_bld = OWENS.getSectPropsFromPreComp(spanpos,numadIn_bld,bld_precompoutput;GX=true)

# thickness_flap is distance from shear center x to top
# thickness_lag is distance from shear center y to trailing edge
# shear center is relative to the blade reference axis
# blade reference axis is from leading edge to le_loc*chord and along the chord line
# reference axes Y is along the chord, and X is perpendicular to the chord

thickness_precomp_lag = zeros(length(bld_precompinput))
thickness_precomp_flap = zeros(length(bld_precompinput))
for ipc = 1:length(bld_precompinput)
    refY = bld_precompinput[ipc].le_loc*bld_precompinput[ipc].chord
                                # Negative distance for lag, to align with SAND-88-1144
    thickness_precomp_lag[ipc] = -(bld_precompinput[ipc].chord-(refY+bld_precompoutput[ipc].y_sc))
    thickness_precomp_flap[ipc] = maximum(bld_precompinput[ipc].ynode)*bld_precompinput[ipc].chord - bld_precompoutput[ipc].x_sc
end
# PyPlot.figure()
# PyPlot.plot(bld_precompinput[1].ynode,bld_precompinput[1].ynode)
spanposmid = cumsum(diff(spanpos))
thickness = FLOWMath.akima(numadIn_bld.span,thickness_precomp_flap,spanposmid)
thickness_lag = FLOWMath.akima(numadIn_bld.span,thickness_precomp_lag,spanposmid)
# thickness = thicknessGX[1:end-1]


NuMad_geom_xlscsv_file = "$path/data/NuMAD_34m_TowerGeom.csv"
numadIn = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

for (i,airfoil) in enumerate(numadIn.airfoil)
    numadIn.airfoil[i] = "$path/airfoils/$airfoil"
end

NuMad_mat_xlscsv_file = "$path/data/NuMAD_34m_TowerMaterials.csv"
plyprops = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

precompoutput,precompinput = OWENS.getPreCompOutput(numadIn;plyprops)
sectionPropsArray_twr = OWENS.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn,precompoutput;precompinputs=precompinput)

stiff_twr, mass_twr = OWENS.getSectPropsFromPreComp(LinRange(0,1,nTwrElem),numadIn,precompoutput;GX=true)
#Struts
# They are the same as the end properties of the blades

# Combined Section Props
bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
Nremain = mymesh.numEl-length(sectionPropsArray_twr)-length(bldssecprops) #strut elements remain
sectionPropsArray = [sectionPropsArray_twr;bldssecprops;fill(sectionPropsArray_bld[end],Nremain)]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]

rotationalEffects = ones(mymesh.numEl)

for i = 1:length(sectionPropsArray)
    # sectionPropsArray[i].rhoA .*= 0.25
    # sectionPropsArray[i].EIyy .*= 5.0
    # sectionPropsArray[i].EIzz .*= 5.0
    # sectionPropsArray[i].EIyz .*= 5.0
    # sectionPropsArray[i].GJ .*= 5.0
    # sectionPropsArray[i].EA .*= 5.0
    # sectionPropsArray[i].rhoIyy .*= 0.1
    # sectionPropsArray[i].rhoIzz .*= 0.1
    # sectionPropsArray[i].rhoIyz .*= 0.1
    # sectionPropsArray[i].rhoJ .*= 0.1
end


#store data in element object
myel = GyricFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

top_idx = 23#Int(myjoint[7,2])
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0
top_idx 1 0
top_idx 2 0
top_idx 3 0
top_idx 4 0
top_idx 5 0]
# top_idx 6 0]

##############################################
# Setup Aero
#############################################


shapeX_spline = FLOWMath.Akima(SNL34Z, SNL34X)
RefArea_half, error = QuadGK.quadgk(shapeX_spline, 0, height, atol=1e-10)
RefArea = RefArea_half*2*0.95

B = 2
Nslices = 35

T1 = round(Int,(5.8/height)*Nslices)
T2 = round(Int,(11.1/height)*Nslices)
T3 = round(Int,(29.0/height)*Nslices)
T4 = round(Int,(34.7/height)*Nslices)

airfoils = fill("$(path)/airfoils/NACA_0021.dat",Nslices)
airfoils[T1:T4] .= "$(path)/airfoils/Sandia_001850.dat"

chord = fill(1.22,Nslices)
chord[T1:T4] .= 1.07
chord[T2:T3] .= 0.9191
#TODO: when twist introduced, aero should pull it from the precomp input data
rho = 0.94 # for texas site (3880 ft) at 80F

RPM = 34.0
omega = RPM*2*pi/60

# filename = "$(path)/data/legacyfiles/SNL34m"
# OWENS.saveOWENSfiles(filename,mymesh,myort,myjoint,myel,pBC,numadIn_bld)

# system, assembly, sections, frames_ow = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,mass_twr, mass_bld, stiff_twr, stiff_bld;VTKmeshfilename="$path/vtk/SNL34m")

function runmeOWENS()
    mass = 0.0
    for (i,sectionProp) in enumerate(sectionPropsArray)
        lenEl = 0
        try
            lenEl = myort.Length[i]
        catch
            lenEl = myort.Length[i-1]
        end
        rhoA = sectionProp.rhoA[1]
        mass += lenEl*rhoA
    end
return mass
end

massOwens = runmeOWENS()
#
# function runmeGX()
#     mass = 0.0
#     for element in assembly.elements
#         mass += element.L*element.mass[1,1]
#     end
# return mass
# end
#
# massGX = runmeGX()
#
println("Mass")
println("massOwens $massOwens")
# println("massGX $massGX")
#
# ei_flap_exp_base = 1.649e7
# ei_flap_exp_mid = 5.197e6
# ei_flap_exp_center = 3.153e6
#
# ei_lag_exp_base = 2.744e8
# ei_lag_exp_mid = 1.125e8
# ei_lag_exp_center = 6.674e7
#
# ea_exp_base = 2.632e9
# ea_exp_mid = 1.478e9
# ea_exp_center = 1.185e9
#
#
# println("1 ei_flap: $(bld_precompoutput[1].ei_flap) exp: $ei_flap_exp_base")
# println("1 ei_lag: $(bld_precompoutput[1].ei_lag) exp: $ei_lag_exp_base")
# println("1 ea: $(bld_precompoutput[1].ea) exp: $ea_exp_base")
# println()
# println("5 ei_flap: $(bld_precompoutput[5].ei_flap) exp: $ei_flap_exp_mid")
# println("5 ei_lag: $(bld_precompoutput[5].ei_lag) exp: $ei_lag_exp_mid")
# println("5 ea: $(bld_precompoutput[5].ea) exp: $ea_exp_mid")
# println()
# println("12 ei_flap: $(bld_precompoutput[12].ei_flap) exp: $ei_flap_exp_center")
# println("12 ei_lag: $(bld_precompoutput[12].ei_lag) exp: $ei_lag_exp_center")
# println("12 ea: $(bld_precompoutput[12].ea) exp: $ea_exp_center")
# println()

function runowens(model,feamodel,mymesh,myel,aeroForcesDMS,deformTurb;steady=true,system=nothing,assembly=nothing,VTKFilename="./outvtk")

    if !steady
        println("running unsteady")

        t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
        FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
        epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist = OWENS.Unsteady(model;
        topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForcesDMS,deformAero=deformTurb,system,assembly)

        meanepsilon_z_hist = mean(epsilon_z_hist,dims=1)
        meanepsilon_y_hist = mean(epsilon_y_hist,dims=1)

        # if system != nothing
        #     println("Saving VTK time domain files")
        #     OWENS.gyricFEA_VTK(VTKFilename,t,uHist,system,assembly,sections;scaling=1,azi=aziHist)

        # end

    else
        println("here")
        println("running steady")

        feamodel.analysisType = "S"

        displ=zeros(mymesh.numNodes*6)
        elStorage = GyricFEA.initialElementCalculations(feamodel,myel,mymesh)
        displ,elStrain,staticAnalysisSuccessful,FReaction = GyricFEA.staticAnalysis(feamodel,mymesh,myel,displ,model.OmegaInit,model.OmegaInit,elStorage)

        # format to match the unsteady method
        eps_x = [elStrain[i].epsilon_x[1] for i = 1:length(elStrain)]
        epsilon_x_hist = zeros(1,length(eps_x),2)
        epsilon_x_hist[1,:,1] = eps_x
        epsilon_x_hist[1,:,2] = eps_x

        eps_y1_OW = [elStrain[i].epsilon_y[1] for i = 1:length(elStrain)]
        eps_y2_OW = [elStrain[i].epsilon_y[2] for i = 1:length(elStrain)]
        eps_y3_OW = [elStrain[i].epsilon_y[3] for i = 1:length(elStrain)]
        eps_y4_OW = [elStrain[i].epsilon_y[4] for i = 1:length(elStrain)]
        eps_y = (eps_y1_OW.+eps_y2_OW.+eps_y3_OW.+eps_y4_OW).*0.25#0.34785484513745385
        meanepsilon_y_hist = zeros(1,length(eps_x),2)
        meanepsilon_y_hist[1,:,1] = eps_y
        meanepsilon_y_hist[1,:,2] = eps_y

        eps_z1_OW = [elStrain[i].epsilon_z[1] for i = 1:length(elStrain)]
        eps_z2_OW = [elStrain[i].epsilon_z[2] for i = 1:length(elStrain)]
        eps_z3_OW = [elStrain[i].epsilon_z[3] for i = 1:length(elStrain)]
        eps_z4_OW = [elStrain[i].epsilon_z[4] for i = 1:length(elStrain)]
        eps_z = (eps_z1_OW.+eps_z2_OW.+eps_z3_OW.+eps_z4_OW).*0.25#0.34785484513745385
        meanepsilon_z_hist = zeros(1,length(eps_x),2)
        meanepsilon_z_hist[1,:,1] = eps_z
        meanepsilon_z_hist[1,:,2] = eps_z

        kappa_x = [elStrain[i].kappa_x[1] for i = 1:length(elStrain)]
        kappa_x_hist = zeros(1,length(eps_x),2)
        kappa_x_hist[1,:,1] = kappa_x
        kappa_x_hist[1,:,2] = kappa_x

        kappa_y = [elStrain[i].kappa_y[1] for i = 1:length(elStrain)]
        kappa_y_hist = zeros(1,length(eps_x),2)
        kappa_y_hist[1,:,1] = kappa_y
        kappa_y_hist[1,:,2] = kappa_y

        kappa_z = [elStrain[i].kappa_z[1] for i = 1:length(elStrain)]
        kappa_z_hist = zeros(1,length(eps_x),2)
        kappa_z_hist[1,:,1] = kappa_z
        kappa_z_hist[1,:,2] = kappa_z

        FReactionHist = zeros(2,6)
        FReactionHist[1,:] = FReaction
        FReactionHist[2,:] = FReaction

        OmegaHist = [model.OmegaInit,model.OmegaInit]
        genTorque = FReactionHist[:,6]
        t = [0.0,1.0]
        torqueDriveShaft = [0.0]
        aziHist = [0.0]
        uHist = [0.0]
    end


    # Interpolate the mesh strains onto the composite layup
    # TODO: or should we interpolate the composite stations onto the mesh?  It would be much more challenging

    N_ts = length(epsilon_x_hist[1,1,:])
    eps_x = zeros(Nbld,N_ts,mymesh.meshSeg[2]+1)
    eps_z = zeros(Nbld,N_ts,mymesh.meshSeg[2]+1)
    eps_y = zeros(Nbld,N_ts,mymesh.meshSeg[2]+1)
    kappa_x = zeros(Nbld,N_ts,mymesh.meshSeg[2]+1)
    kappa_y = zeros(Nbld,N_ts,mymesh.meshSeg[2]+1)
    kappa_z = zeros(Nbld,N_ts,mymesh.meshSeg[2]+1)

    for ibld = 1:Nbld
        start = Int(mymesh.structuralNodeNumbers[ibld,1])
        stop = Int(mymesh.structuralNodeNumbers[ibld,end])
        x = mymesh.z[start:stop]
        x = x.-x[1] #zero
        x = x./x[end] #normalize
        # samplepts = numadIn_bld.span./maximum(numadIn_bld.span) #normalize #TODO: this is spanwise, while everything else is vertical-wise
        for its = 1:N_ts
            #TODO: there are strain values at each quad point, should be better than just choosing one
            eps_x[ibld,its,:] = epsilon_x_hist[1,start:stop,its]#FLOWMath.akima(x,epsilon_x_hist[1,start:stop,its],samplepts)
            eps_z[ibld,its,:] = meanepsilon_z_hist[1,start:stop,its]#FLOWMath.akima(x,meanepsilon_z_hist[1,start:stop,its],samplepts)
            eps_y[ibld,its,:] = meanepsilon_y_hist[1,start:stop,its]#FLOWMath.akima(x,meanepsilon_y_hist[1,start:stop,its],samplepts)
            kappa_x[ibld,its,:] = kappa_x_hist[1,start:stop,its]#FLOWMath.akima(x,kappa_x_hist[1,start:stop,its],samplepts)
            kappa_y[ibld,its,:] = kappa_y_hist[1,start:stop,its]#FLOWMath.akima(x,kappa_y_hist[1,start:stop,its],samplepts)
            kappa_z[ibld,its,:] = kappa_z_hist[1,start:stop,its]#FLOWMath.akima(x,kappa_z_hist[1,start:stop,its],samplepts)
        end
    end

    # PyPlot.figure()
    # PyPlot.plot(t[1:end-1],eps_x[1,:,15],label="eps_x")
    # PyPlot.plot(t[1:end-1],eps_z[1,:,15],label="eps_z")
    # PyPlot.plot(t[1:end-1],eps_y[1,:,15],label="eps_y")
    # PyPlot.plot(t[1:end-1],kappa_x[1,:,15],label="kappa_x")
    # PyPlot.plot(t[1:end-1],kappa_y[1,:,15],label="kappa_y")
    # PyPlot.plot(t[1:end-1],kappa_z[1,:,15],label="kappa_z")
    #
    # PyPlot.plot(t[1:end-1],eps_x[2,:,15],":",label="eps_x2")
    # PyPlot.plot(t[1:end-1],eps_z[2,:,15],":",label="eps_z2")
    # PyPlot.plot(t[1:end-1],eps_y[2,:,15],":",label="eps_y2")
    # PyPlot.plot(t[1:end-1],kappa_x[2,:,15],":",label="kappa_x2")
    # PyPlot.plot(t[1:end-1],kappa_y[2,:,15],":",label="kappa_y2")
    # PyPlot.plot(t[1:end-1],kappa_z[2,:,15],":",label="kappa_z2")
    # PyPlot.legend()

    return eps_x,eps_z,eps_y,kappa_x,kappa_y,kappa_z,t,FReactionHist,OmegaHist,genTorque,torqueDriveShaft,aziHist,uHist
end
