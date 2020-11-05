function initialElementCalculations(model,el,mesh)
    #initialElementCalculations  performs intitial element calculations
    #   [elStorage] = initialElementCalculations(model,el,mesh)
    #
    #   This function performs initial element calculation for use later in
    #   analysis for efficiency gains.
    #
    #      input:
    #      model               = object containing model information
    #      el                  = object containing element information
    #      mesh                = object containing mesh information
    #
    #      output:
    #      elStorage           = object containing stored element data

    #initial element calculation
    numNodesPerEl = 2

    # elStorage = repmat(single_elStorage,mesh.numEl,1)
    elStorage = Array{ElStorage, 1}(undef, mesh.numEl)
    # elStorage = {}

    for i=1:mesh.numEl
        #Calculate Ke and Fe for element i
        elementOrder = model.elementOrder #assign for element i
        modalFlag = true
        xloc = [0.0 el.elLen[i]]
        sectionProps = el.props[i]
        sweepAngle = el.psi[i]
        coneAngle = el.theta[i]
        rollAngle = el.roll[i]
        aeroSweepAngle = 0.0

        elx = zeros(numNodesPerEl)
        ely = zeros(numNodesPerEl)
        elz = zeros(numNodesPerEl)
        for j=1:numNodesPerEl
            #get element cooridnates
            elx[j] = mesh.x[Int(mesh.conn[i,j])]
            ely[j] = mesh.y[Int(mesh.conn[i,j])]
            elz[j] = mesh.z[Int(mesh.conn[i,j])]
        end

        #get concentrated terms associated with elemetn
        massConc,_,_,_,_,_ = ConcMassAssociatedWithElement(mesh.conn[i,:],model.joint,model.nodalTerms.concMass,model.nodalTerms.concStiff,model.nodalTerms.concLoad)

        concMassFlag = !isempty(findall(x->x!=0,massConc))

        Omega = 0.0

        mat"$elStoragein = calculateTimoshenkoElementInitialRun($elementOrder,$modalFlag,$xloc,$sectionProps,$sweepAngle,$coneAngle,$rollAngle,$aeroSweepAngle,$elx,$ely,$elz,$concMassFlag,$massConc,$Omega)" #initial element calculations for storage
        elStorage[i] = ElStorage(elStoragein["K11"],
        elStoragein["K12"],
        elStoragein["K13"],
        elStoragein["K14"],
        elStoragein["K15"],
        elStoragein["K16"],
        elStoragein["K22"],
        elStoragein["K23"],
        elStoragein["K24"],
        elStoragein["K25"],
        elStoragein["K26"],
        elStoragein["K33"],
        elStoragein["K34"],
        elStoragein["K35"],
        elStoragein["K36"],
        elStoragein["K44"],
        elStoragein["K45"],
        elStoragein["K46"],
        elStoragein["K55"],
        elStoragein["K56"],
        elStoragein["K66"],
        elStoragein["M11"],
        elStoragein["M15"],
        elStoragein["M16"],
        elStoragein["M22"],
        elStoragein["M24"],
        elStoragein["M33"],
        elStoragein["M34"],
        elStoragein["M44"],
        elStoragein["M55"],
        elStoragein["M56"],
        elStoragein["M66"],
        elStoragein["S11"],
        elStoragein["S12"],
        elStoragein["S13"],
        elStoragein["S15"],
        elStoragein["S16"],
        elStoragein["S22"],
        elStoragein["S23"],
        elStoragein["S25"],
        elStoragein["S26"],
        elStoragein["S33"],
        elStoragein["S35"],
        elStoragein["S36"],
        elStoragein["S55"],
        elStoragein["S56"],
        elStoragein["S66"],
        elStoragein["S14_1"],
        elStoragein["S14_2"],
        elStoragein["S24_1"],
        elStoragein["S24_2"],
        elStoragein["S34_1"],
        elStoragein["S34_2"],
        elStoragein["S45_1"],
        elStoragein["S45_2"],
        elStoragein["S46_1"],
        elStoragein["S46_2"],
        elStoragein["S44_1"],
        elStoragein["S44_2"],
        elStoragein["S44_3"],
        elStoragein["C12"],
        elStoragein["C13"],
        elStoragein["C23"],
        elStoragein["C24"],
        elStoragein["C25"],
        elStoragein["C26"],
        elStoragein["C34"],
        elStoragein["C35"],
        elStoragein["C36"],
        elStoragein["C14_1"],
        elStoragein["C14_2"],
        elStoragein["C45_1"],
        elStoragein["C45_2"],
        elStoragein["C46_1"],
        elStoragein["C46_2"],
        elStoragein["mel"],
        elStoragein["moiel"],
        elStoragein["xmel"])
    end
    return elStorage
end

function ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,nodalStiffnessTerms,nodalLoads)
    #ConcMassAssociatedWithElement gets concentrated terms associated w/ el
    #   [mass,stiff,load,modJoint,modNodalMassTerms,...
    #    modNodalStiffnessTerms,modNodalLoads] = ...
    #     ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,...
    #     nodalStiffnessTerms,nodalLoads)
    #
    #   This function compiles concentrated mass, stiffness, and load
    #   associated with a node from both ndl and joint files. The mod*
    #   variables are passed back with these terms removed to prevent
    #   duplicate application of shared nodal terms between elements
    #
    #      input:
    #      conn                = connectivity list for element
    #      joint               = joint array for nodal terms
    #      nodalMassTerms      = listing of concentrated nodal mass terms
    #      nodalStiffnessTerms = listing of concentrated nodal stiffness terms
    #      nodalLoads          = listing of concentrated nodal loads terms
    #
    #
    #      output:
    #      mass                = array of concentrated mass associated with element
    #      stiff               = array of concentrated stiffness associated with
    #                            element
    #      load                = array of concentrated loads associated with element
    #      modJoint            = modified joint object removing nodal terms that
    #                            have/will be applied to the element calculations
    #      modNodalMassTerms   = modified nodal mass object removing nodal terms that
    #                            have/will be applied to the element calculations
    #      modalStiffnessTerms = modified nodal stiffness object removing nodal terms that
    #                            have/will be applied to the element calculations
    #      modNodalLoads       = modified nodal loads object removing nodal terms that
    #                            have/will be applied to the element calculations

    node1 = conn[1] #define node #1 and node #2
    node2 = conn[2]

    mass1=0  #initialize concentrated mass amd moi for nodes
    mass2=0
    moix1=0
    moiy1=0
    moiz1=0
    moix2=0
    moiy2=0
    moiz2=0

    stiff1x=0 #initialize concentrated stifness for nodes
    stiff2x=0
    stiff1y=0
    stiff2y=0
    stiff1z=0
    stiff2z=0
    stiff1mx=0
    stiff2mx=0
    stiff1my=0
    stiff2my=0
    stiff1mz=0
    stiff2mz=0

    f1x = 0   #initialize concentrated loads/moments
    f2x = 0
    f1y = 0
    f2y = 0
    f1z = 0
    f2z = 0
    m1x =0
    m2x =0
    m1y =0
    m2y =0
    m1z =0
    m2z =0

    modJoint = joint                         #create copies of joint, and nodal mass, stiffness, loads arrays
    modNodalMassTerms = nodalMassTerms
    modNodalStiffnessTerms = nodalStiffnessTerms
    modNodalLoads = nodalLoads

    numJoints,_=size(joint)    #get number of joints in model

    if numJoints > 0
        node1flag=joint[:,2].==node1  #see if nodes are associated with a joint constraint as a master node
        node2flag=joint[:,2].==node2
    else
        node1flag = false
        node2flag = false
        mass1 = 0.0
        mass2 = 0.0
    end

    for i=1:numJoints           #if nodes are associated with joint constraint, use (if any) mass and stiffness specification from the joint file
        if node1flag[i]==1
            mass1 = mass1+joint[i,5]
            #             stiff1x = stiff1x + joint[i,6]
            #             stiff1y = stiff1y + joint[i,6]
            #             stiff1z = stiff1z + joint[i,6]
            #             stiff1mx = stiff1mx + joint[i,6]
            #             stiff1my = stiff1my + joint[i,6]
            #             stiff1mz = stiff1mz + joint[i,6]
            #             modJoint(i,5)=0.0
            #             modJoint(i,6)=0.0
        end
        if node2flag[i]==1
            mass2 = mass2+joint[i,5]
            #             stiff2x = stiff2x + joint[i,6]
            #             stiff2y = stiff2y + joint[i,6]
            #             stiff2z = stiff2z + joint[i,6]
            #             stiff2mx = stiff2mx + joint[i,6]
            #             stiff2my = stiff2my + joint[i,6]
            #             stiff2mz = stiff2mz + joint[i,6]
            #             modJoint(i,5)=0.0
            #             modJoint(i,6)=0.0
        end
    end

    #apply concentrated mass/stiffness from NDL file

    for i=1:length(nodalMassTerms)   #if node is specified in nodal mass terms file add to mass properties for this node
        node1flagM=nodalMassTerms[i].nodeNum.==node1
        node2flagM=nodalMassTerms[i].nodeNum.==node2
        if node1flagM==1
            if nodalMassTerms[i].dof == 1
                mass1 = mass1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 4
                moix1 = moix1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 5
                moiy1 = moiy1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 6
                moiz1 = moiz1+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
        end
        if node2flagM==1
            mass2 = mass2+nodalMassTerms[i].val
            modNodalMassTerms[i].val = 0.0

            if nodalMassTerms[i].dof == 4
                moix2 = moix2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 5
                moiy2 = moiy2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
            if nodalMassTerms[i].dof == 6
                moiz2 = moiz2+nodalMassTerms[i].val
                modNodalMassTerms[i].val = 0.0
            end
        end
    end



    for i=1:length(nodalStiffnessTerms)     #if node is specified in nodal stiffness terms file add to stiffness properties for this node
        node1flagK=nodalStiffnessTerms[i].nodeNum.==node1
        node2flagK=nodalStiffnessTerms[i].nodeNum.==node2
        if node1flagK==1
            if nodalStiffnessTerms[i].dof==1
                stiff1x = stiff1x+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==2
                stiff1y = stiff1y+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==3
                stiff1z = stiff1z+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==4
                stiff1mx = stiff1mx+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==5
                stiff1my = stiff1my+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==6
                stiff1mz = stiff1mz+nodalStiffnessTerms[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalStiffnessTerms[i].val = 0.0
        end
        if node2flagK==1
            if nodalStiffnessTerms[i].dof==1
                stiff2x = stiff2x+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==2
                stiff2y = stiff2y+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==3
                stiff2z = stiff2z+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==4
                stiff2mx = stiff2mx+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==5
                stiff2my = stiff2my+nodalStiffnessTerms[i].val
            elseif nodalStiffnessTerms[i].dof==6
                stiff2mz = stiff2mz+nodalStiffnessTerms[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalStiffnessTerms[i].val = 0.0
        end
    end

    for i=1:length(nodalLoads)  #if node is specified in nodal forces terms file add to concentrated force for this node
        node1flagF=nodalLoads[i].nodeNum.==node1
        node2flagF=nodalLoads[i].nodeNum.==node2
        if node1flagF==1
            if nodalLoads[i].dof==1
                f1x = f1x+nodalLoads[i].val
            elseif nodalLoads[i].dof==2
                f1y = f1y+nodalLoads[i].val
            elseif nodalLoads[i].dof==3
                f1z = f1z+nodalLoads[i].val
            elseif nodalLoads[i].dof==4
                m1x = m1x+nodalLoads[i].val
            elseif nodalLoads[i].dof==5
                m1y = m1y+nodalLoads[i].val
            elseif nodalLoads[i].dof==6
                m1z = m1z+nodalLoads[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalLoads[i].val = 0.0
        end
        if node2flagF==1
            if nodalLoads[i].dof==1
                f2x = f2x+nodalLoads[i].val
            elseif nodalLoads[i].dof==2
                f2y = f2y+nodalLoads[i].val
            elseif nodalLoads[i].dof==3
                f2z = f2z+nodalLoads[i].val
            elseif nodalLoads[i].dof==4
                m2x = m2x+nodalLoads[i].val
            elseif nodalLoads[i].dof==5
                m2y = m2y+nodalLoads[i].val
            elseif nodalLoads[i].dof==6
                m2z = m2z+nodalLoads[i].val
            else
                error("DOF not valid for  concentrated stiffness term.")
            end
            modNodalLoads[i].val = 0.0
        end
    end


    #compile nodal concentrated terms into mass, stiffness, and load arrays
    mass = [mass1 mass2
    moix1 moix2
    moiy1 moiy2
    moiz1 moiz2]

    stiff = [stiff1x stiff2x
    stiff1y stiff2y
    stiff1z stiff2z
    stiff1mx stiff2mx
    stiff1my stiff2my
    stiff1mz stiff2mz]

    load = [f1x f2x
    f1y f2y
    f1z f2z
    m1x m2x
    m1y m2y
    m1z m2z]

    return mass,stiff,load,modJoint,modNodalMassTerms,modNodalStiffnessTerms,modNodalLoads

end
