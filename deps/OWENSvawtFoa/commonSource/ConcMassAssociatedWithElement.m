function [mass,stiff,load,modJoint,modNodalMassTerms,modNodalStiffnessTerms,modNodalLoads] = ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,nodalStiffnessTerms,nodalLoads)
%ConcMassAssociatedWithElement gets concentrated terms associated w/ el
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [mass,stiff,load,modJoint,modNodalMassTerms,...
%    modNodalStiffnessTerms,modNodalLoads] = ...
%     ConcMassAssociatedWithElement(conn,joint,nodalMassTerms,...
%     nodalStiffnessTerms,nodalLoads)
%                    
%   This function compiles concentrated mass, stiffness, and load 
%   associated with a node from both ndl and joint files. The mod* 
%   variables are passed back with these terms removed to prevent
%   duplicate application of shared nodal terms between elements
%
%      input:
%      conn                = connectivity list for element
%      joint               = joint array for nodal terms
%      nodalMassTerms      = listing of concentrated nodal mass terms
%      nodalStiffnessTerms = listing of concentrated nodal stiffness terms
%      nodalLoads          = listing of concentrated nodal loads terms
%
%
%      output:
%      mass                = array of concentrated mass associated with element
%      stiff               = array of concentrated stiffness associated with 
%                            element
%      load                = array of concentrated loads associated with element
%      modJoint            = modified joint object removing nodal terms that
%                            have/will be applied to the element calculations
%      modNodalMassTerms   = modified nodal mass object removing nodal terms that
%                            have/will be applied to the element calculations
%      modalStiffnessTerms = modified nodal stiffness object removing nodal terms that
%                            have/will be applied to the element calculations
%      modNodalLoads       = modified nodal loads object removing nodal terms that
%                            have/will be applied to the element calculations

    node1 = conn(1); %define node #1 and node #2
    node2 = conn(2);
    
    mass1=0;  %initialize concentrated mass amd moi for nodes
    mass2=0;
    moix1=0;
    moiy1=0;
    moiz1=0;
    moix2=0;
    moiy2=0;
    moiz2=0;
    
    stiff1x=0; %initialize concentrated stifness for nodes
    stiff2x=0;
    stiff1y=0;
    stiff2y=0;
    stiff1z=0;
    stiff2z=0;
    stiff1mx=0;
    stiff2mx=0;
    stiff1my=0;
    stiff2my=0;
    stiff1mz=0;
    stiff2mz=0;
    
    f1x = 0;   %initialize concentrated loads/moments
    f2x = 0;
    f1y = 0;
    f2y = 0;
    f1z = 0;
    f2z = 0;
    m1x =0;
    m2x =0;
    m1y =0;
    m2y =0;
    m1z =0;
    m2z =0;
    
    modJoint = joint;                         %create copies of joint, and nodal mass, stiffness, loads arrays
    modNodalMassTerms = nodalMassTerms;
    modNodalStiffnessTerms = nodalStiffnessTerms;
    modNodalLoads = nodalLoads;
    
    [numJoints,dum]=size(joint);    %get number of joints in model
    
    if(numJoints > 0)
    node1flag=ismember(joint(:,2),node1);  %see if nodes are associated with a joint constraint as a master node
    node2flag=ismember(joint(:,2),node2);
    end
    
    for i=1:numJoints           %if nodes are associated with joint constraint, use (if any) mass and stiffness specification from the joint file
        if(node1flag(i)==1)
            mass1 = mass1+joint(i,5);
%             stiff1x = stiff1x + joint(i,6);
%             stiff1y = stiff1y + joint(i,6);
%             stiff1z = stiff1z + joint(i,6);
%             stiff1mx = stiff1mx + joint(i,6);
%             stiff1my = stiff1my + joint(i,6);
%             stiff1mz = stiff1mz + joint(i,6);
%             modJoint(i,5)=0.0;
%             modJoint(i,6)=0.0;
        end
        if(node2flag(i)==1)
            mass2 = mass2+joint(i,5);
%             stiff2x = stiff2x + joint(i,6);
%             stiff2y = stiff2y + joint(i,6);
%             stiff2z = stiff2z + joint(i,6);
%             stiff2mx = stiff2mx + joint(i,6);
%             stiff2my = stiff2my + joint(i,6);
%             stiff2mz = stiff2mz + joint(i,6);
%             modJoint(i,5)=0.0;
%             modJoint(i,6)=0.0;
        end
    end

    
    
    
    %apply concentrated mass/stiffness from NDL file
    
   
    
    for i=1:length(nodalMassTerms)   %if node is specified in nodal mass terms file add to mass properties for this node
         node1flagM=ismember(nodalMassTerms(i).nodeNum,node1);
         node2flagM=ismember(nodalMassTerms(i).nodeNum,node2);
        if(node1flagM==1)
            if(nodalMassTerms(i).dof == 1)
                mass1 = mass1+nodalMassTerms(i).val;
                modNodalMassTerms(i).val = 0.0;
            end
            if(nodalMassTerms(i).dof == 4)
                moix1 = moix1+nodalMassTerms(i).val;
                modNodalMassTerms(i).val = 0.0;
            end
            if(nodalMassTerms(i).dof == 5)
                moiy1 = moiy1+nodalMassTerms(i).val;
                modNodalMassTerms(i).val = 0.0; 
            end
            if(nodalMassTerms(i).dof == 6)
                moiz1 = moiz1+nodalMassTerms(i).val;
                modNodalMassTerms(i).val = 0.0;  
            end
        end
        if(node2flagM==1)
            mass2 = mass2+nodalMassTerms(i).val;
            modNodalMassTerms(i).val = 0.0;

            if(nodalMassTerms(i).dof == 4)
                moix2 = moix2+nodalMassTerms(i).val;
                modNodalMassTerms(i).val = 0.0;
            end
            if(nodalMassTerms(i).dof == 5)
                moiy2 = moiy2+nodalMassTerms(i).val;
                modNodalMassTerms(i).val = 0.0; 
            end
            if(nodalMassTerms(i).dof == 6)
                moiz2 = moiz2+nodalMassTerms(i).val;
                modNodalMassTerms(i).val = 0.0;  
            end
        end
    end
    
    
    
    for i=1:length(nodalStiffnessTerms)     %if node is specified in nodal stiffness terms file add to stiffness properties for this node
         node1flagK=ismember(nodalStiffnessTerms(i).nodeNum,node1);
         node2flagK=ismember(nodalStiffnessTerms(i).nodeNum,node2);
        if(node1flagK==1)
            if(nodalStiffnessTerms(i).dof==1)
                stiff1x = stiff1x+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==2)
                stiff1y = stiff1y+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==3)
                stiff1z = stiff1z+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==4)
                stiff1mx = stiff1mx+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==5)
                stiff1my = stiff1my+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==6)
                stiff1mz = stiff1mz+nodalStiffnessTerms(i).val;
            else
                error('DOF not valid for  concentrated stiffness term.');
            end
            modNodalStiffnessTerms(i).val = 0.0;
        end
        if(node2flagK==1)
            if(nodalStiffnessTerms(i).dof==1)
                stiff2x = stiff2x+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==2)
                stiff2y = stiff2y+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==3)
                stiff2z = stiff2z+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==4)
                stiff2mx = stiff2mx+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==5)
                stiff2my = stiff2my+nodalStiffnessTerms(i).val;
            elseif(nodalStiffnessTerms(i).dof==6)
                stiff2mz = stiff2mz+nodalStiffnessTerms(i).val;
            else
                error('DOF not valid for  concentrated stiffness term.');
            end
            modNodalStiffnessTerms(i).val = 0.0;
        end
    end
    
    for i=1:length(nodalLoads)  %if node is specified in nodal forces terms file add to concentrated force for this node
         node1flagF=ismember(nodalLoads(i).nodeNum,node1);
         node2flagF=ismember(nodalLoads(i).nodeNum,node2);
        if(node1flagF==1)
            if(nodalLoads(i).dof==1)
                f1x = f1x+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==2)
                f1y = f1y+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==3)
                f1z = f1z+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==4)
                m1x = m1x+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==5)
                m1y = m1y+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==6)
                m1z = m1z+nodalLoads(i).val;
            else
                error('DOF not valid for  concentrated stiffness term.');
            end
            modNodalLoads(i).val = 0.0;
        end
        if(node2flagF==1)
            if(nodalLoads(i).dof==1)
                f2x = f2x+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==2)
                f2y = f2y+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==3)
                f2z = f2z+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==4)
                m2x = m2x+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==5)
                m2y = m2y+nodalLoads(i).val;
            elseif(nodalLoads(i).dof==6)
                m2z = m2z+nodalLoads(i).val;
            else
                error('DOF not valid for  concentrated stiffness term.');
            end
            modNodalLoads(i).val = 0.0;
        end
    end
    
    
    %compile nodal concentrated terms into mass, stiffness, and load arrays
    mass = [mass1 mass2;
            moix1 moix2
            moiy1 moiy2
            moiz1 moiz2];
        
    stiff = [stiff1x stiff2x;
             stiff1y stiff2y;
             stiff1z stiff2z;
             stiff1mx stiff2mx;
             stiff1my stiff2my;
             stiff1mz stiff2mz];
         
   load = [f1x f2x;
           f1y f2y;
           f1z f2z;
           m1x m2x;
           m1y m2y;
           m1z m2z];

end